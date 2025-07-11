/****
 * @date Created on 2025-07-07 at 10:03:22 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the SparseComplexMatrix methods. This contains the interface to UMFPACK.
 ***/
#include "SparseComplexMatrix.hpp"
#include <suitesparse/umfpack.h>

/**
 * Constructor of the SparseComplexMatrix object.
 */
SparseComplexMatrix::SparseComplexMatrix() {
    std::cout << TAG_INFO << "Creating SparseComplexMatrix (no allocation).\n";
    allocated = false; // At the beginning, the matrix is not allocated.
    mapped = false; // At the beginning, the matrix has not been mapped.
}

/**
 * Destructor of the SparseComplexMatrix object.
 */
SparseComplexMatrix::~SparseComplexMatrix() {
    if (allocated) {
        std::cout << TAG_INFO << "Deleting SparseComplexMatrix.\n";
        delete[] pcol;
        delete[] irow;
        delete[] map;
        delete[] real;
        delete[] imag;
        delete[] control;
        delete[] worki;
        delete[] work;
    }
}

/**
 * Returns the number of rows of the sparse matrix.
 */
int64_t SparseComplexMatrix::getNrow() const {
    if (not allocated) {
        throw std::logic_error("In getNrow(): Matrix is not allocated. There is no number of rows.");
    }
    return nrow;
}

/**
 * Returns the number of columns of the sparse matrix.
 */
int64_t SparseComplexMatrix::getNcol() const {
    if (not allocated) {
        throw std::logic_error("In getNcol(): Matrix is not allocated. There is no number of columns.");
    }
    return ncol;
}

/**
 * Returns the number of nonzero elements of the sparse matrix.
 */
int64_t SparseComplexMatrix::getNnz() const {
    if (not allocated) {
        throw std::logic_error("In getNnz(): Matrix is not allocated. There is no number of nonzero elements.");
    }
    return nnz;
}

/**
 * Setup the matrix elements from the given triplet form. This function assumes the input arrays (ti, tj, treal, timag) are
 * allocated with the correct size "nnz" stored in the present SparseComplexMatrix object.
 * This function assumes the triplets (ti, tj) are always the same, i.e., the order in which the elements are computed remains the same.
 * Furthermore, it assumes that there is no duplicate entries in the triplet (ti, tj, treal, timag).
 * In case of doubt, i.e., the triplet may have changed, then pass same=False. It is slower but safer.
 */
void SparseComplexMatrix::set(const int64_t nrow, const int64_t ncol, const int64_t nnz, 
                              const int64_t* ti, const int64_t* tj, const double* treal, const double* timag, const bool same) {
    
    if (not allocated) {// If the Jacobian is not allocated, then allocate memory.
        
        if (nrow <= 0 || ncol <= 0 || nnz <= 0) {
            throw std::invalid_argument("In SparseComplexMatrix(): Received negative size.");
        }
        else if (nnz > nrow*ncol) {
            throw std::invalid_argument("In SparseComplexMatrix(): Number of nonzeros is larger than what the matrix can contain.");
        }
        this->nrow = nrow;
        this->ncol = ncol;
        this->nnz = nnz;
        
        pcol = new int64_t[ncol+1];  // Indices of 'pointers' to columns.
        irow = new int64_t[nnz];     // Row indices of nonzero elements.
        map  = new int64_t[nnz];     // Map array relating the triplet format to the compressed column format (accelerates the construction of the matrix).
        real = new double[nnz];      // Real parts of nonzeor elements.
        imag = new double[nnz];      // Imaginary parts of nonzeor elements.
        worki = new int64_t[ncol];   // Workspace used by UMFPACK (move to class field).
        work = new double[4*ncol];   // Workspace array. Size=4*ncol without iterative refinement, and size=10*ncol with iterative refinement (for complex matrix type).
        
        control = new double[UMFPACK_CONTROL]; // UMFPACK's control parameters (input).
        umfpack_zl_defaults(control); // Setup the default UMFPACK parameters for complex arrays (and long indices).
        control[UMFPACK_PRL] = 1;     // Control UMFPACK's default printing level (see UMFPACK's Manual).
                                      // 0 or less: no output, even when an error occurs
                                      // 1: error messages only
                                      // 2 or more: print status, whether or not an error occurred
                                      // 4 or more: also print the UMFPACK Copyright
                                      // 6 or more: also print the UMFPACK License
        control[UMFPACK_IRSTEP] = 0;  // Set number of iterative refinement to zero. 
                                      // Since Newton-Raphson iteration solves many times the system, it is not necessary to be that accurate.
        allocated = true;             // Sets the allocation flag to true.
    }
    
    if (not mapped || not same) {// If the matrix has not been mapped yet, then call UMFPACK's triplet_to_col().
        
        int status = umfpack_zl_triplet_to_col(nrow, ncol, nnz, ti, tj, treal, timag, pcol, irow, real, imag, map);
        
        if (status != UMFPACK_OK) {// Check for possible UMFPACK errors.
            std::string info = "In setTriplet(): UMFPACK's triplet_to_col() failed. Error code = " + std::to_string(status) + ".";
            throw std::runtime_error(info);
        }
        mapped = true;  // After this, the matrix gets mapped, i.e., "map" is full with proper indices in this->real and this->imag.
    }
    else {// If the matrix is already mapped, then it is faster to fill it manually (see UMFPACK Manual).
        for (int k = 0; k < nnz; k++) {
            real[map[k]] = treal[k];
            imag[map[k]] = timag[k];
        }
    }
}

/**
 * Solve the linear system A.x = b using UMFPACK where the matrix A is the current SparseComplexMatrix object.
 */
void solveUmfpack(const SparseComplexMatrix& a, const ComplexVector& b, ComplexVector& x) {
    
    // Check for possible errors:
    if (b.size != a.nrow) {
        throw std::invalid_argument("In solveUmfpack(): Invalid size of 'b'.");
    }
    else if (x.size != a.ncol) {
        throw std::invalid_argument("In solveUmfpack(): Invalid size of 'x'.");
    }
    else if (not a.mapped) {
        throw std::logic_error("In solveUmfpack(): Matrix is not initialized. Please call set() to build it.");
    }
    
    // Prepare for calling UMFPACK:
    int sys = UMFPACK_A;  // Tells UMFPACK that the system to be solved has the form A.x = b.
    double info[UMFPACK_INFO];  // Output information (including UMFPACK's returned value "status").
    
    void *symbolic, *numeric;  // Symbolic and numeric factorization of the sparse matrix.
    
    umfpack_zl_symbolic(a.nrow, a.ncol, a.pcol, a.irow, a.real, a.imag, &symbolic, a.control, info);  // Perform symbolic reording to minimize fill-in.
    umfpack_zl_numeric(a.pcol, a.irow, a.real, a.imag, symbolic, &numeric, a.control, info); // Perform numericl LU factorization.
    umfpack_zl_free_symbolic(&symbolic); // Free the memory allocated by the symbolic factorization.
    umfpack_zl_wsolve(sys, a.pcol, a.irow, a.real, a.imag, x.real, x.imag, b.real, b.imag, numeric, a.control, info, a.worki, a.work);
    umfpack_zl_free_numeric(&numeric); // Free the memory allocated by the numeric factorization.
}

/**
 * Prints the SparseComplexMatrix object to standard output.
 * 
 * Argument:
 * 
 * verbosity: Printing level for UMFPACK (see UMFPACK Manual).
 *            2 or less: no output. returns silently without checking anything.
 *            3: fully check input, and print a short summary of its status
 *            4: as 3, but print first few entries of the input
 *            5: as 3, but print all of the input
 */
void SparseComplexMatrix::print(const int verbosity) const {
    if (not mapped) {
        throw std::logic_error("In print(): Matrix is not initialized. Nothing to print.");
    }
    double prl_tmp = control[UMFPACK_PRL]; // Save current UMFPACK's printing level.
    control[UMFPACK_PRL] = (double)verbosity;  // Set the printing level for UMFPACK (see UMFPACK Manual).
    int colform = 1; // Declares the matrix as column-oriented format (see UMFPACK Manual).
    
    umfpack_zl_report_matrix(nrow, ncol, pcol, irow, real, imag, colform, control);
    
    control[UMFPACK_PRL] = prl_tmp;  // Restore UMFPACK's default verbosity.
}

/**
 * Compare the present sparse matrix to a dense version.
 * This function should be used for testing purpose only.
 * Note that the "dense" array is modified. The output value is the difference dense - sparse.
 */
int SparseComplexMatrix::compareToDense(dcomplex* dense, const double tolerr) const {
    if (not mapped) {
        throw std::logic_error("In compareToDense(): Matrix is not initialized. Comparison is not possible.");
    }
    int fail = 0;  // Status flag: 0=Success, 1=Failure.
    int64_t i, j = 0;
    dcomplex sparse_elem, dense_elem;
    double error, maxerr = 0.;
    
    // 1. First check that all sparse element are located in the dense matrix:
    for (int k = 0; k < nnz; k++) {// Loop over the nonzero elements.
        
        // Find the (i,j) position of the k-th nonzero element:
        i = irow[k];
        while (pcol[j] <= k) j++;  // Compute the column index from the column pointer array.
        j--;
        
        // Compare with the dense matrix:
        sparse_elem = dcomplex(real[k], imag[k]);
        dense_elem  = dense[i + j*nrow];
        error = std::abs(sparse_elem - dense_elem)/std::max(std::abs(dense_elem), MEPS);
        dense[i + j*nrow] = 0.; // Remove the dense element so that we can check afterwhile that there is no missing element in sparse.
        
        if (error > maxerr) maxerr = error;  // Track the largest errors.
        
        if (error > tolerr) {// If the error exceeds the tolerance, then flag.
            fail = 1;
            std::cout << "[FAIL] At (i, j) = (" << i << ", " << j << "), "
                        << "sparse element = " << sparse_elem << ", " 
                        << "dense element = " << dense_elem << ", " 
                        << "error = " << error << " (> tol = " << tolerr << ").\n";
        }
    }
    
    // 2. Then check that there is no element in the dense matrix that is missing in sparse:
    for (i = 0; i < nrow; i++) {// Loop over all dense elements.
        for (j = 0; j < ncol; j++) {
            
            dense_elem = dense[i + j*nrow];  // At this step, every dense element is an error.
            error = std::abs(dense_elem);    // Compute the error.
            if (error > maxerr) maxerr = error;  // Track the largest errors.
            
            if (error > tolerr) {// If the error exceeds the tolerance, then flag.
                fail = 1;
                std::cout << "[FAIL] At (i, j) = (" << i << ", " << j << "), "
                            << "dense element = " << dense_elem << ", " 
                            << "error = " << error << " (> tol = " << tolerr << ").\n";
            }
        }
    }
    if (fail) {
        std::cout << "[FAIL] Detected errors in sparse Jacobian (largest errors = " << maxerr << ").\n";
    }
    else {
        std::cout << "[PASS] OK, no error found in sparse Jacobian (largest errors = " << maxerr << ").\n";
    }
    return fail;
}
