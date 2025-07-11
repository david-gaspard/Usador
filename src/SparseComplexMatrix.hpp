/****
 * @date Created on 2025-07-07 at 09:46:51 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SparseComplexMatrix object, an interface to UMFPACK.
 * See also UMFPACK Manual: https://fossies.org/linux/SuiteSparse/UMFPACK/Doc/UMFPACK_UserGuide.pdf.
 ***/
#ifndef _SPARSECOMPLEXMATRIX_H
#define _SPARSECOMPLEXMATRIX_H
#include "ComplexVector.hpp"

/**
 * Class defining the SparseComplexMatrix object based on the compressed column format (CSC format) used by UMFPACK.
 */
class SparseComplexMatrix {
    
    private:
    
    int64_t nrow;   // Number of rows.
    int64_t ncol;   // Number of columns.
    int64_t nnz;    // Number of nonzero elements.
    
    bool allocated; // Flag indicating if the arrays of the sparse matrix have been allocated.
    bool mapped;    // Flag indicating if the matrix has already been mapped. If true, then use "map" to accelerate construction. 
                    // Note that this assumes the matrix elements do not change position in the matrix.
    
    int64_t* pcol;  // Array of 'pointers' to the columns of the matrix (see UMFPACK Manual). Size = ncol+1.
                    // Column j of the matrix A is held in Ai[ Ap[j], Ap[j]+1, ..., Ap[j+1]-1 ].
                    // The first entry, Ap[0], must be zero, and the last one, nz=Ap[ncol] is thus the total number of entries in the pattern of the matrix A.
    int64_t* irow;  // Array of row indices of matrix elements. Size = nnz.
    int64_t* map;   // Array holding the position of the triplets in the column-form matrix, that is Ax[Map[k]] = Tx[k] (see UMFPACK Manual). Size = nnz.
                    // This array is used to accelerate the construction of the matrix with setTriplet(). 
                    // Note that this assumes the matrix elements do not change position in the matrix.
    
    double* real;   // Array of real parts of matrix elements. Size = nnz.
    double* imag;   // Array of imaginary part of matrix elements. Size = nnz.
    
    double* control;  // Control array used by UMFPACK. Size = UMFPACK_CONTROL. 
    int64_t* worki;   // Workspace array used by UMFPACK's *_wsolve(). See UMFPACK Manual. Size = ncol.
    double* work;     // Workspace array used by UMFPACK's *_wsolve(). See UMFPACK Manual. This avoids UMFPACK allocating memory many times.
                      // Size = 4*ncol without iterative refinement. Size = 10*ncol with iterative refinement.
    
    public:
    
    // Constructors/Destructor:
    SparseComplexMatrix();
    ~SparseComplexMatrix();
    
    // Getters:
    int64_t getNrow() const;
    int64_t getNcol() const;
    int64_t getNnz() const;
    
    // Methods:
    void set(const int64_t nrow, const int64_t ncol, const int64_t nnz, const int64_t* ti, 
             const int64_t* tj, const double* treal, const double* timag, const bool same);
    void print(const int verbosity) const;
    int compareToDense(dcomplex* dense, const double tolerr) const;
    friend void solveUmfpack(const SparseComplexMatrix& a, const ComplexVector& b, ComplexVector& x);  // Promote access to private matrix field to the linear solver.
    
};

#endif
