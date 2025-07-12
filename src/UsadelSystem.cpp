/****
 * @date Created on 2025-07-01 at 14:27:05 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the Usadel system methods which containing the algorithm for the Usadel equation.
 ***/
#include "UsadelSystem.hpp"
#include "BaseTools.hpp"
#include <cstring>
#include <random>
#include <chrono>
#include <iomanip>

/**
 * Constructor of the Usadel System object.
 */
UsadelSystem::UsadelSystem(SquareMesh& mesh, std::vector<Contact>& contact, const double holscat, const double holabso) {
    std::cout << TAG_INFO << "Creating UsadelSystem.\n";
    this->mesh = &mesh;
    this->contact = &contact;
    this->holscat = holscat;
    this->holabso = holabso;
    npoint = mesh.getNPoint();
    field = new ComplexVector(2*npoint);
}

/**
 * Destructor of the Usadel System object.
 */
UsadelSystem::~UsadelSystem() {
    std::cout << TAG_INFO << "Deleting UsadelSystem.\n";
    delete field;
}

/***********************************************************
 * GETTERS
 ***********************************************************/

/**
 * Returns the total number of points in the mesh.
 */
int UsadelSystem::getNPoint() const {
    return npoint;
}

/**
 * Returns the qvector at the point of index "ipoint".
 */
QVector UsadelSystem::getQVector(int ipoint) const {
    if (ipoint >= 0 && ipoint < npoint) {// If Q is inside the mesh, then simply convert the parameters to qvector.
        return QVector(field[0][2*ipoint], field[0][2*ipoint+1]);
    }
    else {
        throw std::out_of_range("In getQVector(): Index out of range.");
    }
}

/**
 * Returns the value of the transmission eigenvalue density according to the current solution of the Q field.
 * This function assumes the Q field given by the (theta, eta) parameters is solved.
 * 
 * @todo This function remains to be implemented...........
 */
double UsadelSystem::getRho() const {
    
    //TODO: To be implemented......
    
    return 0.;
}

/***********************************************************
 * INTERNAL COMPUTATIONAL METHODS
 ***********************************************************/

/**
 * Returns the Q(x) field in vector form at a given point "ipoint" in the mesh, assuming a given (theta, eta) field "someField".
 * This function takes into account the boundary conditions for negative values of "i".
 * If "i" is BND_MIRROR, then returns 0, and if "i" is BND_OPEN, then returns the vacuum Q field.
 */
QVector UsadelSystem::qvectorAtPos(const ComplexVector& someField, int ipoint) const {
    if (ipoint >= 0 && ipoint < npoint) {// If Q is inside the mesh, then simply convert the parameters to qvector.
        return QVector(someField[2*ipoint], someField[2*ipoint+1]);
    }
    else if (ipoint == BND_MIRROR) {// Mirror boundary conditions. If Q is outside the mesh, Q is simply zero.
        return QVector(0., 0., 0.);
    }
    else if (ipoint == BND_OPEN) {// Open boundary conditions. If Q is outside the mesh, then Q is replaced by its vacuum value.
        return QVector(0., 0., holscat/EXTRAPOLEN);
    }
    else {// If non of the above, then throw an error.
        throw std::out_of_range("In getQVector(): Index out of range.");
    }
}

/**
 * Returns the two residuals of the Usadel equation at a given point of index "i" in the mesh and assuming a given (theta, eta) field "someField".
 * This method gives nothing but a numerical definition of the Usadel equation, and thus lies at the core of the program.
 */
void UsadelSystem::residualAtPos(const ComplexVector& someField, int ipoint, dcomplex& res1, dcomplex& res2) const {
    
    // Get the Q field (vector form) of the 4-point neighborhood (boundary conditions are included):
    QVector qvn, qvs, qve, qvw, l3v, lav, sqv;
    MeshPoint p = mesh->getPoint(ipoint);  // If index "ipoint" is not in the mesh, then this should throw an error.
    qvn = qvectorAtPos(someField, p.north);
    qvs = qvectorAtPos(someField, p.south);
    qve = qvectorAtPos(someField, p.east);
    qvw = qvectorAtPos(someField, p.west);
    lav = QVector(0., 0., DIMENSION*holscat*holabso);  // Independent term of the Usadel equation.
    
    // Deal with contact interactions:
    Vector2D pv(p); // Get the central point in vector form.
    
    for (unsigned int ict = 0; ict < contact->size(); ict++) {// Loop over all contact interactions.
        contact->at(ict).applyTransfo(pv, pv + Vector2D(0., +1.), qvn);
        contact->at(ict).applyTransfo(pv, pv + Vector2D(0., -1.), qvs);
        contact->at(ict).applyTransfo(pv, pv + Vector2D(+1., 0.), qve);
        contact->at(ict).applyTransfo(pv, pv + Vector2D(-1., 0.), qvw);
    }
    
    // Compute the full residual q-vector of the Usadel equation:
    sqv = (qvn + qvs + qve + qvw + lav).rotate_2(-someField[2*ipoint+1]).rotate_1(-someField[2*ipoint]);
    
    res1 = -sqv.q2;  // The residual is given by 1_z.cross(sqv), resulting in this expression.
    res2 =  sqv.q1;
}

/**
 * Returns the two residuals of the Usadel equation at a given point of index "i" in the mesh and assuming a given (theta, eta) field "someField".
 * 
 * @deprecated Old version of the residual. Should only be used for testing purposes.
 */
void UsadelSystem::residualAtPosOld(const ComplexVector& someField, int ipoint, dcomplex& res1, dcomplex& res2) const {
    
    // Get the Q field (vector form) at the five points of the stencil (boundary conditions are included):
    QVector qv0, qvn, qvs, qve, qvw, l3v, resv;
    MeshPoint p0 = mesh->getPoint(ipoint);  // If index "ipoint" is not in the mesh, then this should throw an error.
    qv0 = qvectorAtPos(someField, ipoint);
    qvn = qvectorAtPos(someField, p0.north);
    qvs = qvectorAtPos(someField, p0.south);
    qve = qvectorAtPos(someField, p0.east);
    qvw = qvectorAtPos(someField, p0.west);
    l3v = QVector(0., 0., DIMENSION*holscat*holabso);
    
    // Deal with contact interactions:
    Vector2D p0v(p0); // Get the central point in vector form.
    
    for (unsigned int ict = 0; ict < contact->size(); ict++) {// Loop over all contact interactions.
        contact->at(ict).applyTransfo(p0v, p0v + Vector2D(0., +1.), qvn);
        contact->at(ict).applyTransfo(p0v, p0v + Vector2D(0., -1.), qvs);
        contact->at(ict).applyTransfo(p0v, p0v + Vector2D(+1., 0.), qve);
        contact->at(ict).applyTransfo(p0v, p0v + Vector2D(-1., 0.), qvw);
    }
    
    // Compute the full residual q-vector of the Usadel equation:
    resv = qv0.cross(qvn + qvs + qve + qvw + l3v).rotate_2(-someField[2*ipoint+1]).rotate_1(-someField[2*ipoint]);
    // The final rotations R_y(-eta).R_x(-theta).resv are intended to cancel the third component of resv exactly.
    // Therefore, in principle, resv.q3 should be exactly zero (within numerical accuracy):
    //std::cout << TAG_INFO << "Residual[i=" << i << "] = " << resv << "\n";
    
    // Compute the two components that will be used for residual minimization (note the rotation):
    res1 = resv.q1;
    res2 = resv.q2;
}

/**
 * Computes the whole residual vector that should be minimized using the Newton solver.
 * This method assumes "resvector" is already allocated to the size 2*npoint.
 */
void UsadelSystem::residualVector(const ComplexVector& someField, ComplexVector& resvec) const {
    dcomplex res1, res2;
    for (int ipoint = 0; ipoint < npoint; ipoint++) {
        residualAtPos(someField, ipoint, res1, res2);
        resvec.set(2*ipoint,   res1);
        resvec.set(2*ipoint+1, res2);
    }
}

/**
 * Compute the Jacobian matrix element corresponding to the derivative of the residual at the point of index "ipoint" with respect to the angle parameter of index "jparam".
 * This function applies central finite difference method using the given perturbed fields "fieldPlus" and "fieldMinus".
 * The resulting derivatives are stored in the sparse-matrix arrays: irow, icol, jacreal, jacimag.
 * This function assumes the arrays irow, icol, jacreal, jacimag are allocated enough space: 20*npoint (upper bound).
 */
void UsadelSystem::computeJacobianAtPos(const ComplexVector& fieldPlus, const ComplexVector& fieldMinus, int64_t ipoint, int64_t jparam, 
                                        int64_t& nnz, int64_t* irow, int64_t* icol, double* jacreal, double* jacimag) const {
    
    if (ipoint >= 0 && ipoint < npoint) {// Only proceeds if the point is in the mesh.
        dcomplex res1p, res2p, res1m, res2m, der1, der2, htrue;
        
        // First compute the residuals:
        residualAtPos(fieldPlus,  ipoint, res1p, res2p);
        residualAtPos(fieldMinus, ipoint, res1m, res2m);
        
        // Compute the derivatives of the two residuals "res1" and "res2":
        htrue = fieldPlus[jparam] - fieldMinus[jparam];  // Compute the true finite step.
        der1 = (res1p - res1m)/htrue;
        der2 = (res2p - res2m)/htrue;
        
        // Add the derivative of "res1" with respect to parameter of index "jparam" to the sparse Jacobian:
        icol[nnz] = jparam;
        irow[nnz] = 2*ipoint;
        jacreal[nnz] = der1.real();
        jacimag[nnz] = der1.imag();
        nnz++;
        
        // Add the derivative of "res2" with respect to parameter of index "jparam" to the sparse Jacobian:
        icol[nnz] = jparam;
        irow[nnz] = 2*ipoint+1;
        jacreal[nnz] = der2.real();
        jacimag[nnz] = der2.imag();
        nnz++;
    }
}

/**
 * Compute the sparse Jacobian matrix of the residual vector corresponding to the Usadel equation.
 * This function uses central difference formula to approximate the Jacobian.
 * This function assumes the (theta, eta) parameters are given by "someField".
 */
void UsadelSystem::computeJacobian(const ComplexVector& someField, SparseComplexMatrix& jac) const {
    
    int64_t jpoint, nparam = 2*npoint;
    double h;  // Discrete step used in the finite difference formula.
    MeshPoint p;
    
    // Allocate space the sparse Jacobian (triplet format):
    int64_t nnz = 0;  // Initialize the number of nonzero element, used as an index in "irow", "icol", "jacreal", "jacimag".
    int64_t nnzub = 20*npoint;  // Expected number of nonzero elements (safe upper bound).
    auto irow = new int64_t[nnzub];
    auto icol = new int64_t[nnzub];
    auto jacreal = new double[nnzub];
    auto jacimag = new double[nnzub];
    
    // Performs a deep copy of the field using copy constructor:
    ComplexVector fieldPlus  = someField;
    ComplexVector fieldMinus = someField;
    
    for (int64_t jparam = 0; jparam < nparam; jparam++) {// Loop over the columns of the jacobian (all parameters to be varied).
        
        // Compute the array of parameters:
        h = CBRTEPS*std::max(std::abs(someField[jparam]), 1.);  // Compute the discrete step (a real value).
        fieldPlus .add(jparam, +h);  // Increment the parameter at index "jparam".
        fieldMinus.add(jparam, -h);  // Increment the parameter at index "jparam".
        
        if (jparam != 0) {// If this is not the first variation, then cancel the previous variation.
            fieldPlus .set(jparam-1, someField[jparam-1]);
            fieldMinus.set(jparam-1, someField[jparam-1]);
        }
        
        // Look for neighboring points of "j":
        jpoint = jparam/2; // Index of the point whose parameter is varied.
        p = mesh->getPoint(jpoint); // Get the neighborhood of the point which is varied.
        computeJacobianAtPos(fieldPlus, fieldMinus, jpoint,  jparam, nnz, irow, icol, jacreal, jacimag);
        computeJacobianAtPos(fieldPlus, fieldMinus, p.north, jparam, nnz, irow, icol, jacreal, jacimag);
        computeJacobianAtPos(fieldPlus, fieldMinus, p.south, jparam, nnz, irow, icol, jacreal, jacimag);
        computeJacobianAtPos(fieldPlus, fieldMinus, p.east,  jparam, nnz, irow, icol, jacreal, jacimag);
        computeJacobianAtPos(fieldPlus, fieldMinus, p.west,  jparam, nnz, irow, icol, jacreal, jacimag);
        
    }
    
    //jacreal[nnzub/3] += 0.001;  // Introduces a malicious error for tesing.
    
    jac.set(nparam, nparam, nnz, irow, icol, jacreal, jacimag, true); // Call the initializer of the sparse Jacobian object.
                                                                      // TODO: Maybe remove the "same" parameter by storing "ti", "tj" in the sparse object ??
    
    // Free memory held for the triplet form:
    delete[] irow;
    delete[] icol;
    delete[] jacreal;
    delete[] jacimag;
    
}

/**
 * Compute the dense Jacobian matrix of the residual vector corresponding to the Usadel equation.
 * This version uses central difference formula to approximate the Jacobian.
 * This function assumes the array "jacdense" is allocated enough space: (2*npoint)^2.
 * This matrix is stored in column-major format.
 * 
 * @deprecated This function should only be used for test comparison with computeJacobian() because it is far less efficient and consumes much memory.
 */
void UsadelSystem::computeJacobianDense(const ComplexVector& someField, dcomplex* jacdense) const {
    
    double h;  // Discrete step used in the finite difference formula.
    int nrow = 2*npoint, ncol = 2*npoint;  // Number of rows/columns in the Jacobian.
    dcomplex res1p, res2p, res1m, res2m, der1, der2, htrue;
    
    ComplexVector fieldPlus  = field[0];  // Performs a deep copy of the field using copy constructor.
    ComplexVector fieldMinus = field[0];
    
    for (int jparam = 0; jparam < ncol; jparam++) {// Loop over the columns of the jacobian (all parameters to be varied).
        
        // Compute the array of parameters:
        h = CBRTEPS*std::max(std::abs(someField[jparam]), 1.);  // Compute the discrete step (a real value).
        fieldPlus .add(jparam, +h);  // Increment the parameter at index "jparam".
        fieldMinus.add(jparam, -h);  // Increment the parameter at index "jparam".
        
        if (jparam != 0) {// If this is not the first variation, then cancel the previous variation.
            fieldPlus .set(jparam-1, field[0][jparam-1]);
            fieldMinus.set(jparam-1, field[0][jparam-1]);
        }
        
        for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh to compute the residuals.
            
            residualAtPos(fieldPlus,  ipoint, res1p, res2p);
            residualAtPos(fieldMinus, ipoint, res1m, res2m);
            
            htrue = fieldPlus[jparam] - fieldMinus[jparam];  // Compute the true finite step.
            jacdense[2*ipoint + jparam*nrow]     = (res1p - res1m)/htrue;
            jacdense[2*ipoint + 1 + jparam*nrow] = (res2p - res2m)/htrue;
            
        }
    }
}

/***********************************************************
 * TESTING METHODS:
 ***********************************************************/

/**
 * Test the computation of the residual by comparison with an older version.
 * This function assumes the "field" is already initialized.
 * Return 0 on success, 1 on failure.
 */
int UsadelSystem::testResidual(const double tolerr) const {
    
    int fail = 0; // Status flag (0=Success, 1=Failure).
    dcomplex res1, res2, res1_expc, res2_expc;
    double error_res1, error_res2, maxerr1 = 0., maxerr2 = 0.;
    
    for (int i = 0; i < npoint; i++) {// Loop on the points of the mesh.
        
        residualAtPos(field[0], i, res1, res2);  // Compute the two residuals on each point with the best method.
        residualAtPosOld(field[0], i, res1_expc, res2_expc);  // Compute the two residuals on each point with the best method.
        error_res1 = std::abs(res1 - res1_expc)/std::max(std::abs(res1_expc), MEPS);  // Relative error for res1.
        error_res2 = std::abs(res2 - res2_expc)/std::max(std::abs(res2_expc), MEPS);  // Relative error for res2.
        
        if (error_res1 > maxerr1) maxerr1 = error_res1;  // Track the largest errors.
        if (error_res2 > maxerr2) maxerr2 = error_res2;
        
        if (error_res1 > tolerr || error_res2 > tolerr) {// If the error exceeds the tolerance, then flag.
            fail = 1;
            std::cout << "[FAIL] At pos i=" << i << ", found residuals [res1, res2] = [" << res1 << ", " << res2 << "], " 
                      << "expected = [" << res1_expc << ", " << res2_expc << "], errors = [" << error_res1 << ", " << error_res2 << "].\n";
        }
    }
    if (fail) {
        std::cout << "[FAIL] Detected errors in residuals (largest errors = [" << maxerr1 << ", " << maxerr2 << "]).\n";
    }
    else {
        std::cout << "[PASS] OK, no error found in residuals (largest errors = [" << maxerr1 << ", " << maxerr2 << "]).\n";
    }
    return fail;
}

/**
 * Test the sparse Jacobian by comparison with a safer (but far less efficient) dense Jacobian.
 * This function assumes the "field" is already initialized.
 * @warning This function should only be called on small meshes to avoid memory overflow.
 */
int UsadelSystem::testJacobian(const double tolerr) const {
    
    // 1. Compute the sparse Jacobian:
    SparseComplexMatrix jac;
    int nnz_expc = 20*npoint;  // Expected number of nonzero elements (safe upper bound).
    int nparam = 2*npoint;     // Expected size of the Jacobian (number of rows, or number of coolumns).
    
    std::cout << TAG_INFO << "Computing sparse Jacobian... ";
    auto start = std::chrono::steady_clock::now();
    computeJacobian(field[0], jac);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Done in " << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
    std::cout << TAG_INFO << "Found nnz = " << jac.getNnz() << ", expected nnz = " << nnz_expc << " (20*npoint with npoint = " << npoint << ").\n";
    
    // 2. Compute the reference dense Jacobian:
    auto jacdense = new dcomplex[nparam*nparam]();  // Initialize dense Jacobian (column-major ordering) with zero intialization.
    
    std::cout << TAG_INFO << "Computing dense reference Jacobian... ";
    start = std::chrono::steady_clock::now();
    computeJacobianDense(field[0], jacdense);
    end = std::chrono::steady_clock::now();
    std::cout << "Done in " << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
    
    //jacdense[99 + 101*nparam] += 0.001; // Introduces a malicious error for testing.
    
    int fail = jac.compareToDense(jacdense, tolerr);
    
    delete[] jacdense;
    
    // 3. Comparison test:
    return fail;
}


/***********************************************************
 * PUBLIC COMPUTATIONAL METHODS:
 ***********************************************************/

/**
 * Initialize the Q field to random values according to a given "seed".
 * This method ensures the Q values are picked up isotropically on the Q^2 = 1 manifold.
 */
void UsadelSystem::initRandom(const uint64_t seed) {
    double q1r, q1i, q2r, q2i, q3r, q3i;
    dcomplex theta, eta;
    QVector qv;
    
    std::mt19937_64 rng;  // Instantiate the standard Mersenne Twister random number generator (64-bit version).
    rng.seed(seed);       // Initialize the random generator with the given seed.
    std::normal_distribution<double> random_normal(0., 1.);
    
    for (int i = 0; i < npoint; i++) {
        
        q1r = random_normal(rng);  // This separation forces the random generator to be called in this specific order.
        q1i = random_normal(rng);
        q2r = random_normal(rng);
        q2i = random_normal(rng);
        q3r = random_normal(rng);
        q3i = random_normal(rng);
        
        qv = QVector(dcomplex(q1r, q1i), dcomplex(q2r, q2i), dcomplex(q3r, q3i));
        //qv.getThetaEta(theta, eta);  // Convert the qvector into angle parameters (theta, eta).
        field[0].set(2*i,   qv.getTheta());
        field[0].set(2*i+1, qv.getEta());
    }
}

/**
 * Initialize the Q field assuming theta = 0 and eta = arctan(-i*sqrt(gamma_a*gamma_b)).
 * This method should be called before running the Newton solver.
 */
void UsadelSystem::initConstant() {
    
    dcomplex theta0, eta0, ga_avg;
    
    // 1. Compute the geometric average of all gamma parameters (this is coarse):
    ga_avg = dcomplex(1., 0.);
    unsigned int nct = contact->size();
    for (unsigned int ict = 0; ict < nct; ict++) {// Loop over the contact interactions.
        ga_avg *= contact->at(ict).getGamma();
    }
    ga_avg = pow(ga_avg, 1./nct);  // Geometric average (in the case gamma_a != gamma_b).
    
    // 2. Assign constant values to all points in the mesh:
    theta0 = dcomplex(0., 0.); // Very coarse. Theta always varies in space.
    eta0 = std::atan(-I*ga_avg);    // Exact only in the absence of losses (no absorption, complete channel control).
    
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        field[0].set(2*ipoint,   theta0);
        field[0].set(2*ipoint+1, eta0);
    }
}

/**
 * Solve the Usadel System using the Newton-Raphson iterative algorithm.
 * This method does not call initializers (parameters (theta, eta) are zeros), but the usage of an appropriate ansatz
 * is highly recommended in order to accelerate convergence while avoiding improper solutions.
 * 
 * Arguments:
 * 
 * maxit   = Maximum number of iterations. Typically: 100-500.
 * nsub    = Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
 *           The smallest reduction factor is thus 2^-nsub and should be larger than the machine epsilon (2^-53).
 * tolp    = Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
 * tolr    = Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
 * verbose = Verbosity level in standard output. 0=No output, 1=Display each iteration.
 * 
 * Returns:
 * 
 * found   = Flag equal to 0 if the Newton-Raphson algorithm did not convergence within the prescribed number of iterations,
 *           and 1 if a solution has been found, i.e., the convergence criteria have been met.
 */
int UsadelSystem::solveNewton(const int maxit, const int nsub, const double tolp, const double tolr, const int verbose) {
    
    int found = 0;  // Return status. 0=No convergence (solution not found), 1=Success (found solution).
    int iter, s;    // Iteration indices.
    double fac, resnorm, resnorm0, newresnorm, deltanorm, fieldnorm;
    
    const int nparam = 2*npoint;          // Size of the vectors and matrix (number of columns, or number of rows).
    ComplexVector res(nparam), newfield(nparam), newres(nparam), delta(nparam);  // Allocate vectors: residual, new field, new residual, displacement.
    SparseComplexMatrix jac;  // Declare sparse Jacobian (no allocation).
    
    residualVector(field[0], res); // Compute the residual from the current field.
    
    resnorm = res.norm(); // Compute the norm of the residual.
    resnorm0 = resnorm;   // Store the norm of the first residual for the stopping criterion.
    
    for (iter = 0; iter < maxit; iter++) {// Loop over the allowed Newton-Raphson iterations.
        
        // 1. Solves the Jacobian system using UMFPACK to find the Newton-Raphson displacement:
        computeJacobian(field[0], jac);  // Compute the sparse Jacobian using finite differences. This allocates memory only at the first iteration.
        solveUmfpack(jac, res, delta);  // Solve the system Jac.delta = resv, for the Newton-Raphson displacement "delta".
        
        // 2. Backtracking line search along the direction to minimize the residual (ensure convergence of Newton-Raphson):
        fac = -1.;  // Negative because the actual Newton-Raphson displacement is: newfield = field[0] - delta.
        for (s = 0; s < nsub; s++) {
            newfield = field[0] + fac*delta;
            residualVector(newfield, newres); // Compute the new residual corresponding to "newfield".
            newresnorm = newres.norm();
            if (newresnorm < resnorm) break;
            fac /= 2.;  // Reduces the line search parameter.
        }
        
        // 3. Saves the solution with minimal residual (and print the status):
        field[0] = newfield;   // The new field becomes the current one (deep copy).
        res = newres;          // The new residual becomes the current one (deep copy).
        resnorm = newresnorm;  // Update the residual norm with the new one.
        deltanorm = delta.norm();
        fieldnorm = field[0].norm();
        
        if (verbose >= 1) {// Print the current status.
            std::cout << TAG_INFO << "#" << iter << " | fieldn = " << fieldnorm << ", deltan = " << deltanorm 
                      << ", resn = " << resnorm << ", sub = " << s+1 << " (fac=" << fac << ")\n";
        }
        
        // 4. Stopping criterion:
        if (deltanorm < tolp*std::max(fieldnorm, 1.) && resnorm < tolr*std::max(resnorm0, 1.)) {
            if (verbose >= 1) {
                std::cout << TAG_INFO << "Found solution\n";
            }
            found = 1;
            break;
        }
    }
    
    return found;
}

/**
 * Save the present field to a file of given "filename" using the string "sep" as a separator.
 * Columns are: x, y, theta_re, theta_im, eta_re, eta_im, q11_re, q11_im, q12_re, q12_im, q21_re, q21_im
 */
void UsadelSystem::saveField(const char* filename, const char* sep, const int prec) const {
    
    MeshPoint p;
    dcomplex theta, eta;
    QVector qv;
    
    std::ofstream ofs;  // Declare output stream object.
    ofs.open(filename); // Open the file in write mode.
    
    const auto default_precision = ofs.precision(); // Saves the default precision.
    ofs << std::setprecision(prec); // Set the printing precision.
    
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    ofs << "%% Parameters: holscat = " << holscat << ", holabso = " << holabso << "\n";
    ofs << "x, y, theta_re, theta_im, eta_re, eta_im, q11_re, q11_im, q12_re, q12_im, q21_re, q21_im\n";
    
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        
        p = mesh->getPoint(ipoint);
        theta = field[0][2*ipoint]; // Extract the angle parameters.
        eta   = field[0][2*ipoint + 1];
        qv = QVector(theta, eta);
        ofs << p.x << sep << p.y << sep
            << theta.real() << sep << theta.imag() << sep 
            << eta.real() << sep << eta.imag() << sep 
            << qv.getQ11().real() << sep << qv.getQ11().imag() << sep
            << qv.getQ12().real() << sep << qv.getQ12().imag() << sep
            << qv.getQ21().real() << sep << qv.getQ21().imag() << "\n";
    }
    
    ofs.close();  // Close the file.
    
    ofs << std::setprecision(default_precision); // Restore default precision for printing.
}
