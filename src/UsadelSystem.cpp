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
 * 
 * Arguments:
 * 
 * name    = Name of the UsadelSystem (typically describing the system geometry) used for generating output files.
 * mesh    = Given SquareMesh object, properly initialized (see: SquareMesh).
 * holscat = Scattering strength, h/lscat, where "h" is the lattice step, and "lscat" is the mean free path.
 * holabso = Absorption strength, h/labso, where "h" is the lattice step, and "labso" is the ballistic absorption length.
 * tval    = Transmission eigenvalue, between 0 and 1.
 */
UsadelSystem::UsadelSystem(const std::string& name, SquareMesh& mesh, const double holscat, const double holabso, const double tval) {
    //std::cout << TAG_INFO << "Creating UsadelSystem.\n";
    this->mesh = new SquareMesh(mesh);  // Call copy constructor (implicitly defined).
    npoint = mesh.getNPoint();
    field = new ComplexVector(2*npoint);
    this->holscat = holscat;
    this->holabso = holabso;
    setTransmission(tval);
    this->name = name;
}

/**
 * Constructor creating a disordered rectangular waveguide of given length "length", width "width", scattering thickness "dscat", and absorption thickness "dabso".
 * This constructor is mainly used for tests because the waveguide solution is accurately known for dabso=0.
 * 
 * Arguments:
 * 
 * name   = Name of the UsadelSystem (typically describing the system geometry) used for generating output files.
 * length = Horizontal length of the waveguide in units of the lattice step.
 * width  = Vertical width of the waveguide in units of the lattice step.
 * dscat  = Scattering thickness, L/lscat, where "L" is the length of the waveguide and "lscat" is the scattering mean free path.
 * dabso  = Absorption thickness, L/labso, where "L" is the length of the waveguide and "labso" is the ballistic absorption length.
 * tval   = Transmission eigenvalue, between 0 and 1.
 */
UsadelSystem::UsadelSystem(const std::string& name, const int length, const int width, const double dscat, const double dabso, const double tval) {
    //std::cout << TAG_INFO << "Creating UsadelSystem (waveguide template).\n";
    mesh = new SquareMesh();
    mesh->addRectangle(0, length, 0, width, BND_MIRROR);
    mesh->setBoundaryRectangle(0, 0, 0, width, DIR_WEST, BND_INPUT);
    mesh->setBoundaryRectangle(length, length, 0, width, DIR_EAST, BND_OUTPUT);
    mesh->finalize();
    
    npoint = mesh->getNPoint();
    field = new ComplexVector(2*npoint);
    holscat = dscat/length;
    holabso = dabso/length;
    setTransmission(tval);
    this->name = name;
}

/**
 * Define explicit copy constructor for the UsadelSystem object (because of pointers). Used especially in parallel computing.
 */
UsadelSystem::UsadelSystem(const UsadelSystem& usys) {
    //std::cout << TAG_INFO << "Creating UsadelSystem (copy constructor).\n";
    this->mesh = new SquareMesh(*usys.mesh);  // Call copy constructor (implicitly defined).
    npoint = usys.mesh->getNPoint();
    field = new ComplexVector(2*npoint);
    this->name = usys.name;
    copy(usys); // Copy the field contained in the given UsadelSystem to the present UsadelSystem.
}

/**
 * Destructor of the Usadel System object.
 */
UsadelSystem::~UsadelSystem() {
    //std::cout << TAG_INFO << "Deleting UsadelSystem.\n";
    delete mesh;
    delete field;
}

/**
 * Assigns the transmission value to the present Contact interaction object.
 */
void UsadelSystem::setTransmission(const double tval) {
    if (tval <= 0. || tval > 1.) {
        throw std::invalid_argument("In setTransmission(): Invalid transmission value, expected in (0, 1].");
    }
    this->tval = tval;
    this->ga = dcomplex(std::sqrt(1./tval), SQRTEPS);
    this->gb = this->ga;
}

/**
 * Copy the solution of the Usadel equation stored in the given UsadelSystem to the present UsadelSystem.
 */
void UsadelSystem::copy(const UsadelSystem& usys) {
    if (this->name != usys.name) {
        throw std::invalid_argument("In copy(): Copying two UsadelSystem with different names is forbidden.");
    }
    setTransmission(usys.tval);
    this->holscat = usys.holscat;
    this->holabso = usys.holabso;
    field[0] = usys.field[0];  // Deep copy the field.
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
 * Returns the i-th point of the mesh.
 */
MeshPoint UsadelSystem::getPoint(const int ipoint) const {
    return mesh->getPoint(ipoint);
}

/**
 * Returns the index of the given point (x, y).
 * If (x, y) does not belong the mesh, then return BND_DEFAULT, a negative integer.
 */
int UsadelSystem::indexOf(const int x, const int y) const {
    return mesh->indexOf(x, y);
}

/**
 * Returns the value of "holscat" which is defined by h/lscat, where "h" is the lattice step
 * (the unit length) and "lscat" is the scattering mean free path.
 */
double UsadelSystem::getHolscat() const {
    return holscat;
}

/**
 * Returns the value of "holabso" which is defined by h/labso, where "h" is the lattice step
 * (the unit length) and "labso" is the ballistic absorption length.
 */
double UsadelSystem::getHolabso() const {
    return holabso;
}

/**
 * Returns the transmission eigenvalue T between 0 and 1.
 */
double UsadelSystem::getTransmission() const {
    return tval;
}

/**
 * Returns a deep copy of the name of the UsadelSystem.
 */
std::string UsadelSystem::getName() const {
    return name;
}

/**
 * Returns the value of the Q field (in vector representation) at the point of index "ipoint".
 */
QVector UsadelSystem::getQVector(const int ipoint) const {
    return qvectorAtPos(field[0], ipoint);
}

/**
 * Returns the value of the matrix current J defined by J = -(lscat/d) * Q * grad Q.
 * between the point of index "ipoint" and the point of index "jpoint", in the direction ipoint -> jpoint.
 * If "ipoint" or "jpoint" is equal to one of the boundary values BND_*, then uses the Q field from the boundary condition.
 * This function assumes the Q field given by the (theta, eta) parameters is solved.
 */
QVector UsadelSystem::getJVector(const int ipoint, const int jpoint) const {
    return qvectorAtPos(field[0], jpoint).cross(qvectorAtPos(field[0], ipoint)) * (I/(DIMENSION * holscat));
}

/**
 * Returns the value of the transmission eigenvalue density according to the current solution of the Q field.
 * This function assumes the Q field given by the (theta, eta) parameters is solved.
 */
double UsadelSystem::getRho() const {
    
    MeshPoint p;
    QVector lpv(dcomplex(1., 0.), dcomplex(0., +1.), dcomplex(0., 0.));
    QVector lmv(dcomplex(1., 0.), dcomplex(0., -1.), dcomplex(0., 0.));
    dcomplex flux, fun;
    int len_input, len_output, ninput, noutput;
    
    flux = dcomplex(0., 0.);
    len_input = 0;
    len_output = 0;
    
    for (int i = 0; i < npoint; i++) {// Loop on all the contact interactions.
        
        p = mesh->getPoint(i);
        
        ninput  = (p.north == BND_INPUT)  + (p.south == BND_INPUT)  + (p.east == BND_INPUT)  + (p.west == BND_INPUT);
        noutput = (p.north == BND_OUTPUT) + (p.south == BND_OUTPUT) + (p.east == BND_OUTPUT) + (p.west == BND_OUTPUT);
        
        if (ninput != 0) {
            flux += dcomplex(ninput, 0.) * getJVector(BND_INPUT, i).dot(lpv);  // Note that the input current must be computed from the boundary to the edge point (and not in the other direction).
            len_input += ninput;
        }
        if (noutput != 0) {
            flux += dcomplex(noutput, 0.) * getJVector(i, BND_OUTPUT).dot(lmv);  // Note that the output current must be computed from the edge to the boundary.
            len_output += noutput;
        }
    }
    
    //std::cout << TAG_INFO << "In getRho(): len_input = " << len_input << ", len_output = " << len_output << "\n";
    //std::cout << TAG_INFO << "In getRho(): Total flux = " << flux << "\n";
    
    fun = DEXTPOL * (I*std::sqrt(tval)/(2.*std::min(len_input, len_output))) * flux;
    
    return fun.imag()/(PI*tval*tval);
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
        return QVector(0., 0., 1.) * (holscat/EXTRAPOLEN);
    }
    else if (ipoint == BND_INPUT) {// Input boundary condition. If Q is outside the mesh, then Q is replaced by Q_a.
        return QVector(-I*ga, +ga, 1.) * (holscat/EXTRAPOLEN);
    }
    else if (ipoint == BND_OUTPUT) {// Output boundary condition. If Q is outside the mesh, then Q is replaced by Q_b.
        return QVector(-I*gb, -gb, 1.) * (holscat/EXTRAPOLEN);
    }
    else {// If non of the above, then throw an error.
        throw std::out_of_range("In qvectorAtPos(): Index out of range.");
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
 * 
 * @deprecated This method provides a very bad ansatz to initialize the Newton-Raphson solver. It is only used for testing purposes.
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
    
    ga_avg = std::sqrt(ga*gb); // Geometric average of gamma values (equal to 1/sqrt(T)).
    
    // Assign constant values to all points in the mesh:
    theta0 = dcomplex(0., 0.); // Very coarse. Theta always varies in space.
    eta0 = std::atan(-I*ga_avg);    // Exact only in the absence of losses (no absorption, complete channel control).
    
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        field[0].set(2*ipoint,   theta0);
        field[0].set(2*ipoint+1, eta0);
    }
}

/**
 * Solve the Usadel System using the Newton-Raphson iterative algorithm.
 * This method does not call initializers (parameters (theta, eta) are not initialized), but the usage of an appropriate ansatz
 * is highly recommended in order to accelerate convergence while avoiding improper solutions.
 * 
 * Arguments:
 * 
 * maxit   = Maximum number of iterations. Typically: 50-500.
 * nsub    = Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
 *           The smallest reduction factor is thus 2^-nsub and should be larger than the machine epsilon (2^-53).
 * toldf   = Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7 or SQRTEPS=1e-8.
 * tolr    = Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10 or SQRTEPS=1e-8.
 * verbose = Verbosity level in standard output. 0=No output, 1=Display each iteration.
 * 
 * Returns:
 * 
 * niter   = Number of iterations used to meet the convergence criteria. If niter = maxit+1, then the iteration failed to reach convergence.
 */
int UsadelSystem::solveNewton(const int maxit, const int nsub, const double toldf, const double tolr, const int verbose) {
    
    int iter, s;    // Iteration indices.
    double fac, resnorm, resnorm0, newresnorm, fieldnorm, dfof, ror0;
    
    const int nparam = 2*npoint;          // Size of the vectors and matrix (number of columns, or number of rows).
    ComplexVector res(nparam), newfield(nparam), newres(nparam), delta(nparam);  // Allocate vectors: residual, new field, new residual, displacement.
    SparseComplexMatrix jac;  // Declare sparse Jacobian (no allocation).
    
    residualVector(field[0], res); // Compute the residual from the current field.
    
    resnorm = res.norm(); // Compute the norm of the residual.
    resnorm0 = std::max(resnorm, 1.);   // Store the norm of the first residual for the stopping criterion.
    
    for (iter = 1; iter <= maxit; iter++) {// Loop over the allowed Newton-Raphson iterations.
        
        // 1. Solves the Jacobian system using UMFPACK to find the Newton-Raphson displacement:
        computeJacobian(field[0], jac);  // Compute the sparse Jacobian using finite differences. This allocates memory only at the first iteration.
        solveUmfpack(jac, res, delta);  // Solve the system Jac.delta = resv, for the Newton-Raphson displacement "delta".
        
        // 2. Backtracking line search along the direction to minimize the residual (ensure convergence of Newton-Raphson):
        fac = -1.;  // Negative because the actual Newton-Raphson displacement is: newfield = field[0] - delta.
        for (s = 1; s <= nsub; s++) {
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
        fieldnorm = field[0].norm();
        dfof = delta.norm()/std::max(fieldnorm, 1.);
        ror0 = resnorm/resnorm0;
        
        if (verbose >= 1) {// Print the current status.
            std::cout << std::setprecision(4) 
                      << TAG_INFO << "#" << iter << "\t| f = " << fieldnorm << ",\tdf/f = " << dfof
                      << ",\tr/r0 = " << ror0 << ",\tfac = " << (fac!=-1 ? "1/" : "") << 1./std::abs(fac) << "\n";
        }
        
        // 4. Stopping criterion:
        if (dfof < toldf && ror0 < tolr) {
            if (verbose >= 1) {
                std::cout << TAG_INFO << "#" << iter << "\t| \033[92mFound solution\033[0m: df/f = "
                          << dfof << " < toldf = " << toldf << ", and r/r0 = " << ror0 << " < tolr = " << tolr << "\n";
            }
            break;
        }
    }
    
    if (iter > maxit && verbose >= 1) {
        std::cout << TAG_WARN << "#" << iter << "\t| \033[91mSolution not found\033[0m: df/f = " 
                  << dfof << " > toldf = " << toldf << ", and r/r0 = " << ror0 << " > tolr = " << tolr << "\n";
    }
    
    return iter;
}

/**
 * Save the mesh contained in the present UsadelSystem object.
 * @note The output of saveMesh() is now merged with that of saveField(). This function is thus no longer necessary.
 */
void UsadelSystem::saveMesh(const std::string& filename, const char* sep) const {
    mesh->saveMesh(filename, sep);
}

/**
 * Save the present field to a file of given "filename" using the string "sep" as a separator and printing "prec" decimals.
 * 
 * Columns of the output file are: x, y, north, south, east, west, theta_re, theta_im, eta_re, eta_im, q11_re, q11_im, q12_re, q12_im, q21_re, q21_im, I_a, I_b, C_ab
 * - Columns (x, y) are the coordinates of the points in the mesh.
 * - Columns (north, south, east, west) contains the indices of the neighbors if any. Otherwise, it contains the type of boundary condition for plotting.
 * - Columns (theta_re, theta_im, eta_re, eta_im) are the real and imaginary parts of the standard angular parameters (theta, eta).
 * - Columns (q11_re, q11_im, q12_re, q12_im, q21_re, q21_im) are the real and imaginary parts of the components of the Q field.
 * - Columns (I_a, I_b, C_ab) are the normalized intensity profiles deduced from Re(Q_12), Re(Q_21), and Im(Q_11), respectively.
 *   Theses intensity profiles are normalized such that the incident intensities in the input and output channels are equal to 1.
 *   This fixes the absolute normalization of the intensity profiles and allow for comparison with simulations based on the wave equation.
 * 
 * Of course, this function also assumes the Q field is found.
 * 
 * Arguments:
 * 
 * filename = The output filename, preferably with the '.csv' extension.
 * sep      = The separator string to use. Typically ',' for a CSV file.
 * prec     = Number of decimal places desired. Typically 16 in double precision.
 */
void UsadelSystem::saveField(const std::string& filename, const char* sep, const int prec) const {
    
    MeshPoint p;
    dcomplex theta, eta;
    QVector qv;
    int Na, Nb, Nv;
    double rho, Ia, Ib, Cab;
    
    Na = mesh->getNBoundary(BND_INPUT);   // Number of input channels.
    Nb = mesh->getNBoundary(BND_OUTPUT);  // Number of output channels.
    Nv = std::min(Na, Nb);  // The number of eigenvalues is the minimum between the number of channels at the two ports (input/output).
    rho = getRho();  // Extract the density of eigenvalues.
    
    if (rho < 0.) {// Warn the user about improper values of the transmission eigenvalue density rho(T).
        std::cout << TAG_WARN << "Density rho=" << rho << " is negative. This may indicate an improper solution.\n";
    } else if (rho < RHOMIN) {
        std::cout << TAG_WARN << "Density rho=" << rho << " is very small (<" << RHOMIN << "). This may indicate the gap of rho(T).\n";
    }
    
    std::ofstream ofs;  // Declare output stream object.
    ofs.open(filename.c_str()); // Open the file in write mode.
    
    const auto default_precision = ofs.precision(); // Saves the default precision.
    ofs << std::setprecision(prec); // Set the printing precision.
    
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    ofs << "%% Parameters: name=" << name << ", Npoint=" << npoint << ", h/lscat=" << holscat << ", h/labso=" << holabso 
        << ", Na=" << Na << ", Nb=" << Nb << ", Tval=" << tval << ", rho=" << rho << "\n"
        << "x" << sep << "y" << sep << "north" << sep << "south" << sep << "east" << sep << "west" << sep 
        << "theta_re" << sep << "theta_im" << sep << "eta_re" << sep << "eta_im" << sep 
        << "q11_re" << sep << "q11_im" << sep << "q12_re" << sep << "q12_im" << sep << "q21_re" << sep << "q21_im" << sep 
        << "I_a" << sep << "I_b" << sep << "C_ab\n";
    
    for (int ipoint = 0; ipoint < npoint; ipoint++) {// Loop over the points of the mesh.
        
        p = mesh->getPoint(ipoint); // Extract the point to get its coordinates.
        theta = field[0][2*ipoint]; // Extract the angle parameters.
        eta   = field[0][2*ipoint + 1];
        qv = QVector(theta, eta);  // Compute the Q field (in vector form) at the given point.
        
        Ia = Na*gb.real()*qv.getQ12().real()/(PI*Nv*rho);           // Eigenchannel intensity profile starting from input port.
        Ib = Nb*ga.real()*qv.getQ21().real()/(PI*Nv*rho);           // Eigenchannel intensity profile starting from output port.
        Cab = std::sqrt(Na*Nb/tval)*qv.getQ11().imag()/(PI*Nv*rho); // Correlator between input and output eigenchannels.
        
        ofs << p.x << sep << p.y << sep
            << boundaryTypeString(p.north) << sep << boundaryTypeString(p.south) << sep 
            << boundaryTypeString(p.east)  << sep << boundaryTypeString(p.west)  << sep
            << theta.real() << sep << theta.imag() << sep 
            << eta.real() << sep << eta.imag() << sep 
            << qv.getQ11().real() << sep << qv.getQ11().imag() << sep
            << qv.getQ12().real() << sep << qv.getQ12().imag() << sep
            << qv.getQ21().real() << sep << qv.getQ21().imag() << sep
            << Ia << sep << Ib << sep << Cab << "\n";
    }
    
    ofs.close();  // Close the file.
    
    ofs << std::setprecision(default_precision); // Restore default precision for printing.
}

/**
 * Save all the contextual data of the present UsadelSystem object in view of plotting the field.
 * This method also checks if the given path already exists. If directories are missing, then create them.
 * If target files already exist, then increase the counter until unique filenames are found.
 */
void UsadelSystem::savePlot(const std::string& path) const {
    
    std::string filename_field; // Target filename.
    
    uniqueFilename(path, ".csv", filename_field);  // Create a unique filename "filename_field". The result is of the form "<path><number><suffix>".
    
    std::cout << TAG_INFO << "Saving fields to file: '" << filename_field << "'...\n";
    saveField(filename_field, ", ", 16);
    
    std::string cmd("plot/plot_map.py lin I_a " + filename_field);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}
