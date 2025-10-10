/****
 * @date Created on 2025-07-01 at 13:47:33 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the Usadel object which describes the system to be solved.
 ***/
#ifndef _USADEL_SYSTEM_H
#define _USADEL_SYSTEM_H
#include "SquareMesh.hpp"
#include "QVector.hpp"
#include "SparseComplexMatrix.hpp"

/**
 * Class defining the Square Mesh object.
 */
class UsadelSystem {
    
    private:
    
    SquareMesh* mesh; // Pointer to the square mesh.
    ComplexVector* field;  // Array containing the two complex parameters (theta, eta) at each corresponding point in the mesh. Size: 2*npoint.
    int npoint;       // Total number of points in the "mesh".
    
    double holscat;   // Mesh step divided by the scattering mean free path (total length is not a well defined unit).
    double holabso;   // Mesh step divided by the absorption length (total length is not a well defined unit).
    double tval;      // Transmission eigenvalue, between 0 and 1.
    dcomplex ga;      // Input contact parameter, gamma_a, given by sqrt(1/tval) + I*0^+.
    dcomplex gb;      // Output contact parameter, gamma_b, given by sqrt(1/tval) + I*0^+. They must be related by ga*gb = 1/tval + I*0^+.
    
    std::string name; // String containing the name of the system (typically describing the system geometry) which is used to generate file output.
    double solvetime;  // Computation time of the last call to solveNewton() in seconds.
    
    public:
    
    // Constructors/Destructors:
    UsadelSystem(const std::string& name, SquareMesh& mesh, const double holscat, const double holabso, const double tval);
    UsadelSystem(const std::string& name, const int length, const int width, const double dscat, const double dabso, const double tval);
    UsadelSystem(const UsadelSystem& usys); // Explicit copy constructor.
    ~UsadelSystem();
    
    // Getters:
    int getNPoint() const;
    MeshPoint getPoint(const int ipoint) const;
    int indexOf(const int x, const int y) const;
    double getHolscat() const;
    double getHolabso() const;
    double getTransmission() const;
    std::string getName() const;
    QVector getQVector(const int ipoint) const;
    QVector getJVector(const int ipoint, const int jpoint) const;
    double getRho() const;
    
    // Initializers/Setters:
    void setTransmission(const double tval);
    void copy(const UsadelSystem& usys);
    void initRandom(const uint64_t seed);
    void initConstant();
    void initWaveguide();
    
    // Testing functions:
    int testResidual(const double tolerr) const;
    int testJacobian(const double tolerr) const;
    
    // Public computational functions:
    int solveNewton(const int maxit, const int nsub, const double tolp, const double tolr, const int verbose);
    
    // Output functions:
    void saveMesh(const std::string& filename, const char* sep) const;
    void saveField(const std::string& filename, const char* sep, const int prec) const;
    void savePlot(const std::string& path) const;
    
    private:
    
    // Private computational methods:
    QVector qvectorAtPos(const ComplexVector& someField, int ipoint) const;
    void residualAtPos   (const ComplexVector& someField, int ipoint, dcomplex& res1, dcomplex& res2) const;
    void residualAtPosOld(const ComplexVector& someField, int ipoint, dcomplex& res1, dcomplex& res2) const;
    void residualVector(const ComplexVector& someField, ComplexVector& resvec) const;
    void computeJacobianAtPos(const ComplexVector& fieldPlus, const ComplexVector& fieldMinus, int64_t ipoint, int64_t jparam, 
                              int64_t& nnz, int64_t* irow, int64_t* icol, double* jacreal, double* jacimag) const;
    void computeJacobian(const ComplexVector& someField, SparseComplexMatrix& jac) const;
    void computeJacobianDense(const ComplexVector& someField, dcomplex* jacdense) const;
    
};

#endif
