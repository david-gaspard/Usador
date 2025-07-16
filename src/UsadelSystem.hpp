/****
 * @date Created on 2025-07-01 at 13:47:33 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the Usadel object which describes the system to be solved.
 ***/
#ifndef _USADEL_SYSTEM_H
#define _USADEL_SYSTEM_H
#include "SquareMesh.hpp"
#include "Contact.hpp"
#include "SparseComplexMatrix.hpp"

/**
 * Class defining the Square Mesh object.
 */
class UsadelSystem {
    
    private:
    
    SquareMesh* mesh;  // Pointer to the mesh.
    std::vector<Contact>* contact;  // Pointer to a vector of contact interactions.
    double holscat;   // Mesh step divided by the scattering mean free path (total length is not a well defined unit).
    double holabso;   // Mesh step divided by the absorption length (total length is not a well defined unit).
    int npoint;       // Total number of points in the "mesh".
    ComplexVector* field;  // Array containing the parameters (theta, eta) at each corresponding point in the mesh. Size: 2*npoint.
    
    public:
    
    // Constructors/Destructors:
    UsadelSystem(SquareMesh& mesh, std::vector<Contact>& contact, const double holscat, const double holabso);
    UsadelSystem(const int length, const int width, const double dscat, const double dabso);
    ~UsadelSystem();
    
    // Getters:
    int getNPoint() const;
    MeshPoint getPoint(const int ipoint) const;
    int indexOf(const int x, const int y) const;
    double getHolscat() const;
    double getHolabso() const;
    QVector getQVector(const int ipoint) const;
    QVector getJVector(const int ipoint, const Direction dir) const;
    double getRho() const;
    
    // Initializers/Setters:
    void setTransmission(const double tval);
    void initRandom(const uint64_t seed);
    void initConstant();
    void initWaveguide();
    
    // Testing functions:
    int testResidual(const double tolerr) const;
    int testJacobian(const double tolerr) const;
    
    // Public computational functions:
    int solveNewton(const int maxit, const int nsub, const double tolp, const double tolr, const int verbose);
    
    // Output functions:
    void saveMesh(const char* filename, const char* sep, const int verbosity) const;
    void saveField(const char* filename, const char* sep, const int prec) const;
    // TODO: Plot methods for boundaries and contact interactions......
    
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
