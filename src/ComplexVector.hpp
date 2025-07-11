/****
 * @date Created on 2025-07-06 at 19:11:48 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the complex vector object.
 ***/
#ifndef _COMPLEXVECTOR_H
#define _COMPLEXVECTOR_H
#include <Constants.hpp>
#include <iostream>

class SparseComplexMatrix;  // Forward declaration is needed by friend function.

/**
 * Class defining the ComplexVector object.
 */
class ComplexVector {
    
    private:
    
    int64_t size;  // Number of complex elements.
    double* real;  // Array containing the real parts. Size: n.
    double* imag;  // Array containing the imaginary parts. Size: n.
    
    public:
    
    // Constructors/Destructors:
    ComplexVector(const int64_t size);
    ComplexVector(const ComplexVector& vec);
    ComplexVector(const std::initializer_list<dcomplex> c);
    ~ComplexVector();
    
    // Setters:
    ComplexVector& operator=(const ComplexVector& vec);
    void set(const int64_t i, const dcomplex z);
    void add(const int64_t i, const dcomplex z);
    
    // Getters:
    dcomplex operator[](const int64_t i) const;
    int64_t getSize() const;
    double norm() const;
    
    // Other operations:
    ComplexVector operator+(const ComplexVector& vec) const;
    ComplexVector operator-(const ComplexVector& vec) const;
    ComplexVector operator*(const dcomplex scalar) const;
    friend ComplexVector operator*(const dcomplex scalar, const ComplexVector& vec);
    friend std::ostream& operator<<(std::ostream& os, const ComplexVector& vec);
    friend void solveUmfpack(const SparseComplexMatrix& a, const ComplexVector& b, ComplexVector& x); // Promote access to the private fields to the linear solver.
    
};

#endif
