/****
 * @date Created on 2025-07-01 at 17:40:56 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the QVector structure.
 ***/
#ifndef _QVECTOR_H
#define _QVECTOR_H
#include "Constants.hpp"
#include <iostream>

/**
 * Vector representation of the Q Field.
 */
struct QVector {
    
    dcomplex q1;  // Component of the Q(x) field multiplying the Pauli matrix Sigma_1 = [[0, 1],[1, 0]].
    dcomplex q2;  // Component of the Q(x) field multiplying the Pauli matrix Sigma_2 = [[0, -i], [i, 0]].
    dcomplex q3;  // Component of the Q(x) field multiplying the Pauli matrix Sigma_3 = [[1, 0], [0, 1]].
    
    // Constructors:
    QVector(const dcomplex q1, const dcomplex q2, const dcomplex q3);
    QVector(const dcomplex theta, const dcomplex eta);
    QVector();
    
    // Methods:
    dcomplex getTheta() const;
    dcomplex getEta() const;
    dcomplex getQ11() const;
    dcomplex getQ12() const;
    dcomplex getQ21() const;
    QVector operator+(const QVector& v) const;
    QVector operator-(const QVector& v) const;
    dcomplex norm() const;
    dcomplex dot(const QVector& v) const;
    QVector cross(const QVector& v) const;
    QVector rotate_1(const dcomplex alpha) const;
    QVector rotate_2(const dcomplex alpha) const;
    QVector rotate_3(const dcomplex alpha) const;
    QVector shear_plus(const dcomplex a) const;
    QVector shear_minus(const dcomplex b) const;
    
    friend std::ostream& operator<<(std::ostream& os, const QVector& v);
    
};

#endif
