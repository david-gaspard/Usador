/****
 * @date Created on 2025-07-01 at 17:37:18 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the methods for the QVector structure.
 ***/
#include "QVector.hpp"

/**
 * Constructor of the QVector structure.
 */
QVector::QVector(const dcomplex q1, const dcomplex q2, const dcomplex q3) : q1(q1), q2(q2), q3(q3) {}

/**
 * Constructor of the QVector structure.
 */
QVector::QVector() : q1(0.), q2(0.), q3(0.) {}

/**
 * Returns the 3-vector representation of the Q field from the given standard angular parameters (theta, eta).
 * This function returns the vector: q(x) = [sin(eta)*cos(theta), -sin(theta), cos(eta)*cos(theta)].
 */
QVector::QVector(const dcomplex theta, const dcomplex eta) {
    dcomplex costheta = std::cos(theta);
    q1 = std::sin(eta) * costheta;
    q2 = -std::sin(theta);
    q3 = std::cos(eta) * costheta;
}

/**
 * Returns the polar angle parameter "theta" of the current QVector.
 * This function does not assume the QVector is normalized.
 */
dcomplex QVector::getTheta() const {
    return std::asin(-q2/std::sqrt(q1*q1 + q2*q2 + q3*q3));
}

/**
 * Returns the azimuthal angle parameter "eta" of the current QVector.
 * This function does not assume the QVector is normalized.
 */
dcomplex QVector::getEta() const {
    return -I*std::log((q3 + I*q1)/std::sqrt(q1*q1 + q3*q3));
}

/**
 * Return the Q_11 component of Q in 2-by-2 matrix form.
 * The Q matrix is given by: Q = q1 * Sigma_1 + q2 * Sigma_2 + q3 * Sigma_3 ,
 * where Sigma_1, Sigma_2, Sigma_3 are the three standard Pauli matrices.
 */
dcomplex QVector::getQ11() const {
    return q3;
}

/**
 * Return the Q_12 component of Q in 2-by-2 matrix form.
 * The Q matrix is given by: Q = q1 * Sigma_1 + q2 * Sigma_2 + q3 * Sigma_3 ,
 * where Sigma_1, Sigma_2, Sigma_3 are the three standard Pauli matrices.
 */
dcomplex QVector::getQ12() const {
    return q1 - I*q2;
}

/**
 * Return the Q_21 component of Q in 2-by-2 matrix form.
 * The Q matrix is given by: Q = q1 * Sigma_1 + q2 * Sigma_2 + q3 * Sigma_3 ,
 * where Sigma_1, Sigma_2, Sigma_3 are the three standard Pauli matrices.
 */
dcomplex QVector::getQ21() const {
    return q1 + I*q2;
}

/**
 * Addition of two qvectors.
 */
QVector QVector::operator+(const QVector& v) const {
    return QVector(q1 + v.q1, q2 + v.q2, q3 + v.q3);
}

/**
 * Difference of two qvectors.
 */
QVector QVector::operator-(const QVector& v) const {
    return QVector(q1 - v.q1, q2 - v.q2, q3 - v.q3);
}

/**
 * Complex norm of a qvector.
 */
dcomplex QVector::norm() const {
    return std::sqrt(q1*q1 + q2*q2 + q3*q3);
}

/**
 * Scalar (dot) product of two qvectors.
 */
dcomplex QVector::dot(const QVector& v) const {
    return q1*v.q1 + q2*v.q2 + q3*v.q3;
}

/**
 * Vector (cross) product of two qvectors.
 */
QVector QVector::cross(const QVector& v) const {
    return QVector(q2 * v.q3 - q3 * v.q2, q3 * v.q1 - q1 * v.q3, q1 * v.q2 - q2 * v.q1);
}

/**
 * Rotate the qvector around the Sigma_1 axis.
 * In matrix form, this applies the transformation Q -> exp(-i * Sigma_1 * alpha/2) * Q * exp(i * Sigma_1 * alpha/2) .
 */
QVector QVector::rotate_1(const dcomplex alpha) const {
    dcomplex cosa, sina, r2, r3;
    cosa = std::cos(alpha);
    sina = std::sin(alpha);
    r2 = cosa * q2 - sina * q3;
    r3 = sina * q2 + cosa * q3;
    return QVector(q1, r2, r3);
}

/**
 * Rotate the qvector around the Sigma_2 axis.
 * In matrix form, this applies the transformation Q -> exp(-i * Sigma_2 * alpha/2) * Q * exp(i * Sigma_2 * alpha/2) .
 */
QVector QVector::rotate_2(const dcomplex alpha) const {
    dcomplex cosa, sina, r1, r3;
    cosa = std::cos(alpha);
    sina = std::sin(alpha);
    r3 = cosa * q3 - sina * q1;
    r1 = sina * q3 + cosa * q1;
    return QVector(r1, q2, r3);
}

/**
 * Rotate the qvector around the Sigma_3 axis.
 * In matrix form, this applies the transformation Q -> exp(-i * Sigma_3 * alpha/2) * Q * exp(i * Sigma_3 * alpha/2) .
 */
QVector QVector::rotate_3(const dcomplex alpha) const {
    dcomplex cosa, sina, r1, r2;
    cosa = std::cos(alpha);
    sina = std::sin(alpha);
    r1 = cosa * q1 - sina * q2;
    r2 = sina * q1 + cosa * q2;
    return QVector(r1, r2, q3);
}

/**
 * In matrix form, this applies the transformation Q -> exp(a * Sigma_+) * Q * exp(-a * Sigma_+) ,
 * where Sigma_+ is the raising Pauli matrix.
 */
QVector QVector::shear_plus(const dcomplex a) const {
    dcomplex s1, s2, s3, aq21, aaq;
    aq21 = a*(q1 + I*q2);  // Invariant matrix component.
    aaq = a*(q3 + aq21/2.);
    s1 = q1 -   aaq;
    s2 = q2 - I*aaq;
    s3 = q3 + aq21;
    return QVector(s1, s2, s3);
}

/**
 * In matrix form, this applies the transformation Q -> exp(b * Sigma_-) * Q * exp(-b * Sigma_-) ,
 * where Sigma_- is the lowering Pauli matrix.
 */
QVector QVector::shear_minus(const dcomplex b) const {
    dcomplex s1, s2, s3, bq12, bbq;
    bq12 = b*(q1 - I*q2);  // Invariant matrix element.
    bbq = b*(q3 - bq12/2.);
    s1 = q1 +   bbq;
    s2 = q2 - I*bbq;
    s3 = q3 - bq12;
    return QVector(s1, s2, s3);
}

/**
 * Overloads the stream operator for printing the QVector.
 */
std::ostream& operator<<(std::ostream& os, const QVector& v) {
    os << "[ " << v.q1 << "  " << v.q2 << "  " << v.q3 << " ]";
    return os;
}
