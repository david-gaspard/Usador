/****
 * @date Created on 2025-07-03 at 13:38:38 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the Vector2D methods.
 ***/
#include "Vector2D.hpp"
#include <cmath>

/**
 * Constructor of the Vector2D object.
 */
Vector2D::Vector2D(int x, int y) : x(x), y(y) {}

/**
 * Constructor of the Vector2D structure.
 */
Vector2D::Vector2D(double x, double y) : x(x), y(y) {}

/**
 * Constructor of the Vector2D structure from a mesh point.
 */
Vector2D::Vector2D(MeshPoint p) : x(p.x), y(p.y) {}

/**
 * Default constructor of the Vector2D structure.
 */
Vector2D::Vector2D() : x(0.), y(0.) {}

/**
 * Addition of two vectors.
 */
Vector2D Vector2D::operator+(const Vector2D& v) const {
    return Vector2D(x + v.x, y + v.y);
}

/**
 * Subtraction of two vectors.
 */
Vector2D Vector2D::operator-(const Vector2D& v) const {
    return Vector2D(x - v.x, y - v.y);
}

/**
 * Norm of a vector.
 */
double Vector2D::norm() const {
    return sqrt(x*x + y*y);
}

/**
 * Dot (scalar) product of two vectors.
 */
double Vector2D::dot(const Vector2D& v) const {
    return x * v.x + y * v.y;
}

/**
 * Cross (pseudo-scalar) product of two vectors.
 */
double Vector2D::cross(const Vector2D& v) const {
    return x * v.y - y * v.x;
}

/**
 * Overloads the stream operator for printing the vector.
 */
std::ostream& operator<<(std::ostream& os, const Vector2D& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}
