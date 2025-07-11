/****
 * @date Created on 2025-07-06 at 19:44:15 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the ComplexVector methods.
 ***/
#include "ComplexVector.hpp"
#include <cstring>

/**
 * Constructor of the ComplexVector object.
 */
ComplexVector::ComplexVector(const int64_t size) {
    if (size < 0) {
        throw std::invalid_argument("In ComplexVector(): Received invalid negative size.");
    }
    //std::cout << TAG_INFO << "Creating ComplexVector (allocated but not initialized).\n";
    this->size = size;
    real = new double[size]();  // Zero initialization.
    imag = new double[size]();
}

/**
 * Copy constructor for the ComplexVector object.
 * Used implicitly in expressions such as c = a + b for instance.
 */
ComplexVector::ComplexVector(const ComplexVector& vec) {
    //std::cout << TAG_INFO << "Creating ComplexVector using copy constructor.\n";
    size = vec.size;
    real = new double[size];
    imag = new double[size];
    std::memcpy(real, vec.real, size*sizeof(double)); // Deep copy of the array (should be faster than loop).
    std::memcpy(imag, vec.imag, size*sizeof(double)); // Deep copy of the array (should be faster than loop).
}

/**
 * Constructor of ComplexVector object taking in argument an initializer list.
 */
ComplexVector::ComplexVector(const std::initializer_list<dcomplex> list) {
    //std::cout << TAG_INFO << "Creating ComplexVector from initializer list.\n";
    size = list.size();
    real = new double[size];
    imag = new double[size];
    
    int i = 0;
    for (auto p = list.begin(); p < list.end(); p++) {
        real[i] = (*p).real();
        imag[i] = (*p).imag();
        i++;
    }
}

/**
 * Destructor of the ComplexVector object.
 */
ComplexVector::~ComplexVector() {
    //std::cout << TAG_INFO << "Deleting ComplexVector.\n";
    delete[] real;
    delete[] imag;
}

/**
 * Overloads the assignement operator to perform a deep copy.
 */
ComplexVector& ComplexVector::operator=(const ComplexVector& vec) {
    if (this == &vec) {
        throw std::invalid_argument("In ComplexVector operator=(): Invalid self assignment.");
    }
    if (this->size != vec.size) {
        std::string info = "In ComplexVector operator=(): Vectors have different sizes. LHS = " 
                         + std::to_string(this->size) + ", RHS = " + std::to_string(vec.size) + ".";
        throw std::invalid_argument(info);
    }
    //std::cout << TAG_INFO << "ComplexVector operator=(): Performs memory copy.\n";
    std::memcpy(this->real, vec.real, size*sizeof(double)); // Deep copy of the array (should be faster than loop).
    std::memcpy(this->imag, vec.imag, size*sizeof(double)); // Deep copy of the array (should be faster than loop).
    
    return *this;
}

/**
 * Returns the i-th complex element of the ComplexVector for read only.
 */
dcomplex ComplexVector::operator[](const int64_t i) const {
    if (i < 0 || i >= size) {
        std::string info = "In ComplexVector operator[]: Index out of range. Received " + std::to_string(i) 
                         + ", expected in [0, " + std::to_string(size-1) + "]";
        throw std::out_of_range(info);
    }
    return dcomplex(real[i], imag[i]);
}

/**
 * Assigns an element of the ComplexVector.
 */
void ComplexVector::set(const int64_t i, const dcomplex z) {
    if (i < 0 || i >= size) {
        std::string info = "In ComplexVector set(): Index out of range. Received " + std::to_string(i) 
                         + ", expected in [0, " + std::to_string(size-1) + "]";
        throw std::out_of_range(info);
    }
    real[i] = z.real();
    imag[i] = z.imag();
}

/**
 * Add the given constant at some position in the ComplexVector.
 * This is equivalent to: vec.set(i, vec[i] + z)
 */
void ComplexVector::add(const int64_t i, const dcomplex z) {
    if (i < 0 || i >= size) {
        std::string info = "In ComplexVector set(): Index out of range. Received " + std::to_string(i) 
                         + ", expected in [0, " + std::to_string(size-1) + "]";
        throw std::out_of_range(info);
    }
    real[i] += z.real();
    imag[i] += z.imag();
}

/**
 * Returns the size of the ComplexVector.
 */
int64_t ComplexVector::getSize() const {
    return size;
}

/**
 * Returns the norm of the ComplexVector.
 */
double ComplexVector::norm() const {
    double n2 = 0.;
    for (int64_t i = 0; i < size; i++) {
        n2 += real[i]*real[i] + imag[i]*imag[i];
    }
    return std::sqrt(n2);
}

/**
 * Overload the addition operator between two ComplexVectors.
 */
ComplexVector ComplexVector::operator+(const ComplexVector& vec) const {
    if (size != vec.size) {
        throw std::invalid_argument("In ComplexVector operator+: Vectors have different sizes.");
    }
    ComplexVector result(*this);
    for (int i = 0; i < size; i++) {
        result.real[i] += vec.real[i];
        result.imag[i] += vec.imag[i];
    }
    return result;
}

/**
 * Overload the subtraction operator between two ComplexVectors.
 */
ComplexVector ComplexVector::operator-(const ComplexVector& vec) const {
    if (size != vec.size) {
        throw std::invalid_argument("In ComplexVector operator-: Vectors have different sizes.");
    }
    ComplexVector result(*this);
    for (int i = 0; i < size; i++) {
        result.real[i] -= vec.real[i];
        result.imag[i] -= vec.imag[i];
    }
    return result;
}

/**
 * Overload the multiply operator to perform the multiplication of a ComplexVector by a scalar number.
 */
ComplexVector ComplexVector::operator*(const dcomplex scalar) const {
    dcomplex product;
    ComplexVector result(*this);
    for (int i = 0; i < size; i++) {
        product = scalar * dcomplex(result.real[i], result.imag[i]);
        result.real[i] = product.real();
        result.imag[i] = product.imag();
    }
    return result;
}

/**
 * Multiplication of a scalar number by a vector.
 */
ComplexVector operator*(const dcomplex scalar, const ComplexVector& vec) {
    return vec * scalar;
}

/**
 * Overload the stream operator to print the ComplexVector.
 */
std::ostream& operator<<(std::ostream& os, const ComplexVector& vec) {
    os << "[";
    for (int64_t i = 0; i < vec.size; i++) {
        os << " " << vec.real[i] << (vec.imag[i] >= 0. ? "+" : "") << vec.imag[i] << "i ";
    }
    os << "]";
    return os;
}

