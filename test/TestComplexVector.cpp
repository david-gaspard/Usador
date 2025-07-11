/****
 * @date Created on 2025-07-06 at 20:17:03 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing test for the ComplexVector class.
 ***/
#include "ComplexVector.hpp"
#include <iostream>

/**
 * Perform some basic test on the ComplexVector object.
 */
int basicTest() {
    
    ComplexVector a = {dcomplex(1., 2.), dcomplex(3., 4.), dcomplex(5., 6.), dcomplex(7., 8.), dcomplex(9., 10.)};  // Assignment from array.
    
    std::cout << TAG_INFO << "a = " << a << "\n";
    std::cout << TAG_INFO << "norm(a) = " << a.norm() << "\n";
    
    ComplexVector b = a;  // Call copy constructor.
    
    b.set(1, dcomplex(-2.5, -3.5)); // Assign a given vector element.
    
    std::cout << TAG_INFO << "b = " << b << "\n";
    std::cout << TAG_INFO << "norm(b) = " << b.norm() << "\n";
    
    std::cout << TAG_INFO << "a + b = " << a + b << "\n";
    std::cout << TAG_INFO << "a - b = " << a - b << "\n";
    std::cout << TAG_INFO << "(2-I)*a = " << (2. - I)*a << "\n";
    std::cout << TAG_INFO << "a*(2-I) = " << a*(2. - I) << "\n";
    
    return 0;
}

/**
 * Main function of the test of the ComplexVector class.
 */
int main(int argc, char** argv) {
    
    basicTest();
    
    return 0;
}
