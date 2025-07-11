/****
 * @date Created on 2025-07-07 at 14:09:58 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing test for the SparseComplexMatrix object.
 ***/
#include "SparseComplexMatrix.hpp"

/**
 * Test solving a simple system of equations. First version fuond in UMFPACK's Manual (https://fossies.org/linux/SuiteSparse/UMFPACK/Doc/UMFPACK_UserGuide.pdf).
 */
int testSolve1() {
    std::cout << "====== TEST SOLVE UMFPACK #1 ======\n";
    
    // 1. Initialize a sparse complex matrix:
    int64_t n = 5;
    
    int64_t ti[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
    int64_t tj[] = {0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4};
    double treal[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
    double timag[] = {0., 0., 0.,  0., 0., 0.,  0., 0., 0., 0., 0., 0.};
    
    int64_t nnz = sizeof(ti)/sizeof(int64_t);
    
    SparseComplexMatrix a;
    
    a.set(n, n, nnz, ti, tj, treal, timag, true); // Initialize the matrix using triplets.
    //a.print(5); // Print full matrix using UMFPACK.
    
    // 2. Initialize an independent vector term:
    ComplexVector b = {8., 45., -3., 3., 19.};
    ComplexVector x(n);
    
    solveUmfpack(a, b, x); // Solve the system A.x = b.
    
    // 3. Check the validity of the solution using the expected solution:
    ComplexVector x_expc = {1., 2., 3., 4., 5.};  // Expected solution.
    
    std::cout << TAG_INFO << "Found solution: x = " << x << "\n";
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm()/x_expc.norm() << "\n";
    
    return 0;
}

/**
 * Test solving a simple system of equations. Second version customized to test for complex numbers.
 */
int testSolve2() {
    std::cout << "====== TEST SOLVE UMFPACK #2 ======\n";
    
    // 1. Initialize a sparse complex matrix:
    int64_t n = 5;
    
    int64_t ti[] = {0, 1, 1, 2, 2, 2, 3, 3, 4};
    int64_t tj[] = {0, 0, 1, 1, 2, 3, 3, 4, 4};
    double treal[] = {-1.6, -8.7, 8.2,  7.5, -4.2, -8.2, 1.3, 9.0, 8.7};
    double timag[] = {-9.7, -0.3, 2.9, -7.3, -1.4, -7.5, 2.4, 6.3, 2.5};
    
    int64_t nnz = sizeof(ti)/sizeof(int64_t);
    
    SparseComplexMatrix a;
    
    a.set(n, n, nnz, ti, tj, treal, timag, true); // Initialize the matrix using triplets.
    //a.print(5); // Print full matrix using UMFPACK.
    
    // 2. Initialize an independent vector term:
    ComplexVector b = {dcomplex(-3.1,-0.7), dcomplex(0.4,2.9), dcomplex(3.7,1.4), dcomplex(3.3,-1.0), dcomplex(-1.2,1.4)};
    ComplexVector x(n);
    
    solveUmfpack(a, b, x); // Solve the system A.x = b.
    
    // 3. Check the validity of the solution using the expected solution:
    ComplexVector x_expc = {dcomplex(0.1215726849456803,-0.2995344024831868), dcomplex(0.18041460985307842,-0.02349754447487198), dcomplex(-3.363564620069479,4.046269756344257), dcomplex(0.22513379408406553,-2.0569639267887947), dcomplex(-0.08469611911154504,0.1852575054918233)};  // Expected solution.
    
    std::cout << TAG_INFO << "Found solution: x = " << x << "\n";
    std::cout << TAG_INFO << "Error = " << (x - x_expc).norm()/x_expc.norm() << "\n";
    
    return 0;
}

/**
 * Main function of the tests of the SparseComplexMatrix object.
 */
int main(int argc, char** argv) {
    
    testSolve1();
    testSolve2();
    
    return 0;
}
