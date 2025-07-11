/****
 * @date Created on 2025-07-02 at 09:34:47 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the Usadel System.
 ***/
#include "UsadelSystem.hpp"
#include "BaseTools.hpp"

/**
 * Creates the reference UsadelSystem used for testing purposes.
 */
UsadelSystem usadelSystemRef1(SquareMesh& mesh, std::vector<Contact>& contact) {
    
    double tm, holscat, holabso;
    
    mesh.addRectangle(-10, 10, -5, 5);
    mesh.setBoundaryRegion(-10, -10, -5, 5, BND_OPEN);
    mesh.setBoundaryRegion( 10,  10, -5, 5, BND_OPEN);
    mesh.fixNeighbors();
    
    tm = 0.5; // Transmission probability.
    dcomplex ga = dcomplex(sqrt(1./tm), 1.e-7);  // Value of gamma_a = gamma_b.
    Contact cta = Contact(Vector2D(-9.5, 5.5), Vector2D(-9.5, -5.5), +1, ga);
    Contact ctb = Contact(Vector2D(9.5, 5.5), Vector2D(9.5, -5.5), -1, ga);
    contact = {cta, ctb};  // Create contact interactions.
    
    holscat = 0.1;    // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.01;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    
    return UsadelSystem(mesh, contact, holscat, holabso);
}

/**
 * Try to extract the Q field at some point in the middle of the mesh.
 */
int extractQPoint() {
    std::cout << "====== TEST EXTRACT Q POINT ======\n";
    
    SquareMesh mesh;
    std::vector<Contact> contact;
    UsadelSystem usys = usadelSystemRef1(mesh, contact);
    
    const uint64_t seed = 1;
    usys.initRandom(seed);  // Set the Q field to random to avoid trivial output.
    
    int i0 = usys.getNPoint()/2;
    
    std::cout << TAG_INFO << "At pos i=" << i0 << ", QVector = " << usys.getQVector(i0) << "\n";
    
    return 0;
}

/**
 * Try to extract the whole Q field from the UsadelSystem object.
 */
int extractQField() {
    std::cout << "====== TEST EXTRACT Q FIELD ======\n";
    
    SquareMesh mesh;
    std::vector<Contact> contact;
    UsadelSystem usys = usadelSystemRef1(mesh, contact);
    
    const uint64_t seed = 1;
    usys.initRandom(seed);  // Set the Q field to random to avoid trivial output.
    
    int npoint = usys.getNPoint();
    
    for (int i = 0; i < npoint; i++) {
        std::cout << TAG_INFO << "At pos i=" << i << ", QVector = " << usys.getQVector(i) << "\n";
    }
    
    return 0;
}

/**
 * Test the residual of the Usadel equation.
 */
int testResidual() {
    std::cout << "====== TEST USADEL RESIDUAL ======\n";
    
    SquareMesh mesh;
    std::vector<Contact> contact;
    UsadelSystem usys = usadelSystemRef1(mesh, contact);
    
    const uint64_t seed = 1; // Initialize the Q field to random values to avoid trivial output.
    usys.initRandom(seed);
    
    const double tolerr = 1e-13; // Tolerance used for comparison.
    
    return usys.testResidual(tolerr);
}

/**
 * Test the residual of the Usadel equation.
 */
int testJacobian() {
    std::cout << "====== TEST USADEL JACOBIAN ======\n";
    
    SquareMesh mesh;
    std::vector<Contact> contact;
    UsadelSystem usys = usadelSystemRef1(mesh, contact);
    
    const uint64_t seed = 1; // Initialize the Q field to random values to avoid trivial output.
    usys.initRandom(seed);
    
    const double tolerr = 1e-13; // Tolerance used for comparison.
    
    return usys.testJacobian(tolerr);
}

/**
 * Test solving hte Usadel equation using Newton-Raphson iterative method.
 * In a reference situation: Quasi-one-dimensional and translationally invariant in y direction.
 */
int testSolve() {
    std::cout << "====== TEST SOLVE ======\n";
    
    int maxit, nsub, verbose;
    double tm, holscat, holabso, tolp, tolr;
    
    // Create a square mesh:
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    mesh.setBoundaryRegion(-30, -30, -14, 14, BND_OPEN);
    mesh.setBoundaryRegion( 30,  30, -14, 14, BND_OPEN);
    mesh.fixNeighbors();
    
    // Contact interactions:
    tm = 0.99; // Transmission probability.
    dcomplex ga = dcomplex(sqrt(1./tm), 1.e-7);  // Value of gamma_a = gamma_b.
    Contact cta = Contact(Vector2D(-30.5, 15.5), Vector2D(-30.5, -15.5), +1, ga);
    Contact ctb = Contact(Vector2D(30.5, 15.5), Vector2D(30.5, -15.5), -1, ga);
    std::vector<Contact> contact = {cta, ctb};
    
    holscat = 0.083; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    UsadelSystem usys(mesh, contact, holscat, holabso);
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    usys.initConstant();
    usys.solveNewton(maxit, nsub, tolp, tolr, verbose);
    
    const char* filename = "test_field_x.csv";
    std::cout << TAG_INFO << "Saving data to file '" << filename << "'...\n";
    usys.saveField(filename, ", ", 16);
    
    return 0;
}

/**
 * Test saving the field.
 */
int testSaveField() {
    std::cout << "====== TEST SAVE FIELD ======\n";
    
    SquareMesh mesh;
    std::vector<Contact> contact;
    UsadelSystem usys = usadelSystemRef1(mesh, contact);
    
    const uint64_t seed = 1;
    usys.initRandom(seed);  // Set the Q field to random to avoid trivial output.
    
    const char* filename = "test_save_field.csv";
    std::cout << TAG_INFO << "Saving random field to file '" << filename << "'...\n";
    usys.saveField(filename, ", ", 16);
    
    return 0;
}

/**
 * Main function of the test of the UsadelSystem object.
 */
int main(int argc, char** argv) {
    
    std::cout << "This is " << PROGRAM_COPYRIGHT << "\n";
    
    testResidual();
    testJacobian();
    testSolve();
    //testSaveField();
    
    return 0;
}
