/****
 * @date Created on 2025-07-02 at 09:34:47 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing tests for the Usadel System.
 ***/
#include "UsadelSystem.hpp"
#include "BaseTools.hpp"
#include <iomanip>

/**
 * Try to extract the Q field at some point in the middle of the mesh.
 */
int extractQPoint() {
    std::cout << "====== TEST EXTRACT Q POINT ======\n";
    
    const int length = 20;
    const int width = 10;
    
    UsadelSystem usys(length, width, 5., 0.1, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
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
    
    const int length = 20;
    const int width = 10;
    
    UsadelSystem usys(length, width, 5., 0.1, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
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
    
    const int length = 20;
    const int width = 10;
    
    UsadelSystem usys(length, width, 5., 0.1, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
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
    
    const int length = 20;
    const int width = 10;
    
    UsadelSystem usys(length, width, 5., 0.1, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
    const uint64_t seed = 1; // Initialize the Q field to random values to avoid trivial output.
    usys.initRandom(seed);
    
    const double tolerr = 1e-13; // Tolerance used for comparison.
    
    return usys.testJacobian(tolerr);
}

/**
 * Test saving the field.
 */
int testSaveField() {
    std::cout << "====== TEST SAVE FIELD ======\n";
    
    const int length = 20;
    const int width = 10;
    
    UsadelSystem usys(length, width, 5., 0.1, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
    const uint64_t seed = 1;
    usys.initRandom(seed);  // Set the Q field to random to avoid trivial output.
    
    const char* filename = "out/test/test_save_field.csv";
    std::cout << TAG_INFO << "Saving random field to file '" << filename << "'...\n";
    usys.saveField(filename, ", ", 16);
    
    return 0;
}

/**
 * Test solving the Usadel equation using Newton-Raphson iterative method.
 * In a reference situation: Quasi-one-dimensional and translationally invariant in y direction.
 */
int testSolve() {
    std::cout << "====== TEST SOLVE ======\n";
    
    int maxit, nsub, verbose;
    double tval, holscat, holabso, tolp, tolr;
    
    // Create a square mesh:
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    mesh.setBoundaryRegion(-30, -30, -14, 14, WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -14, 14, EAST, BND_OUTPUT);
    mesh.fixNeighbors();
    
    // Contact interactions:
    tval = 0.99; // Transmission probability.
    //Contact cta = Contact(Vector2D(-30.5, 15.5), Vector2D(-30.5, -15.5), +1, tval); // Input contact.
    //Contact ctb = Contact(Vector2D(30.5, 15.5), Vector2D(30.5, -15.5), -1, tval);   // Output contact.
    //std::vector<Contact> contact = {cta, ctb};
    
    holscat = 0.083; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    UsadelSystem usys(mesh, holscat, holabso, tval);
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    usys.initConstant();
    usys.solveNewton(maxit, nsub, tolp, tolr, verbose);
    
    const char* filename = "out/test/test_field_x.csv";
    std::cout << TAG_INFO << "Saving data to file '" << filename << "'...\n";
    usys.saveField(filename, ", ", 16);
    
    return 0;
}

/**
 * Compare the solution of the program to the exact solution known in the waveguide geometry (without absorption, complete channel control).
 */
int testWaveguideSolution() {
    std::cout << "====== TEST WAVEGUIDE SOLUTION ======\n";
    
    int length, width, maxit, nsub, verbose, found;
    double dscat, dabso, tval, tolp, tolr;
    
    length = 20;  // Length of the waveguide (in units of the lattice step).
    width = 10;   // Width of the waveguide (in units of the lattice step).
    dscat = 5.;   // Scattering thickness L/lscat (in the direction of "length").
    dabso = 0.;   // Absorption thickness L/labso.
    tval = 0.5;   // Transmission eigenvalue.
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    // Solves the Usadel equation:
    UsadelSystem usys(length, width, dscat, dabso, tval);  // Create the reference Usadel System with waveguide geometry.
    usys.initConstant();  // Initialize the Usadel solver (using contact guess).
    found = usys.solveNewton(maxit, nsub, tolp, tolr, verbose); // Solves the Usadel equation using the Newton method.
    
    // Compare with the reference solution in the waveguide geometry:
    ComplexVector field_expc(2*usys.getNPoint());
    
    // TODO: Compute the reference solution.................................
    
    
    
    return 0;
}

/**
 * Checks that we retrieve the bimodal distribution rho(T) = Tavg/(2*T*sqrt(1 - T)) in the diffusive waveguide.
 */
int testRhoBimodal() {
    std::cout << "====== TEST RHO BIMODAL ======\n";
    
    int length, width, found, ntval, maxit, nsub, verbose;
    double tval, rho, rho_expc, dscat, dabso, tavg_expc, holscat, holabso, tolp, tolr;
    
    //SquareMesh mesh;  // Create a square mesh.
    length = 20;  // Length of the waveguide (in units of the lattice step). Must be even number here.
    width  = 10;  // Width of the waveguide (in units of the lattice step). Must be even number here.
    dscat = 5.;      // Scattering thickness.
    dabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    holscat = dscat/length;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    tavg_expc = 1./(1. + dscat/(2.*EXTRAPOLEN)); // Average transmission probability in the diffusive regime.
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    ntval = 20; // Number of points on the curve [T, rho(T)].
    
    UsadelSystem usys(length, width, dscat, dabso, 0.5);  // Create the reference Usadel System with waveguide geometry.
    
    const char* filename_distrib = "out/test/distrib/distrib_x.csv";
    const char* filename_field = "out/test/distrib/distrib_x_field.csv";
    const char* filename_mesh  = "out/test/distrib/distrib_x_mesh.csv";
    const char* sep = ", ";  // Separator used in the CSV file.
    const int prec = 16;  // Desired precision (in number of decimals).
    
    std::ofstream ofs;  // Declare output stream object.
    ofs.open(filename_distrib); // Open the file in write mode.
    const auto default_precision = ofs.precision(); // Saves the default precision.
    ofs << std::setprecision(prec); // Set the printing precision.
    
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    ofs << "%% Parameters: h/lscat = " << holscat << ", h/labso = " << holabso
        << ", L/lscat = " << dscat << ", Tavg_expc = " << tavg_expc << "\n"
        << "tval, rho, rho_expc, found\n";
    
    for (int itval = 0; itval < ntval; itval++) {// Loop on the transmission eigenvalue.
        tval = (1. - std::cos((itval + 0.5)*PI/ntval))/2.;   // Choose Chebyshev nodes as the transmission eigenvalues (they are denser at the edges).
        usys.setTransmission(tval); // Assigns the transmission eigenvalue to all contact interactions.
        usys.initConstant();  // Initialize the Usadel solver (using contact guess).
        found = usys.solveNewton(maxit, nsub, tolp, tolr, verbose); // Solves the Usadel equation using the Newton method.
        rho = usys.getRho();  // Actual distribution given by the program.
        rho_expc = tavg_expc/(2.*tval*std::sqrt(1. - tval));   // Expected distribution according to the bimodal law.
        ofs << tval << sep << rho << sep << rho_expc << sep << found << "\n";
        std::cout << TAG_INFO << "# " << itval << " | tval = " << tval 
                  << ", " << (found == 1 ? "found" : "failed") << ", rho = " << rho << ", rho_expc = " << rho_expc << "\n";
    }
    ofs.close();
    ofs << std::setprecision(default_precision); // Restore default precision for printing.
    std::cout << TAG_INFO << "Distribution saved to file: '" << filename_distrib << "'.\n";
    
    usys.saveField(filename_field, sep, prec);
    std::cout << TAG_INFO << "Field saved to file: '" << filename_field << "'.\n";
    
    usys.saveMesh(filename_mesh, sep);
    std::cout << TAG_INFO << "Mesh saved to file: '" << filename_mesh << "'.\n";
    
    return 0;
}

/**
 * Main function of the test of the UsadelSystem object.
 */
int main(int argc, char** argv) {
    
    std::cout << "This is " << PROGRAM_COPYRIGHT << "\n";
    
    //testResidual();
    //testJacobian();
    //testSaveField();
    //testSolve();
    //testWaveguideSolution();
    
    testRhoBimodal();
    
    return 0;
}
