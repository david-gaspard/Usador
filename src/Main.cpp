/****
 * @date Created on 2025-07-18 at 13:44:59 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the main functions of the Usador program.
 * These essentially high-level functions constructing the actual system and solving for specific quantities (Q field, intensities, distribution rho(T), etc.).
 ***/
#include "UsadelSystem.hpp"
#include "BaseTools.hpp"

/******************************************************************************
 * SYSTEM CREATION FUNCTIONS
 *****************************************************************************/

/**
 * Create the UsadelSystem for a simple waveguide geometry.
 */
UsadelSystem createWaveguide() {
    
    int length, width;
    double dscat, dabso, tval;
    const std::string name("waveguide");
    
    length = 60;  // Length of the waveguide in units of the lattice step.
    width = 30;   // Width of the waveguide in units of the lattice step.
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the length and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the length and labso the ballistic absorption length.
    tval = 0.99;  // Transmission eigenvalue, between 0 and 1.
    
    return UsadelSystem(name, length, width, dscat, dabso, tval);  // Create the reference Usadel System with waveguide geometry.
}

/**
 * Create the UsadelSystem for an asymmetric waveguide using open boundary conditions on one side to break the transverse symmetry.
 */
UsadelSystem createAsymmetricWaveguide1() {
    
    double holscat, holabso, tval;
    const std::string name("asymmetric-waveguide-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, 15, 15, NORTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;      // Transmission probability.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for an asymmetric waveguide using open boundary conditions on one side to break the transverse symmetry.
 */
UsadelSystem createAsymmetricWaveguide2() {
    
    double holscat, holabso, tval;
    const std::string name("asymmetric-waveguide-2");
    
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    mesh.addDisk(0, 15, 22);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-22, 22, 30, 40, NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, EAST,  BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, WEST,  BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;      // Transmission probability.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem with open sides.
 */
UsadelSystem createWaveguideOpenSides() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-open-sides");
    
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, +15, +15, NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-10, 10, -15, -15, SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.3;      // Transmission probability.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a circular disordered cavity.
 */
UsadelSystem createCircularCavity() {
    
    double holscat, holabso, tval;
    const std::string name("circular-cavity");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 50);
    mesh.addRectangle(-52, -40, -15, 15);
    mesh.addRectangle(-15, 15, 40, 52);
    mesh.addRectangle(40, 52, -15, 15);
    
    mesh.setBoundaryRegion(-52, -52, -15, 15, WEST, BND_INPUT);
    mesh.setBoundaryRegion(-15, 15, 52, 52, NORTH, BND_OUTPUT);
    mesh.setBoundaryRegion(52, 52, -15, 15, EAST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.9;     // Transmission probability. It is amazing that Tmax=0.9 although the cavity is open (2/3 of channels are controlled).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a cavity in the shape of the Eiffel Tower (just for fun).
 */
UsadelSystem createEiffelTower() {
    
    double holscat, holabso, tval;
    const std::string name("eiffel-tower");
    
    SquareMesh mesh;
    mesh.addPolygon("shape/eiffel-tower.csv", 50);
    
    mesh.setBoundaryRegion(-52, -20, 0, 0, SOUTH, BND_INPUT);
    mesh.setBoundaryRegion(0, 10, 195, 205, EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(20, 52, 0, 1, SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.65;    // Transmission probability.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/******************************************************************************
 * COMPUTATIONAL FUNCTIONS
 *****************************************************************************/

/**
 * Save the mesh in a file and call an external script to plot it.
 * This function is useful to test a mesh before running a simulation.
 */
void plotMesh(UsadelSystem& usys) {
    
    std::stringstream path;
    path << "out/" << usys.getName() << "/holabso_" << usys.getHolabso() << "/mesh/mesh_";
    std::string filename_mesh;
    
    uniqueFilename(path.str(), ".csv", filename_mesh);  // Create a unique filename. The result is of the form "<path><number><suffix>".
    std::cout << TAG_INFO << "Saving mesh to file: '" << filename_mesh << "'...\n";
    usys.saveMesh(filename_mesh, ", ");
    
    std::string cmd("plot/plot_mesh.py " + filename_mesh);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}

/**
 * Compute the fields for a specific value of the transmission eigenvalue "tval".
 * Note that here the "fields" include:
 *    (1) the standard angular parameters (theta, eta), 
 *    (2) the components of the Q field, and 
 *    (3) the normalized intensity profiles I_a, I_b, C_ab.
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 */
void computeFields(UsadelSystem& usys) {
    
    int maxit, nsub, verbose;
    double tolp, tolr;
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing fields from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Tval=" << usys.getTransmission() << ", maxit=" << maxit << ".\n";
    
    usys.initConstant();                                 // Initialize the UsadelSystem using the currently best ansatz (constants).
    usys.solveNewton(maxit, nsub, tolp, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
    
    // Save the data and plot:
    std::stringstream path;
    path << "out/" << usys.getName() << "/holabso_" << usys.getHolabso() << "/tval_" << usys.getTransmission() << "/result_";
    usys.savePlot(path.str());
}

/**
 * Compute the transmission eigenvalue distribution by solving the Usadel equation for different values of the transmission eigenvalue "tval".
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 */
void computeDistribution(UsadelSystem& usys) {
    
    int ntval, maxit, nsub, verbose;
    double tmin, tmax, tolp, tolr;
    
    ntval = 128;  // Number of samples for the transmission eigenvalue.
    tmin = 0.;    // Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    tmax = 1.;    // Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    
    maxit = 200;  // Maximum number of iterations. Typically: 100-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing rho(T) from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", using NT=" << ntval << " on Trange=" << tmin << ":" << tmax << ", maxit=" << maxit << ".\n";
    
    std::cout << TAG_WARN << "TODO: To be implemented...................................\n";
    
    // TODO: Compute the distribution (possibly using OpenMP ?) and save the results [tval, rho, niter] to a temporary array 'data'.
    //double* data = new double[3*ntval]();
    
    std::stringstream path;
    path << "out/" << usys.getName() << "/holabso_" << usys.getHolabso() << "/distrib/result_";
    std::string filename_distrib;
    uniqueFilename(path.str(), ".csv", filename_distrib);  // Create a unique filename. The result is of the form "<path><number><suffix>".
    
    // TODO: Save the array to file 'filename_distrib'.
    
    // TODO: Call the future script "plot_distrib.py" for plotting....
    
}

/**
 * Main function of the Usador program.
 */
int main(int argc, char** argv) {
    
    std::cout << "****** This is " << PROGRAM_COPYRIGHT << " ******\n";
    
    UsadelSystem usys = createWaveguide();
    //UsadelSystem usys = createAsymmetricWaveguide1();
    //UsadelSystem usys = createAsymmetricWaveguide2();
    //UsadelSystem usys = createWaveguideOpenSides();
    //UsadelSystem usys = createCircularCavity();
    //UsadelSystem usys = createEiffelTower();
    
    //plotMesh(usys);  // Plot the mesh only to check it is as expected.
    
    computeFields(usys); // Compute the fields (theta, eta, and Q) and the intensity profile for the given transmission eigenvalue.
    
    //computeDistribution(usys); // Compute the transmission eigenvalue distribution rho(T) by scanning in T.
    
    return 0;
}
