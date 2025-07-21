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
    mesh.addRectangle(-30, 30, -15, 15, BND_MIRROR);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, 15, 15, DIR_NORTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.6791;   // Transmission probability (Tmax=0.679).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for an asymmetric waveguide using open boundary conditions on one side to break the transverse symmetry.
 */
UsadelSystem createAsymmetricWaveguide2() {
    
    double holscat, holabso, tval;
    const std::string name("asymmetric-waveguide-2");
    
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15, BND_MIRROR);
    mesh.addDisk(0, 15, 22, BND_MIRROR);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-22, 22, 30, 40, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, DIR_EAST,  BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, DIR_WEST,  BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.73;     // Transmission probability (Tmax=0.679).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem with open sides.
 */
UsadelSystem createWaveguideOpenSides1() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-open-sides-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15, BND_MIRROR);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.33;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem with open sides.
 */
UsadelSystem createWaveguideOpenSides2() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-open-sides-2");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    mesh.setBoundaryRegion(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-65, -55,  50,  50, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-5,  5, -50, -50, DIR_SOUTH, BND_OPEN);
    mesh.setBoundaryRegion( 55,  65,  50,  50, DIR_NORTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.62;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a waveguide with absorbers.
 */
UsadelSystem createWaveguideAbsorbers1() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-absorbers-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    mesh.removeDisk(-50,  30, 2);
    mesh.setBoundaryRegion(-60, -40, 20, 40, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 50, -30, 2);
    mesh.setBoundaryRegion(40, 60, -40, -20, DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRegion(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a waveguide with absorbers.
 */
UsadelSystem createWaveguideAbsorbers2() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-absorbers-2");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    mesh.removeDisk(-60,  30, 1.5);
    mesh.setBoundaryRegion(-70, -50, 20, 40, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 0, -30, 1.5);
    mesh.setBoundaryRegion(-10, 10, -40, -20, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 60,  30, 1.5);
    mesh.setBoundaryRegion(50, 70, 20, 40, DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRegion(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRegion(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a circular disordered cavity.
 */
UsadelSystem createCircularCavity() {
    
    double holscat, holabso, tval;
    const std::string name("circular-cavity");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 50, BND_MIRROR);
    mesh.addRectangle(-52, -40, -15, 15, BND_MIRROR);
    mesh.addRectangle(-15, 15, 40, 52, BND_MIRROR);
    mesh.addRectangle(40, 52, -15, 15, BND_MIRROR);
    
    mesh.setBoundaryRegion(-52, -52, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion(-15, 15, 52, 52, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRegion(52, 52, -15, 15, DIR_EAST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.9;     // Transmission probability. It is amazing that Tmax=0.9 although the cavity is open (only 2/3 of channels are controlled).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a bigger circular cavity.
 */
UsadelSystem createCircularCavityBig() {
    
    double holscat, holabso, tval;
    const std::string name("circular-cavity-big");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 100, BND_MIRROR);
    mesh.addRectangle(-105, -80, -30, 30, BND_MIRROR);
    mesh.addRectangle(-30, 30, 80, 105, BND_MIRROR);
    mesh.addRectangle(80, 105, -30, 30, BND_MIRROR);
    
    mesh.setBoundaryRegion(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRegion(105, 105, -30, 30, DIR_EAST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.8;     // Transmission probability. It is amazing that Tmax=0.9 although the cavity is open (only 2/3 of channels are controlled).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a bigger circular cavity.
 */
UsadelSystem createCircularCavityHole() {
    
    double holscat, holabso, tval;
    const std::string name("circular-cavity-hole");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 100, BND_MIRROR);
    mesh.addRectangle(-105, -80, -30, 30, BND_MIRROR);
    mesh.addRectangle(-30, 30, 80, 105, BND_MIRROR);
    mesh.addRectangle(80, 105, -30, 30, BND_MIRROR);
    
    mesh.removeDisk(-49, 49, 5);
    mesh.setBoundaryRegion(-60, -40, 40, 60, DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRegion(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRegion(105, 105, -30, 30, DIR_EAST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.3;     // Transmission probability. It is amazing that Tmax=0.9 although the cavity is open (only 2/3 of channels are controlled).
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for an alternative circular cavity.
 */
UsadelSystem createCircularCavityAlt() {
    
    double holscat, holabso, tval;
    const std::string name("circular-cavity-alt");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 100, BND_MIRROR);
    mesh.addRectangle(-105, -80, -30, 30, BND_MIRROR);
    mesh.addRectangle(-30, 30, 80, 105, BND_MIRROR);
    //mesh.addRectangle(80, 105, -30, 30, BND_MIRROR);
    
    mesh.setBoundaryRegion(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRegion(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRegion(-80, -50, 50, 80, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-80, -50, 50, 80, DIR_WEST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.70;    // Transmission probability. Tmax=0.75 but only Tmax=0.70 converges in a first Newton-Raphson pass.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a cavity in the shape of the Eiffel Tower (just for fun).
 */
UsadelSystem createEiffelTower() {
    
    double holscat, holabso, tval;
    const std::string name("eiffel-tower");
    
    SquareMesh mesh;
    mesh.addPolygon("shape/eiffel-tower.csv", 50, BND_MIRROR);
    mesh.removeDisk(0, 50, 5);
    
    mesh.setBoundaryRegion(-52, -20, 0, 0, DIR_SOUTH, BND_INPUT);
    mesh.setBoundaryRegion(20, 52, 0, 1, DIR_SOUTH, BND_INPUT);
    mesh.setBoundaryRegion(0, 10, 195, 205, DIR_EAST, BND_OUTPUT);
    
    mesh.setBoundaryRegion(-11, 11, 39, 61, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-11, 11, 39, 61, DIR_SOUTH, BND_OPEN);
    mesh.setBoundaryRegion(-11, 11, 39, 61, DIR_WEST, BND_OPEN);
    mesh.setBoundaryRegion(-11, 11, 39, 61, DIR_EAST, BND_OPEN);
    
    mesh.fixNeighbors();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.1;    // Transmission probability. Without absorber: Tmax=0.68.
    
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
 *    (1) the points (x, y) and the boundary conditions (north, south, east, west),
 *    (2) the standard angular parameters (theta, eta) with their real and imaginary parts, 
 *    (3) the components of the Q field (real and imaginary parts), and 
 *    (4) the normalized intensity profiles I_a, I_b, C_ab (positive real quantities).
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 */
void computeFields(UsadelSystem& usys) {
    
    int maxit, nsub, verbose;
    double tolp, tolr;
    
    maxit = 50;  // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing fields from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Tval=" << usys.getTransmission() << ", maxit=" << maxit << ".\n";
    
    //usys.setTransmission(0.70);
    usys.initConstant();                                 // Initialize the UsadelSystem using the currently best ansatz (constants).
    usys.solveNewton(maxit, nsub, tolp, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
    //usys.setTransmission(0.75);
    //usys.solveNewton(maxit, nsub, tolp, tolr, verbose);
    
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
    
    maxit = 200;  // Maximum number of iterations. Typically: 50-500.
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
    
    //UsadelSystem usys = createWaveguide();
    //UsadelSystem usys = createAsymmetricWaveguide1();
    //UsadelSystem usys = createAsymmetricWaveguide2();
    //UsadelSystem usys = createWaveguideOpenSides1();
    //UsadelSystem usys = createWaveguideOpenSides2();
    //UsadelSystem usys = createWaveguideAbsorbers1();
    //UsadelSystem usys = createWaveguideAbsorbers2();
    //UsadelSystem usys = createCircularCavity();
    //UsadelSystem usys = createCircularCavityBig();
    //UsadelSystem usys = createCircularCavityAlt();
    //UsadelSystem usys = createCircularCavityHole();
    UsadelSystem usys = createEiffelTower();
    
    plotMesh(usys);  // Plot the mesh only to check it is as expected.
    
    //computeFields(usys); // Compute the fields (theta, eta, and Q) and the intensity profile for the given transmission eigenvalue.
    
    //computeDistribution(usys); // Compute the transmission eigenvalue distribution rho(T) by scanning in T.
    
    return 0;
}
