/****
 * @date Created on 2025-07-18 at 13:44:59 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the main functions of the Usador program.
 * These essentially high-level functions constructing the actual system and solving for specific quantities (Q field, intensities, distribution rho(T), etc.).
 ***/
#include "UsadelSystem.hpp"
#include "BaseTools.hpp"
#include <iomanip>

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
    dabso = 0.1;   // Absorption thickness, L/labso, where L is the length and labso the ballistic absorption length.
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
    
    mesh.setBoundaryRectangle(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-10, 10, 15, 15, DIR_NORTH, BND_OPEN);
    
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
    
    mesh.setBoundaryRectangle(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-22, 22, 30, 40, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-22, 22, 30, 40, DIR_EAST,  BND_OPEN);
    mesh.setBoundaryRectangle(-22, 22, 30, 40, DIR_WEST,  BND_OPEN);
    
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
    
    mesh.setBoundaryRectangle(-30, -30, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 30,  30, -15, 15, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
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
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-65, -55,  50,  50, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-5,  5, -50, -50, DIR_SOUTH, BND_OPEN);
    mesh.setBoundaryRectangle( 55,  65,  50,  50, DIR_NORTH, BND_OPEN);
    
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
    
    mesh.removeDisk(-50, 30, 2.);
    mesh.setBoundaryDisk(-50, 30, 3., DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 50, -30, 2.);
    mesh.setBoundaryDisk(50, -30, 3., DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
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
    
    mesh.removeDisk(-60, 30, 1.5);
    mesh.setBoundaryDisk(-60, 30, 2.5, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 0, -30, 1.5);
    mesh.setBoundaryDisk(0, -30, 2.5, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk( 60,  30, 1.5);
    mesh.setBoundaryDisk(60, 30, 2.5, DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle(-10, 10, +15, +15, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-10, 10, -15, -15, DIR_SOUTH, BND_OPEN);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
UsadelSystem createDoubleWaveguide1() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    mesh.removeRectangle(-60, 60, -29, 31);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.2/200; // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.65;        // Transmission probability close to Tmax=0.27.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
UsadelSystem createDoubleWaveguide2() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-2");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    mesh.removeRectangle(-60, 60, -29, 31);
    mesh.setBoundaryRectangle(-61, -61, -45, 45, DIR_EAST, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.50;     // Transmission probability close to Tmax=0.27.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide.
 */
UsadelSystem createWaveguideObstacle1() {
    
    double holscat, holabso, tval;
    const std::string name("waveguide-obstacle-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    
    //mesh.removeDisk(0, 0, 29);
    //mesh.setBoundaryRectangle(1, 5, 28, 30, DIR_SOUTH, BND_OPEN);
    
    mesh.removeDisk(0, 0, 3.5);
    mesh.removeDisk(0, 33, 3.5);
    mesh.removeDisk(0, -33, 3.5);
    mesh.setBoundaryRectangle(-10, 10, -40, 40, DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -22, -11, DIR_WEST, BND_INPUT);
    //mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    mesh.setBoundaryRectangle( 100,  100, 11, 22, DIR_EAST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.25;     // Transmission probability close to Tmax=0.27.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double circular waveguide.
 */
UsadelSystem createCircularDoubleWaveguide1() {
    
    double holscat, holabso, tval;
    const std::string name("circular-double-waveguide-1");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -30, 30, BND_MIRROR);
    mesh.addDisk(0, 0, 90, BND_MIRROR);
    mesh.removeDisk(0, 0, 50);
    
    //mesh.removeDisk(0, 0, 10);
    //mesh.setBoundaryDisk(0, 0, 11., DIR_ALL, BND_OPEN);
    
    mesh.removeDisk(0, 70, 2.5);
    mesh.setBoundaryDisk(0, 70, 3.5, DIR_ALL, BND_OPEN);
    
    mesh.removeDisk(0, -70, 1.5);
    mesh.setBoundaryDisk(0, -70, 2.5, DIR_ALL, BND_OPEN);
    
    mesh.addRectangle(-60, 60, -2, 2, BND_MIRROR);
    
    //mesh.addRectangle(-5, 5, 85, 100, BND_MIRROR);
    //mesh.setBoundaryRectangle(-5, 5, 100, 100, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -30, 30, DIR_EAST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax=0.27.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a  circular cavity with opposite leads.
 */
UsadelSystem createCircularOpposite() {
    
    double holscat, holabso, tval;
    const std::string name("circular-opposite");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 100, BND_MIRROR);
    mesh.addRectangle(-110, -70, 5, 65, BND_MIRROR);
    mesh.addRectangle(-110, -70, -65, -5, BND_MIRROR);
    
    mesh.addRectangle(90, 110, -30, 30, BND_MIRROR);
    mesh.setBoundaryRectangle(110, 110, -30, 30, DIR_EAST, BND_INPUT);
    
    mesh.setBoundaryRectangle(-110, -110, 5, 65, DIR_WEST, BND_OPEN);
    mesh.setBoundaryRectangle(-110, -110, -65, -5, DIR_WEST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.8;     // Transmission probability close to Tmax=0.27.
    
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
    
    mesh.setBoundaryRectangle(-52, -52, -15, 15, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(-15, 15, 52, 52, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRectangle(52, 52, -15, 15, DIR_EAST, BND_OPEN);
    
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
    
    mesh.setBoundaryRectangle(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRectangle(105, 105, -30, 30, DIR_EAST, BND_OPEN);
    
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
    mesh.setBoundaryDisk(-49, 49, 6., DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRectangle(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRectangle(105, 105, -30, 30, DIR_EAST, BND_OPEN);
    
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
    
    mesh.setBoundaryRectangle(-105, -105, -30, 30, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(-30, 30, 105, 105, DIR_NORTH, BND_OUTPUT);
    mesh.setBoundaryRectangle(-80, -50, 50, 80, DIR_NORTH, BND_OPEN);
    mesh.setBoundaryRectangle(-80, -50, 50, 80, DIR_WEST, BND_OPEN);
    
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
    mesh.setBoundaryDisk(0, 50, 6., DIR_ALL, BND_OPEN);
    
    mesh.setBoundaryRectangle(-52, -20, 0, 0, DIR_SOUTH, BND_INPUT);
    mesh.setBoundaryRectangle(20, 52, 0, 1, DIR_SOUTH, BND_INPUT);
    mesh.setBoundaryRectangle(0, 10, 195, 205, DIR_EAST, BND_OUTPUT);
    
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
    path << "out/" << usys.getName() << "/mesh/mesh_";
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
    
    usys.initConstant();                                 // Initialize the UsadelSystem using the currently best ansatz (constants).
    usys.solveNewton(maxit, nsub, tolp, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
    //usys.setTransmission(0.98);
    //usys.solveNewton(maxit, nsub, tolp, tolr, verbose);
    
    // Save the data and plot:
    std::stringstream path;
    path << "out/" << usys.getName() << "/holabso_" << usys.getHolabso() << "/tval_" << usys.getTransmission() << "/result_";
    usys.savePlot(path.str());
}

/**
 * Save the distribution rho(T) to a file and then call an external script to plot it.
 */
void plotDistribution(const UsadelSystem& usys, const double* rhodata, const int ntval, const double tmin, const double tmax, const int maxit, const double elapsed_time) {
    
    const char* sep = ", ";  // Separator used in the output file.
    const int prec = 16;     // Precision used in the output file.
    
    std::stringstream path;
    path << "out/" << usys.getName() << "/holabso_" << usys.getHolabso() << "/distrib/result_";
    std::string filename_distrib;
    uniqueFilename(path.str(), ".csv", filename_distrib);  // Create a unique filename. The result is of the form "<path><number><suffix>".
    
    std::cout << TAG_INFO << "Saving distribution to file: '" << filename_distrib << "'...\n";
    std::ofstream ofs;  // Declare output stream object.
    ofs.open(filename_distrib.c_str()); // Open the file in write mode.
    ofs << std::setprecision(prec); // Set the printing precision.
    
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    ofs << "%% Parameters: name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() 
        << ", h/labso=" << usys.getHolabso() << ", Trange=" << tmin << ":" << tmax << ", NT=" << ntval << ", maxit=" << maxit 
        << ", computation_time=" << elapsed_time << "s.\n"
        << "tval" << sep << "rho" << sep << "niter" << sep << "converged\n";
    
    for (int i = 0; i < ntval; i++) {// Loop on the computed couples [T, rho(T)].
        ofs << rhodata[3*i] << sep << rhodata[3*i+1] << sep << rhodata[3*i+2] << sep << (rhodata[3*i+2] <= maxit ? "success" : "failure") << "\n";
    }
    
    ofs.close();  // Close the file.
    
    std::string cmd("plot/plot_distrib.py " + filename_distrib);
    std::cout << TAG_EXEC << cmd << "\n";
    if (std::system(cmd.c_str())) {
        std::cout << TAG_WARN << "The plot script returned an error.\n";
    }
}

/**
 * Compute the transmission eigenvalue distribution by solving the Usadel equation for different values of the transmission eigenvalue "tval".
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 * This serial version can be useful when many simulations fail because each iteration uses the previous solution as a first ansatz.
 */
void computeDistributionSerial(UsadelSystem& usys) {
    
    int ntval, maxit, nsub, verbose, niter;
    double tmin, tmax, tval, rho, tolp, tolr;
    
    ntval = 32;  // Number of samples for the transmission eigenvalue.
    tmin = 0.;    // Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    tmax = 1.;    // Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    
    maxit = 50;   // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing rho(T) from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", NT=" << ntval << ", Trange=" << tmin << ":" << tmax << ", maxit=" << maxit << ".\n";
    
    auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    double* rhodata = new double[3*ntval];  // Temporary array for storing the numerical results [T, rho(T), niter].
    
    usys.setTransmission(tmin + (tmax-tmin) * (1. - std::cos(0.5*PI/ntval))/2.);   // Assigns the 'tmin' before initializing.
    usys.initConstant(); // Initialize the UsadelSystem using the currently best ansatz (constants).
    
    for (int i = 0; i < ntval; i++) {// Loop over the points of the mesh.
        
        tval = tmin + (tmax-tmin) * (1. - std::cos((i + 0.5)*PI/ntval))/2.;  // Choose Chebyshev nodes as the transmission eigenvalues (they are denser at the edges).
        usys.setTransmission(tval);   // Assigns the transmission eigenvalue.
        niter = usys.solveNewton(maxit, nsub, tolp, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
        rho = usys.getRho();  // Compute the transmission eigenvalue distribution at the given 'tval'.
        
        rhodata[3*i] = tval;  // Save the data to the array 'rhodata'.
        rhodata[3*i+1] = rho;
        rhodata[3*i+2] = niter;
        
        std::cout << TAG_INFO << "# " << (i+1) << "\t| tval=" << tval 
                  << ",\trho=" << rho << ",\tniter=" << niter << ",\t" << (niter <= maxit ? "success" : "failure") << "\n";
    }
    
    double elapsed_time = endProgressBar(start);  // Finalize the progress bar and gets the total time (in seconds).
    
    plotDistribution(usys, rhodata, ntval, tmin, tmax, maxit, elapsed_time);
    
    delete[] rhodata;
}


/**
 * Compute the transmission eigenvalue distribution by solving the Usadel equation for different values of the transmission eigenvalue "tval".
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 * This is a parallelized version using OpenMP.
 */
void computeDistributionOMP(UsadelSystem& usys) {
    
    int i, ntval, maxit, nsub, verbose, niter, nthread;
    double tmin, tmax, tval, rho, tolp, tolr;
    
    ntval = 128;  // Number of samples for the transmission eigenvalue.
    tmin = 0.;    // Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    tmax = 1.;    // Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    nthread = 10; // Number of execution threads for OpenMP (typically the number of CPU cores).
    
    maxit = 50;   // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    tolp = 1e-7;  // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing rho(T) from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() 
              << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", NT=" << ntval << ", Trange=" << tmin << ":" << tmax 
              << ", maxit=" << maxit << ", nthread=" << nthread << ".\n";
    
    const std::string msg = "distrib, " + std::to_string(nthread) + " thr";
    int cjob = 0;  // Initialize the number of completed jobs.
    auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    double* rhodata = new double[3*ntval];  // Temporary array for storing the numerical results [T, rho(T), niter].
    
    UsadelSystem* usys_loc;
    
    #pragma omp parallel private(usys_loc, i, tval, niter, rho) num_threads(nthread)
    {
        
        usys_loc = new UsadelSystem(usys); // Deep copy of the UsadelSystem on each thread (OMP-Private). 
                                           // Call the deep copy constructor (explicitly defined elsewhere).
        
        #pragma omp for schedule(dynamic,1)
        for (i = 0; i < ntval; i++) {// Loop over the points of the mesh.
            
            tval = tmin + (tmax-tmin) * (1. - std::cos((i + 0.5)*PI/ntval))/2.;  // Choose Chebyshev nodes as the transmission eigenvalues (they are denser at the edges).
            
            usys_loc->setTransmission(tval);   // Assigns the transmission eigenvalue.
            usys_loc->initConstant();  // Initialize the UsadelSystem using the currently best ansatz (constants).
                                       // WARNING: It may be more appropriate to use the previous solution as an ansatz.
            
            niter = usys_loc->solveNewton(maxit, nsub, tolp, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
            rho = usys_loc->getRho();  // Compute the transmission eigenvalue distribution at the given 'tval'.
            
            rhodata[3*i] = tval;  // Save the data to the array 'rhodata'.
            rhodata[3*i+1] = rho;
            rhodata[3*i+2] = niter;
            
            #pragma omp critical
            {
                cjob++;
                printProgressBar(cjob, ntval, msg, start); // Print the progress bar. Must be critical since 'cjob' is a OMP-Shared variable.
                //std::cout << TAG_INFO << "# " << (i+1) << "\t| tval=" << tval 
                //  << ",\trho=" << rho << ",\tniter=" << niter << ",\t" << (niter <= maxit ? "success" : "failure") << "\n";
            }
            
        }
        
        delete usys_loc;
        
    }
    
    double elapsed_time = endProgressBar(start);  // Finalize the progress bar and gets the total time (in seconds).
    
    plotDistribution(usys, rhodata, ntval, tmin, tmax, maxit, elapsed_time);
    
    delete[] rhodata;
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
    //UsadelSystem usys = createDoubleWaveguide1();
    //UsadelSystem usys = createDoubleWaveguide2();
    //UsadelSystem usys = createWaveguideObstacle1();
    //UsadelSystem usys = createCircularDoubleWaveguide1();
    //UsadelSystem usys = createCircularOpposite();
    //UsadelSystem usys = createCircularCavity();
    //UsadelSystem usys = createCircularCavityBig();
    //UsadelSystem usys = createCircularCavityAlt();
    //UsadelSystem usys = createCircularCavityHole();
    UsadelSystem usys = createEiffelTower();
    
    plotMesh(usys);  // Plot the mesh only to check it is as expected.
    
    //computeFields(usys); // Compute the fields (theta, eta, and Q) and the intensity profile for the given transmission eigenvalue.
    
    //computeDistributionSerial(usys); // Compute the transmission eigenvalue distribution rho(T) by scanning in T.
    //computeDistributionOMP(usys); // Compute the transmission eigenvalue distribution rho(T). Parallelized version.
    
    return 0;
}
