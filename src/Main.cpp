/****
 * @date Created on 2025-07-18 at 13:44:59 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the main functions of the Usador program.
 * These essentially high-level functions constructing the actual system and solving for specific quantities (Q field, intensities, distribution rho(T), etc.).
 ***/
#include "UsadelSystem.hpp"

/**
 * Create the UsadelSystem for a simple waveguide geometry.
 */
UsadelSystem createWaveguide() {
    
    int length, width;
    double dscat, dabso, tval;
    
    length = 60;  // Length of the waveguide in units of the lattice step.
    width = 30;   // Width of the waveguide in units of the lattice step.
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the length and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the length and labso the ballistic absorption length.
    tval = 0.5;   // Transmission eigenvalue, between 0 and 1.
    
    return UsadelSystem(length, width, dscat, dabso, tval);  // Create the reference Usadel System with waveguide geometry.
}

/**
 * Create the UsadelSystem for an asymmetric waveguide.
 */
UsadelSystem createAsymmetricWaveguide() {
    
    double holscat, holabso, tval;
    
    // Create a square mesh:
    SquareMesh mesh;
    mesh.addRectangle(-30, 30, -15, 15);
    mesh.addDisk(0, 15, 22);
    //mesh.removeDisk(0, -15, 12);
    
    mesh.setBoundaryRegion(-22, 22, 30, 40, NORTH, BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, EAST,  BND_OPEN);
    mesh.setBoundaryRegion(-22, 22, 30, 40, WEST,  BND_OPEN);
    //mesh.setBoundaryRegion(-10, 10, -15, -15, SOUTH, BND_OPEN);
    
    mesh.setBoundaryRegion(-30, -30, -15, 15, WEST, BND_INPUT);
    mesh.setBoundaryRegion( 30,  30, -15, 15, EAST, BND_OUTPUT);
    
    mesh.fixNeighbors();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5; // Transmission probability.
    
    return UsadelSystem(mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a circular disordered cavity.
 */
UsadelSystem createCircularCavity() {
    
    double holscat, holabso, tval;
    SquareMesh mesh;
    
    // TODO..........................
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability.
    
    return UsadelSystem(mesh, holscat, holabso, tval);
}

/**
 * Create the UsadelSystem for a cavity in the shape of the Eiffel Tower (just for fun).
 */
UsadelSystem createEiffelTower() {
    
    double holscat, holabso, tval;
    SquareMesh mesh;
    
    // TODO..........................
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability.
    
    return UsadelSystem(mesh, holscat, holabso, tval);
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
    
    std::cout << TAG_INFO << "UsadelSystem with Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Tval=" << usys.getTransmission() << ".\n";
    
    // TODO..................
    
}

/**
 * Compute the transmission eigenvalue distribution by solving the Usadel equation for different values of the transmission eigenvalue "tval".
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 */
void computeDistribution(UsadelSystem& usys) {
    
    // TODO...........
    
}

/**
 * Main function of the Usador program.
 */
int main(int argc, char** argv) {
    
    std::cout << "This is " << PROGRAM_COPYRIGHT << "\n";
    
    UsadelSystem usys = createWaveguide();
    //UsadelSystem usys = createAsymmetricWaveguide();
    //UsadelSystem usys = createCircularCavity();
    //UsadelSystem usys = createEiffelTower();
    
    computeFields(usys);
    
    //computeDistribution(usys);
    
    return 0;
}
