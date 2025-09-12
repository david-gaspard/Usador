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
    
    length = 150;  // Length of the waveguide in units of the lattice step.
    width = 150;   // Width of the waveguide in units of the lattice step.
    dscat = 10.;   // Scattering thickness, L/lscat, where L is the length and lscat the scattering mean free path.
    dabso = 0.;    // Absorption thickness, L/labso, where L is the length and labso the ballistic absorption length.
    tval = 0.999;   // Transmission eigenvalue, between 0 and 1. For dscat=5, dabso=0.1, we have Tmax=0.75.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
UsadelSystem createDoubleWaveguide0() {
    
    const std::string name("double-waveguide-0");
    
    const int tlength = 200;          // Number of lattice points in the longitudinal direction, L/h. Default: length=150
    const int twidth  = 200;          // Number of lattice points in the transverse direction, W/h. Default: width=150
    const int xwidth  = tlength/5;    // Width of the corridors near input/output.
    const int ywidth  = twidth/5;     // Width of the corridors of the two waveguides (before the yshift).
    const int yshift  = twidth/50;    // Transverse shift.
    
    const double dscat = 10.;  // Scattering depth, L/lscat.
    const double dabso = 0.;   // Absorption depth, L/labso.
    
    // Construct the mesh:
    SquareMesh mesh;
    mesh.addRectangle(0, tlength, 0, twidth, BND_MIRROR);
    mesh.removeRectangle(xwidth, tlength-xwidth, ywidth+yshift, twidth-ywidth+yshift);
    
    // Setup boundary conditions:
    mesh.setBoundaryRectangle(0, 0, 0, twidth, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle(tlength, tlength, 0, twidth, DIR_EAST, BND_OUTPUT);
    //mesh.setBoundaryRectangle(xwidth-1, xwidth-1, ywidth+yshift, twidth-ywidth+yshift, DIR_EAST, BND_OPEN);
    
    mesh.finalize();
    
    // Defines the physical parameters:
    const double holscat = dscat/tlength;
    const double holabso = dabso/tlength;
    const double tval = 0.98;
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.50;     // Transmission probability close to Tmax=0.27.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with a minute asymmetry.
 */
UsadelSystem createDoubleWaveguide3() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-3");
    
    SquareMesh mesh;
    mesh.addDisk(0, 0, 100, BND_MIRROR);
    mesh.removeDisk(0, 0, 78);
    mesh.removeHalfDisk(0, 0, 82);
    
    mesh.addRectangle(-90, -75, -15, 15, BND_MIRROR);
    mesh.setBoundaryRectangle(-75, -75, -15, 15, DIR_EAST, BND_OPEN);
    
    mesh.addRectangle(-110, -85, -20, 20, BND_MIRROR);
    mesh.setBoundaryRectangle(-110, -110, -50, 50, DIR_WEST, BND_INPUT);
    
    mesh.addRectangle(85, 110, -20, 20, BND_MIRROR);
    mesh.setBoundaryRectangle( 110,  110, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./60; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;   // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.25;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with a minute asymmetry due to a lossy constriction.
 */
UsadelSystem createDoubleWaveguide4() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-4");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    //mesh.removeRectangle(-59, 80, -30, 32);
    mesh.removeRectangle(-59, 80, -30, 30);
    mesh.removeRectangle(-59, -55, 30, 35);
    
    mesh.setBoundaryRectangle(-60, -60, -35, 35, DIR_EAST, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.55;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry due to a lossy constriction.
 */
UsadelSystem createDoubleWaveguide5() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-5");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    mesh.removeRectangle(-80, 80, -30, 30);
    mesh.removeRectangle(-80, -75, 30, 35);
    mesh.removeRectangle(-80, -75, 45, 50);
    
    mesh.setBoundaryRectangle(-80, -75, 36, 36, DIR_SOUTH, BND_OPEN);
    mesh.setBoundaryRectangle(-80, -75, 44, 44, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.8;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry due to lossy channels.
 */
UsadelSystem createDoubleWaveguide6() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-6");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    mesh.removeRectangle(-80, 80, -30, 30);
    
    mesh.setBoundaryRectangle(60, 80, 50, 50, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.7;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry due to lossy channels.
 */
UsadelSystem createDoubleWaveguide7() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-7");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 50, BND_MIRROR);
    mesh.removeRectangle(-80, 80, -30, 30);
    
    mesh.addRectangle(40, 60, 10, 30, BND_MIRROR);
    //mesh.setBoundaryRectangle(40, 60, 10, 10, DIR_SOUTH, BND_OPEN);
    mesh.addRectangle(-60, 60, -10, 10, BND_MIRROR);
    mesh.setBoundaryRectangle(-60, -60, -10, 10, DIR_WEST, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.3;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry due to lossy channels.
 */
UsadelSystem createDoubleWaveguide8() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-8");
    
    SquareMesh mesh;
    mesh.addRectangle(-100, 100, -50, 70, BND_MIRROR);
    mesh.removeRectangle(-80, 80, -30, 30);
    
    mesh.addRectangle(40, 60, 10, 30, BND_MIRROR);
    //mesh.addRectangle(20, 60, -10, 10, BND_MIRROR);
    //mesh.setBoundaryRectangle(20, 20, -10, 10, DIR_WEST, BND_OPEN);
    
    //mesh.addRectangle(-60, 60, -10, 10, BND_MIRROR);
    //mesh.setBoundaryRectangle(-60, -60, -10, 10, DIR_WEST, BND_OPEN);
    
    mesh.addRectangle(20, 60, -10, 10, BND_MIRROR);
    mesh.setBoundaryRectangle(20, 60, -10, -10, DIR_SOUTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-100, -100, -70, 70, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 100,  100, -70, 70, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.3;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry with one short lossy channel and one long lossless channel.
 */
UsadelSystem createDoubleWaveguide9() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-9");
    
    SquareMesh mesh;
    mesh.addRectangle(-110, 110, -50, 50, BND_MIRROR);
    mesh.removeRectangle(-90, 90, -30, 30);
    
    // Add maze in lossless channel:
    mesh.addRectangle(-70, -10, -30, 10, BND_MIRROR);
    mesh.removeRectangle(-50, -30, -50, -10);
    
    mesh.addRectangle(10, 70, -30, 10, BND_MIRROR);
    mesh.removeRectangle(30, 50, -50, -10);
    
    // Add losses in short channel:
    //mesh.addRectangle(70, 90, 50, 70, BND_MIRROR);
    //mesh.setBoundaryRectangle(70, 90, 70, 70, DIR_NORTH, BND_OPEN);
    
    //mesh.addRectangle(90, 110, 50, 70, BND_MIRROR);
    //mesh.addRectangle(-110, 110, 70, 90, BND_MIRROR);
    //mesh.setBoundaryRectangle(-110, -110, 70, 90, DIR_WEST, BND_OPEN);
    
    mesh.addRectangle(90, 110, 50, 70, BND_MIRROR);
    mesh.setBoundaryRectangle(90, 110, 70, 70, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-110, -110, -50, 50, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 110,  110, -50, 50, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.97;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry with one short lossy channel and one long lossless channel.
 */
UsadelSystem createDoubleWaveguide10() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-10");
    
    SquareMesh mesh;
    
    // Add maze in lossless channel:
    mesh.addRectangle(-110, 110, -110, 110, BND_MIRROR);
    mesh.removeRectangle( -90,  90,  70,  90);
    mesh.removeRectangle( -90,  90, -10,  10);
    mesh.removeRectangle( -90,  90, -90, -70);
    mesh.removeRectangle( -10,  10, -90,  90);
    mesh.removeRectangle(-110, -30,  30,  50);
    mesh.removeRectangle(  30, 110,  30,  50);
    mesh.removeRectangle(-110, -30, -50, -30);
    mesh.removeRectangle(  30, 110, -50, -30);
    
    // Add losses in short channel:
    mesh.addRectangle(90, 110, 110, 130, BND_MIRROR);
    mesh.setBoundaryRectangle(90, 110, 130, 130, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-110, -110, 50, 110, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 110,  110, 50, 110, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 5./200; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.5;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry with one short lossy channel and one long lossless channel.
 */
UsadelSystem createDoubleWaveguide11() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-11");
    
    SquareMesh mesh;
    
    // Add maze in lossless channel:
    mesh.addRectangle(-110, 110, -70, 70, BND_MIRROR);
    mesh.removeRectangle( -90,  90,  30,  50);
    
    mesh.removeRectangle( -90, -70, -50,  50);
    mesh.removeRectangle( -10,  10, -50,  50);
    mesh.removeRectangle(  70,  90, -50,  50);
    
    mesh.removeRectangle( -50, -30, -70,  10);
    mesh.removeRectangle(  30,  50, -70,  10);
    
    // Add losses in short channel:
    mesh.addRectangle(90, 110, 70, 90, BND_MIRROR);
    mesh.setBoundaryRectangle(90, 110, 90, 90, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-110, -110, -70, 70, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 110,  110, -70, 70, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 0.025; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;    // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.833;     // Transmission probability close to Tmax.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a double waveguide with an asymmetry with one short lossy channel and one long lossless channel.
 * Same as createDoubleWaveguide11() but for h/lscat = 0.1 (lscat = half the tube diameter, the recommended value).
 */
UsadelSystem createDoubleWaveguide12() {
    
    double holscat, holabso, tval;
    const std::string name("double-waveguide-12");
    
    SquareMesh mesh;
    
    // Add maze in lossless channel:
    mesh.addRectangle(-110, 110, -70, 70, BND_MIRROR);
    mesh.removeRectangle( -90,  90,  30,  50);
    
    mesh.removeRectangle( -90, -70, -50,  50);
    mesh.removeRectangle( -10,  10, -50,  50);
    mesh.removeRectangle(  70,  90, -50,  50);
    
    mesh.removeRectangle( -50, -30, -70,  10);
    mesh.removeRectangle(  30,  50, -70,  10);
    
    // Add losses in short channel:
    mesh.addRectangle(90, 110, 70, 90, BND_MIRROR);
    mesh.setBoundaryRectangle(90, 110, 90, 90, DIR_NORTH, BND_OPEN);
    
    mesh.setBoundaryRectangle(-110, -110, -70, 70, DIR_WEST, BND_INPUT);
    mesh.setBoundaryRectangle( 110,  110, -70, 70, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
    holscat = 0.1; // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.9;    // Transmission probability close to Tmax.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();  // Do not forget to fix the neighbors to finalize the mesh.
    
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
    
    mesh.finalize();
    
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
    
    mesh.finalize();
    
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
    
    mesh.finalize();
    
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
    
    mesh.finalize();
    
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
    
    mesh.finalize();
    
    holscat = 0.1;  // Value of h/lscat, where "h" is the lattice step and "lscat" the scattering mean free path.
    holabso = 0.0;  // Value of h/labso, where "h" is the lattice step and "labso" the absorption length.
    tval = 0.2;    // Transmission probability. Without absorber: Tmax=0.68.
    
    return UsadelSystem(name, mesh, holscat, holabso, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabTransmission1() {
    
    int width, thick, size_input, size_output, y_output;
    double dscat, dabso, tval;
    const std::string name("finite-slab-transmission-1");
    
    width = 300;        // Transverse width of the slab. Typically = 300.
    thick = 100;        // Thickness of the slab. Typically = 100.
    size_input  = width/2;  // Diameter of the input beam. Typically: width/2.
    size_output = width/8;   // Diameter of the target (output beam). Typically: width/8.
    y_output = 0;       // Position of the output beam. Typically: 0.
    
    SquareMesh mesh;
    mesh.addRectangle(0, thick, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(0, 0, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(thick, y_output, size_output/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 10.;  // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.4;   // Transmission probability. Tmax=0.3 for size_input=160, size_output=20, dscat=5, dabso=0.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabTransmission2() {
    
    int width, thick, size_input, size_output, y_output;
    double dscat, dabso, tval;
    const std::string name("finite-slab-transmission-2");
    
    width = 300;       // Transverse width of the slab. Typically = 300.
    thick = 100;       // Thickness of the slab. Typically = 100.
    size_input  = 20;  // Diameter of the input beam. Typically: 160.
    size_output = 20;  // Diameter of the target (output beam). Typically: 20.
    y_output = 0;      // Position of the output beam. Typically: 0.
    
    SquareMesh mesh;
    mesh.addRectangle(-thick/2, thick/2, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(-thick/2, 0, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(thick/2, y_output, size_output/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.3;  // Transmission probability. Tmax=0.3 for size_input=160, size_output=20, dscat=5, dabso=0.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabTransmission3() {
    
    int width, thick, size_input, size_output, y_output;
    double dscat, dabso, tval;
    const std::string name("finite-slab-transmission-3");
    
    width = 300;       // Transverse width of the slab. Typically = 300.
    thick = 100;       // Thickness of the slab. Typically = 100.
    size_input  = 160;  // Diameter of the input beam. Typically: 160.
    size_output = 160;  // Diameter of the target (output beam). Typically: 20.
    y_output = 0;      // Position of the output beam. Typically: 0.
    
    SquareMesh mesh;
    mesh.addRectangle(-thick/2, thick/2, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(-thick/2, 0, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(thick/2, y_output, size_output/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.3;  // Transmission probability. Tmax=0.3 for size_input=160, size_output=20, dscat=5, dabso=0.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createSlabTransmission() {
    
    int width, thick, size_input, size_output, y_output;
    double dscat, dabso, tval;
    const std::string name("slab-transmission");
    
    width = 300;       // Transverse width of the slab. Typically = 300.
    thick = 100;       // Thickness of the slab. Typically = 100.
    size_input  = 150;  // Diameter of the input beam. Typically: 150.
    size_output = 150;  // Diameter of the target (output beam). Typically: 20.
    y_output = 0;      // Position of the output beam. Typically: 0.
    
    SquareMesh mesh;
    mesh.addRectangle(0, thick, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(0, 0, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(thick, y_output, size_output/2, DIR_EAST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 10.;  // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.8;   // Transmission probability. Tmax=0.3 for size_input=160, size_output=20, dscat=5, dabso=0.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabRemission1() {
    
    int width, thick, size_input, size_output, separation, y_input;
    double dscat, dabso, tval;
    const std::string name("finite-slab-remission-1");
    
    width = 300;       // Transverse width of the slab.
    thick = 100;       // Thickness of the slab.
    size_input = 15;  // Diameter of the input beam.
    size_output = 15;  // Diameter of the target (output beam).
    separation = 20;   // Distance between the input and otput beams.
    y_input = (size_output + separation)/2;  // Position of the input beam ensuring that the distance with respect to the upper and lower edges are equal.
    
    SquareMesh mesh;
    mesh.addRectangle(-thick/2, thick/2, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(-thick/2, y_input, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-thick/2, y_input - size_input/2 - separation - size_output/2, size_output/2, DIR_WEST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.15;  // Transmission probability.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabRemission2() {
    
    int width, thick, size_input, size_output, separation, y_input;
    double dscat, dabso, tval;
    const std::string name("finite-slab-remission-2");
    
    width = 300;       // Transverse width of the slab.
    thick = 100;       // Thickness of the slab.
    size_input = 20;   // Diameter of the input beam.
    size_output = 20;  // Diameter of the target (output beam).
    separation = 40;   // Distance between the input and otput beams.
    y_input = (size_output + separation)/2;  // Position of the input beam ensuring that the distance with respect to the upper and lower edges are equal.
    
    SquareMesh mesh;
    mesh.addRectangle(-thick/2, thick/2, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(-thick/2, y_input, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-thick/2, y_input - size_input/2 - separation - size_output/2, size_output/2, DIR_WEST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.15;  // Transmission probability.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabRemission3() {
    
    int width, thick, winput, woutput, ysep;
    double dscat, dabso, tval;
    const std::string name("finite-slab-remission-3");
    
    width = 300;        // Transverse width of the slab.
    thick = 100;        // Thickness of the slab.
    winput = width/8;   // Diameter of the input beam.
    woutput = width/8;  // Diameter of the target (output beam).
    ysep = width/3;     // Distance between the input and otput beams.
    
    SquareMesh mesh;
    mesh.addRectangle(0, thick, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(0, -ysep/2,  winput/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(0, +ysep/2, woutput/2, DIR_WEST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 10.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.10;  // Transmission probability.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
}

/**
 * Create a finite slab with a focusing point.
 */
UsadelSystem createFiniteSlabDoubleRemission() {
    
    int width, thick, size_input, size_output, separation, y_input;
    double dscat, dabso, tval;
    const std::string name("finite-slab-double-remission");
    
    width = 300;       // Transverse width of the slab.
    thick = 100;       // Thickness of the slab.
    size_input = 30;   // Diameter of the input beam.
    size_output = 30;  // Diameter of the target (output beam).
    separation = 60;   // Distance between the input and otput beams.
    y_input = (size_output + size_input)/2 + separation;  // Position of the center of the input beam.
    
    SquareMesh mesh;
    mesh.addRectangle(-thick/2, thick/2, -width/2, width/2, BND_OPEN);
    
    mesh.setBoundaryDisk(-thick/2, +y_input, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-thick/2, -y_input, size_input/2, DIR_WEST, BND_INPUT);
    mesh.setBoundaryDisk(-thick/2, 0, size_output/2, DIR_WEST, BND_OUTPUT);
    
    mesh.finalize();
    
    dscat = 5.;   // Scattering thickness, L/lscat, where L is the slab thickness and lscat the scattering mean free path.
    dabso = 0.;   // Absorption thickness, L/labso, where L is the slab thickness and labso the ballistic absorption length.
    tval = 0.15;  // Transmission probability.
    
    return UsadelSystem(name, mesh, dscat/thick, dabso/thick, tval);
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
    double toldf, tolr;
    
    maxit = 100;   // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). Typically: 30.
    toldf = 1e-7; // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 1;  // Verbosity level of the Newton-Raphson solver. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing fields from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Tval=" << usys.getTransmission() << ", maxit=" << maxit << ".\n";
    
    //usys.setTransmission(0.7);
    usys.initConstant();                                 // Initialize the UsadelSystem using the currently best ansatz (constants).
    usys.solveNewton(maxit, nsub, toldf, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
    //usys.setTransmission(0.75);
    //usys.solveNewton(maxit, nsub, toldf, tolr, verbose);
    
    // Save the data and plot:
    std::stringstream path;
    path << "out/" << usys.getName() << "/tval_" << usys.getTransmission() << "/field_";
    usys.savePlot(path.str());
}

/**
 * Save the distribution rho(T) to a file and then call an external script to plot it.
 */
void plotDistribution(const UsadelSystem& usys, const double* rhodata, const int ntval, 
                      const double tmin, const double tmax, const int maxit, const double elapsed_time) {
    
    const char* sep = ", ";  // Separator used in the output file.
    const int prec = 16;     // Precision used in the output file.
    
    std::stringstream path;
    path << "out/" << usys.getName() << "/tspectrum/tspectrum_";
    std::string filename_distrib;
    uniqueFilename(path.str(), ".csv", filename_distrib);  // Create a unique filename. The result is of the form "<path><number><suffix>".
    
    std::cout << TAG_INFO << "Saving distribution to file: '" << filename_distrib << "'...\n";
    std::ofstream ofs;  // Declare output stream object.
    ofs.open(filename_distrib.c_str()); // Open the file in write mode.
    ofs << std::setprecision(prec); // Set the printing precision.
    
    writeTimestamp(ofs, "%% "); // Apply a timestamp at the beginning.
    
    ofs << "%% Parameters: name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() 
        << ", h/labso=" << usys.getHolabso() << ", Trange=" << tmin << ":" << tmax << ", Ntval=" << ntval << ", maxit=" << maxit 
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
 * 
 * Arguments:
 * 
 * usys    = Usadel system on which the simulations will be performed.
 * tmin    = Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
 * tmax    = Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
 * ntval   = Number of samples for the transmission eigenvalue. Typically: 32 for quick plots, 256 for final renders (see the serial version).
 */
void computeDistributionSerial(UsadelSystem& usys, const double tmin, const double tmax, int ntval) {
    
    int maxit, nsub, verbose, niter;
    double tval, rho, toldf, tolr;
    
    maxit = 30;   // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    toldf = 1e-7; // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level of the Newton-Raphson solver. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing rho(T) from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Ntval=" << ntval << ", Trange=" << tmin << ":" << tmax << ", maxit=" << maxit << ".\n";
    
    const std::string msg = "distrib, serial";
    int cjob = 0;  // Initialize the number of completed jobs.
    auto start = std::chrono::steady_clock::now(); // Gets the current time.
    
    double* rhodata = new double[3*ntval];  // Temporary array for storing the numerical results [T, rho(T), niter].
    
    usys.setTransmission(tmin + (tmax-tmin) * (1. - std::cos(0.5*PI/ntval))/2.);   // Assigns the 'tmin' before initializing.
    usys.initConstant(); // Initialize the UsadelSystem using the currently best ansatz (constants).
    
    UsadelSystem usys_prev(usys); // Create a deep copy of the UsadelSystem containing the previous result.
    
    for (int i = 0; i < ntval; i++) {// Loop over the points of the mesh.
        
        tval = tmin + (tmax-tmin) * (1. - std::cos((i + 0.5)*PI/ntval))/2.;  // Choose Chebyshev nodes as the transmission eigenvalues (they are denser at the edges).
        usys.setTransmission(tval);   // Assigns the transmission eigenvalue.
        niter = usys.solveNewton(maxit, nsub, toldf, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
        rho = usys.getRho();  // Compute the transmission eigenvalue distribution at the given 'tval'.
        
        rhodata[3*i] = tval;  // Save the data to the array 'rhodata'.
        rhodata[3*i+1] = rho;
        rhodata[3*i+2] = niter;
        
        std::cout << TAG_INFO << "# " << (i+1) << std::string(6-std::floor(std::log10(i+1)), ' ') << std::fixed << std::setprecision(6) 
                  << "| tval=" << tval << "   rho=" << rho << "   niter=" << niter 
                  << "   " << (niter <= maxit ? "\033[92msuccess\033[0m" : "\033[91mfailure\033[0m") << std::string(50, ' ') << "\n";
        
        cjob++;
        printProgressBar(cjob, ntval, msg, start); // Print the progress bar. Must be critical since 'cjob' is a OMP-Shared variable.
        
        if (niter > maxit) {// Exit on failure to avoid loosing time (the subsequent results will fail anyway).
            std::cout << "\n" << TAG_WARN << "Detected failure, exiting now...\n";
            ntval = cjob;
            break;
        }
        else {// If successful convergence, then overwrite the previous
            usys_prev.copy(usys);
        }
        
    }
    
    double elapsed_time = endProgressBar(start);  // Finalize the progress bar and gets the total time (in seconds).
    
    // Save the distribution to a file and plot it:
    plotDistribution(usys, rhodata, ntval, tmin, tmax, maxit, elapsed_time);
    delete[] rhodata;
    
    // Save the last field (with Tmax) and plot it:
    std::stringstream path;
    path << "out/" << usys_prev.getName() << "/tval_" << usys_prev.getTransmission() << "/field_";
    usys_prev.savePlot(path.str());
}

/**
 * Compute the transmission eigenvalue distribution by solving the Usadel equation for different values of the transmission eigenvalue "tval".
 * The results are saved into a CSV file and plotted automatically by calling external scripts.
 * This is a parallelized version using OpenMP.
 * 
 * Arguments:
 * 
 * usys    = Usadel system on which the simulations will be performed.
 * tmin    = Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
 * tmax    = Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
 * ntval   = Number of samples for the transmission eigenvalue. Typically: 32 for quick plots, 256 for final renders (see the serial version).
 * nthread = Number of execution threads for OpenMP (typically the number of CPU cores).
 * 
 */
void computeDistributionOMP(UsadelSystem& usys, const double tmin, const double tmax, const int ntval, const int nthread) {
    
    int i, maxit, nsub, verbose, niter;
    double tval, rho, toldf, tolr;
    
    maxit = 30;   // Maximum number of iterations. Typically: 50-500.
    nsub = 30;    // Maximum number of substep used for backtracking line search (should between 20 and 50 in double precision). 
    toldf = 1e-7; // Tolerance over the relative displacement imposed by the Newton-Raphson step. Typically: 1e-7.
    tolr = 1e-10; // Tolerance over the norm of the residual compared to the norm of the initial residual. Typically: 1e-10.
    verbose = 0;  // Verbosity level in standard output. 0=No output, 1=Display each iteration.
    
    std::cout << TAG_INFO << "Computing rho(T) from UsadelSystem with name=" << usys.getName() << ", Npoint=" << usys.getNPoint() 
              << ", h/lscat=" << usys.getHolscat() << ", h/labso=" << usys.getHolabso() << ", Ntval=" << ntval << ", Trange=" << tmin << ":" << tmax 
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
            
            niter = usys_loc->solveNewton(maxit, nsub, toldf, tolr, verbose);  // Solve the Usadel equation using the Newton-Raphson method.
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
    //UsadelSystem usys = createDoubleWaveguide0();
    //UsadelSystem usys = createDoubleWaveguide1();
    //UsadelSystem usys = createDoubleWaveguide2();
    //UsadelSystem usys = createDoubleWaveguide3();
    //UsadelSystem usys = createDoubleWaveguide4();
    //UsadelSystem usys = createDoubleWaveguide5();
    //UsadelSystem usys = createDoubleWaveguide6();
    //UsadelSystem usys = createDoubleWaveguide7();
    //UsadelSystem usys = createDoubleWaveguide8();
    //UsadelSystem usys = createDoubleWaveguide9();
    //UsadelSystem usys = createDoubleWaveguide10();
    //UsadelSystem usys = createDoubleWaveguide11();
    //UsadelSystem usys = createDoubleWaveguide12();
    //UsadelSystem usys = createWaveguideObstacle1();
    //UsadelSystem usys = createCircularDoubleWaveguide1();
    //UsadelSystem usys = createCircularOpposite();
    //UsadelSystem usys = createCircularCavity();
    //UsadelSystem usys = createCircularCavityBig();
    //UsadelSystem usys = createCircularCavityAlt();
    //UsadelSystem usys = createCircularCavityHole();
    //UsadelSystem usys = createEiffelTower();
    //UsadelSystem usys = createFiniteSlabTransmission1();
    //UsadelSystem usys = createFiniteSlabTransmission2();
    //UsadelSystem usys = createFiniteSlabTransmission3();
    //UsadelSystem usys = createSlabTransmission();
    //UsadelSystem usys = createFiniteSlabRemission1();
    //UsadelSystem usys = createFiniteSlabRemission2();
    //UsadelSystem usys = createFiniteSlabRemission3();
    //UsadelSystem usys = createFiniteSlabDoubleRemission();
    
    double tmin, tmax, dscat, dabso, holscat, holabso;
    int ntval, nthread;
    
    // Constructs the mesh from a PNG file:
    SquareMesh mesh("model/waveguide_102x100.png");
    //SquareMesh mesh("model/slab-transmission-1_101x299.png");
    
    dscat = 10.;  // Scattering depth, L/lscat. Default: dscat=8.6 (in order to get approximately dscat_eff=10).
    dabso = 0.;   // Absorption depth, L/labso.
    
    const std::string sysname = "waveguide_102x100/dscat_" + to_string_prec(dscat, 6);
    
    holscat = dscat/100.;
    holabso = dabso/100.;
    
    UsadelSystem usys(sysname, mesh, holscat, holabso, 0.64);
    
    //plotMesh(usys);  // Plot the mesh only to check it is as expected.
    computeFields(usys); // Compute the fields (theta, eta, and Q) and the intensity profile for the given transmission eigenvalue.
    
    //tmin = 0.10;   // Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    //tmax = 0.001;   // Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    //ntval = 32;    // Number of samples for the transmission eigenvalue. Typically: 64.
    //computeDistributionSerial(usys, tmin, tmax, ntval); // Compute the transmission eigenvalue distribution rho(T) by scanning in T.
    
    //tmin = 0.;    // Minimum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    //tmax = 1.;    // Maximum transmission eigenvalue. Note that this value is never exactly reached due to the Chebyshev nodes.
    //ntval = 40;   // Number of samples for the transmission eigenvalue. Typically: 4*nthread for quick plots.
    //nthread = 10; // Number of execution threads for OpenMP (typically the number of CPU cores).
    //computeDistributionOMP(usys, tmin, tmax, ntval, nthread); // Compute the transmission eigenvalue distribution rho(T). Parallelized version.
    
    return 0;
}
