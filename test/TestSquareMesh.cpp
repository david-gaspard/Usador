/****
 * @date Created on 2025-06-27 at 16:32:13 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code to test the SquareMesh class.
 ***/
#include "SquareMesh.hpp"
#include <iostream>

int main(int argc, char** argv) {
    
    SquareMesh mesh; // Declare the mesh.
    
    // Setup the mesh from geometry:
    mesh.addRectangle(-10, 10, -5, 5);
    mesh.addDisk(10, 5, 7.5);
    mesh.removeRectangle(9, 11, 4, 6);
    
    // Setup boundary conditions (only vacuum because mirror is default):
    mesh.setBoundaryRegion(-10, -10, -5, 5, BND_OPEN);
    mesh.setBoundaryRegion( 10,  10, -5, 5, BND_OPEN);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.fixNeighbors();
    
    // Save the mesh to a file for checking:
    mesh.saveMesh("test_mesh_short.csv", ", ", 0);
    mesh.saveMesh("test_mesh_full.csv", ", ", 1);
    
    // Print the position of a given point:
    const int i0 = 100;
    MeshPoint p = mesh.getPoint(i0);
    std::cout << TAG_INFO << "Point(" << i0 << ") = (" << p.x << ", " << p.y << ").\n";
    
    return 0;
}
