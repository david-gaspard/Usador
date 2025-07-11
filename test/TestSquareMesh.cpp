/****
 * @date Created on 2025-06-27 at 16:32:13 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code to test the SquareMesh class.
 ***/
#include "SquareMesh.hpp"
#include <iostream>

/**
 * Test the basic functions of SquareMesh object.
 */
int testMeshBasic() {
    std::cout << "====== TEST SQUARE MESH BASIC ======\n";
    
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
    const char* filename_short = "test_mesh_short.csv";
    const char* filename_full = "test_mesh_full.csv";
    const char* sep = ", ";
    std::cout << TAG_INFO << "Saving SquareMesh to files '" << filename_short << "' and '" << filename_full << "'...\n";
    mesh.saveMesh(filename_short, sep, 0);
    mesh.saveMesh(filename_full,  sep, 1);
    
    // Print the position of a given point:
    const int i0 = 100;
    MeshPoint p = mesh.getPoint(i0);
    std::cout << TAG_INFO << "Point(" << i0 << ") = (" << p.x << ", " << p.y << ").\n";
    
    return 0;
}

/**
 * Test the polygon generation of SquareMesh object.
 */
int testMeshPolygon1() {
    std::cout << "====== TEST SQUARE MESH POLYGON #1 ======\n";
    
    std::vector<Vector2D> polygon = {Vector2D(10, 0), Vector2D(4, 8), Vector2D(-7, 2), Vector2D(-5, -11), Vector2D(1, -3)};
    
    SquareMesh mesh; // Declare the mesh.
    
    mesh.addPolygon(polygon);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.fixNeighbors();
    
    // Save the mesh to a file for checking:
    const char* filename = "test_mesh_polygon_1.csv";
    std::cout << TAG_INFO << "Saving SquareMesh to file '" << filename << "'...\n";
    mesh.saveMesh(filename, ", ", 0);
    
    return 0;
}

/**
 * Test the polygon generation of SquareMesh object.
 */
int testMeshPolygon2() {
    std::cout << "====== TEST SQUARE MESH POLYGON #2 ======\n";
    
    std::vector<Vector2D> polygon = {
        Vector2D( 25, 0  ), Vector2D( 20, 12 ), Vector2D( 10, 20 ), Vector2D( -5, 18 ), Vector2D(-16, 13 ),
        Vector2D(-27, 4  ), Vector2D(-25, -6 ), Vector2D(-19, -8 ), Vector2D(-12, -15), Vector2D( -3, -32),
        Vector2D(  5, -35), Vector2D( 14, -23), Vector2D( 10, -11), Vector2D(  2, -9 ), Vector2D( -3, -15),
        Vector2D(  1, -24), Vector2D( 12, -21), Vector2D( 22, -10)
    };
    
    SquareMesh mesh; // Declare the mesh.
    
    mesh.addPolygon(polygon);
    
    // Fix the nearest neighbors to finalize. This operation is mandatory:
    mesh.fixNeighbors();
    
    // Save the mesh to a file for checking:
    const char* filename = "test_mesh_polygon_2.csv";
    std::cout << TAG_INFO << "Saving SquareMesh to file '" << filename << "'...\n";
    mesh.saveMesh(filename, ", ", 0);
    
    return 0;
}

/**
 * Main function of the test of the SquareMesh object.
 */
int main(int argc, char** argv) {
    
    //testMesh1();
    //testMeshPolygon1();
    testMeshPolygon2();
    
}
