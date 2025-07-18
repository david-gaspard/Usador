/****
 * @date Created on 2025-06-27 at 16:23:02 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SquareMesh class.
 ***/
#ifndef _SQUARE_MESH_H
#define _SQUARE_MESH_H
#include "Constants.hpp"
#include "MeshPoint.hpp"
#include "Vector2D.hpp"

/**
 * Enumeration specifying the four possible directions on a square lattice:
 */
enum Direction {NORTH, SOUTH, EAST, WEST};

/**
 * Class defining the Square Mesh object.
 */
class SquareMesh {
    
    private:
    
    std::vector<MeshPoint> point; // List of positions of points. The points must be integers (belong to the Cartesian lattice). They also contain neighbor indices.
    bool ready;  // True when the SquareMesh is ready for computations. False, otherwise. This boolean becomes True after fixNeighbors() is called.
    
    public:
    
    // Constructors/Destructors:
    SquareMesh();
    ~SquareMesh();
    
    // Getters:
    MeshPoint getPoint(const int i) const;
    int getNPoint() const;
    int getNBoundary(const int bndtype) const;
    
    // Determine the index of a given point:
    int indexOf(const int x, const int y) const;
    bool containsPoint(const int x, const int y) const;
    
    // Add points:
    void addPoint(const int x, const int y);
    void addRectangle(int xmin, int xmax, int ymin, int ymax);
    void addDisk(const int x0, const int y0, const double radius);
    void addPolygon(const std::vector<Vector2D>& polygon);
    void addPolygon(const char* filename, const double scale);
    
    // Remove points:
    void removePoint(const int x, const int y);
    void removeRectangle(int xmin, int xmax, int ymin, int ymax);
    void removeDisk(const int x0, const int y0, const double radius);
    void removePolygon(const std::vector<Vector2D>& polygon);
    
    // Assign boundary conditions:
    void setBoundaryPoint(const int x, const int y, const Direction dir, const int bndtype);
    void setBoundaryRegion(int xmin, int xmax, int ymin, int ymax, const Direction dir, const int bndtype);
    
    // Finalize the mesh:
    void fixNeighbors();
    
    // Output methods:
    void printSummary() const;
    void saveMesh(const char* filename, const char* sep) const;
    void saveMeshShort(const char* filename, const char* sep) const;
    
};

#endif
