/****
 * @date Created on 2025-06-27 at 16:23:02 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the SquareMesh class.
 ***/
#ifndef _SQUARE_MESH_H
#define _SQUARE_MESH_H
#include "MeshPoint.hpp"
#include "Constants.hpp"
#include <vector>

/**
 * Class defining the Square Mesh object.
 */
class SquareMesh {
    
    private:
    
    std::vector<MeshPoint> point; // List of positions of points. The points must be integers (belong to the Cartesian lattice). They also contain neighbor indices.
    bool ready;  // True when the SquareMesh is ready for computations. False, otherwise. This boolean becomes True after fixNeighbors() is called.
    
    public:
    
    SquareMesh();
    ~SquareMesh();
    MeshPoint getPoint(const int i) const;
    int getNPoint() const;
    int indexOf(const int x, const int y) const;
    bool containsPoint(const int x, const int y) const;
    void addPoint(const int x, const int y);
    void addRectangle(int xmin, int xmax, int ymin, int ymax);
    void addDisk(const int x0, const int y0, const double radius);
    void removePoint(const int x, const int y);
    void removeRectangle(int xmin, int xmax, int ymin, int ymax);
    void setBoundaryPoint(const int x, const int y, const int bndtype);
    void setBoundaryRegion(int xmin, int xmax, int ymin, int ymax, int bndtype);
    void fixNeighbors();
    void printSummary() const;
    void saveMesh(const char* filename, const char* sep, const int verbosity) const;
    
};

#endif
