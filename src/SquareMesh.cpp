/****
 * @date Created on 2025-06-27 at 10:06:14 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing a two-dimensional square mesh object.
 ***/
#include "SquareMesh.hpp"
#include "BaseTools.hpp"
#include <iostream>
#include <cmath>

/**
 * Constructor of the Square Mesh object.
 */
SquareMesh::SquareMesh() {
    //std::cout << TAG_INFO << "Creating SquareMesh.\n";
    ready = false;
}

/**
 * Destructor.
 */
SquareMesh::~SquareMesh() {
    //std::cout << TAG_INFO << "Deleting SquareMesh.\n";
}

/**
 * Return the point at position "i" in the mesh.
 */
MeshPoint SquareMesh::getPoint(const int i) const {
    if (not ready) {
        throw std::logic_error("In getPoint(): SquareMesh is not completely initialized. Please use fixNeighbors().");
    }
    return point.at(i);  // vector objects already make bound checking.
}

/**
 * Returns the number of points in the mesh.
 */
int SquareMesh::getNPoint() const {
    if (not ready) {
        throw std::logic_error("In getNPoint(): SquareMesh is not completely initialized. Please use fixNeighbors().");
    }
    return point.size();
}

/**
 * Return the index of a given point in the mesh.
 * If the point does not exist, then return BND_DEFAULT (a negative value).
 */
int SquareMesh::indexOf(const int x, const int y) const {
    MeshPoint p;
    for (unsigned int i = 0; i < point.size(); i++) {
        p = point.at(i);
        if (p.x == x && p.y == y){
            return i;
        }
    }
    return BND_DEFAULT;
}

/**
 * Returns true if the mesh contains the given point (x, y).
 */
bool SquareMesh::containsPoint(const int x, const int y) const {
    return indexOf(x, y) >= 0;
}

/**
 * Add a single point to the mesh only if it is not already contained in the mesh.
 */
void SquareMesh::addPoint(const int x, const int y) {
    if (not containsPoint(x, y)) {
        point.push_back(MeshPoint(x, y)); // Add point with default neighboring value.
    }
}

/**
 * Add a square region to the mesh.
 */
void SquareMesh::addRectangle(int xmin, int xmax, int ymin, int ymax) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            addPoint(x, y);
        }
    }
}

/**
 * Add a circular region to the mesh.
 */
void SquareMesh::addDisk(int x0, int y0, double radius) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int)  std::ceil(x0 + radius);
    int ymin = (int) std::floor(y0 - radius);
    int ymax = (int)  std::ceil(y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                addPoint(x, y);
            }
        }
    }
}

/**
 * Add a polygon region to the mesh. Uses the even-odd filling rule.
 */
void SquareMesh::addPolygon(const std::vector<Vector2D>& polygon) {
    
    // 1. Determine the bounds of the polygon:
    int xmin, xmax, ymin, ymax;
    polygonBounds(polygon, xmin, xmax, ymin, ymax);
    
    // 2. Loop on points in the rectangular region:
    Vector2D p;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            p = Vector2D(x, y);   // Current point that we attempt to add to the mesh.
            if (p.windingNumber(polygon) % 2 != 0) {// Uses the even-odd rule to fill the polygon (fill if the winding number is odd).
                addPoint(x, y);
            }
        }
    }
}

/**
 * Add a polygon using the coordinates given by a file.
 */
void SquareMesh::addPolygon(const char* filename, const double scale) {
    
    std::vector<Vector2D> polygon;
    std::ifstream ifs(filename);  // Open the file.
    std::string line;
    double x, y;
    size_t sz;
    
    if (ifs.fail()) {// Check if the file is opened.
        std::string info = "In addPolygonFromFile(): Cannot open file '" + std::string(filename) + "'.";
        throw std::runtime_error(info);
    }
    
    while (getline(ifs, line)) {// Loop on each line of the file.
        try {
            x = scale * std::stod(line, &sz);
            y = scale * std::stod(line.substr(sz));
            polygon.push_back(Vector2D(x, y));
        }
        catch (const std::invalid_argument& exc) {
            std::cerr << TAG_WARN << "In addPolygonFromFile(): Argument is not a number, skipping line (what = '" 
                      << exc.what() << "').\n";
        }
    }
    
    ifs.close();
    
    addPolygon(polygon);  // Add the parsed polygon.
}

/**
 * Remove the point at a given position.
 */
void SquareMesh::removePoint(const int x, const int y) {
    int i = indexOf(x, y);
    if (i >= 0) point.erase(point.begin() + i);
}

/**
 * Remove a rectangular region of points.
 */
void SquareMesh::removeRectangle(int xmin, int xmax, int ymin, int ymax) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            removePoint(x, y);
        }
    }
}

/**
 * Remove a circular region from the mesh.
 */
void SquareMesh::removeDisk(const int x0, const int y0, const double radius) {
    int xmin = (int) std::floor(x0 - radius);
    int xmax = (int)  std::ceil(x0 + radius);
    int ymin = (int) std::floor(y0 - radius);
    int ymax = (int)  std::ceil(y0 + radius);
    int r2 = (int) std::ceil(radius*radius);
    int dx, dy;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            dx = x - x0;
            dy = y - y0;
            if (dx*dx + dy*dy <= r2) {
                removePoint(x, y);
            }
        }
    }
}

/**
 * Removes a polygon region from the mesh. Uses the even-odd rule.
 */
void SquareMesh::removePolygon(const std::vector<Vector2D>& polygon) {
    
    // 1. Determine the bounds of the polygon:
    int xmin, xmax, ymin, ymax;
    polygonBounds(polygon, xmin, xmax, ymin, ymax);
    
    // 2. Loop on points in the rectangular region:
    Vector2D p;
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            p = Vector2D(x, y);   // Current point that we attempt to add to the mesh.
            if (p.windingNumber(polygon) % 2 != 0) {// Uses the even-odd rule to remove the polygon.
                removePoint(x, y);
            }
        }
    }
}

/**
 * Setup the boundary condition for the point at the given position (x, y).
 * This function must be called after fix_neighbors() methods.
 */
void SquareMesh::setBoundaryPoint(const int x, const int y, const Direction dir, const int bndtype) {
    int i = indexOf(x, y);
    if (i >= 0) {//Only works if the point is in the mesh.
        MeshPoint& p = point.at(i);
        if (dir == NORTH && p.north < 0) p.north = bndtype;
        else if (dir == SOUTH && p.south < 0) p.south = bndtype;
        else if (dir == EAST && p.east < 0) p.east  = bndtype;
        else if (dir == WEST && p.west < 0) p.west  = bndtype;
    }
}

/**
 * Setup the same boundary condition for a region of points.
 * This function must be called after fix_neighbors() methods.
 */
void SquareMesh::setBoundaryRegion(int xmin, int xmax, int ymin, int ymax, const Direction dir, const int bndtype) {
    
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);
    
    for (int y = ymax; y >= ymin; y--) {// Loop on the square lattice in reading order.
        for (int x = xmin; x <= xmax; x++) {
            setBoundaryPoint(x, y, dir, bndtype);
        }
    }
}

/**
 * Compute the list of nearest neighbors.
 * This method must be called after the mesh construction methods add_*() and remove_*() and before the set_boundary() methods.
 */
void SquareMesh::fixNeighbors() {
    int inorth, isouth, ieast, iwest;
    for (MeshPoint& p : point) {
        inorth = indexOf(p.x, p.y+1);
        isouth = indexOf(p.x, p.y-1);
        ieast  = indexOf(p.x+1, p.y);
        iwest  = indexOf(p.x-1, p.y);
        if (inorth >= 0) p.north = inorth;
        if (isouth >= 0) p.south = isouth;
        if (ieast  >= 0) p.east  = ieast;
        if (iwest  >= 0) p.west  = iwest;
    }
    ready = true; // The SquareMesh gets ready for computations.
}

/**
 * Print a summary of the mesh.
 */
void SquareMesh::printSummary() const {
    std::cout << TAG_INFO << "SquareMesh with " << point.size() << " points.\n";
}

/**
 * Save the points to a file "filename" using the separator "sep".
 */
void SquareMesh::saveMesh(const char* filename, const char* sep) const {
    
    if (not ready) {
        throw std::logic_error("In saveMesh(): SquareMesh is not completely initialized. Please use fixNeighbors().");
    }
    
    std::ofstream ofs;
    ofs.open(filename);
    
    writeTimestamp(ofs, "%% ");
    
    ofs << "%% SquareMesh with " << point.size() << " points.\n"
        << "x" << sep << "y" << sep << "north" << sep << "south" << sep << "east" << sep << "west\n";
    
    for (MeshPoint p : point) {
        ofs << p.x << sep << p.y << sep << boundaryTypeString(p.north) << sep 
                                        << boundaryTypeString(p.south) << sep 
                                        << boundaryTypeString(p.east) << sep 
                                        << boundaryTypeString(p.west) << "\n";
    }
    
    ofs.close();
}

/**
 * Save the points to a file "filename" using the separator "sep".
 * This is a short version for quick plots with less post-processing.
 */
void SquareMesh::saveMeshShort(const char* filename, const char* sep) const {
    
    if (not ready) {
        throw std::logic_error("In saveMeshShort(): SquareMesh is not completely initialized. Please use fixNeighbors().");
    }
    
    std::ofstream ofs;
    ofs.open(filename);
    
    writeTimestamp(ofs, "%% ");
    
    ofs << "%% SquareMesh with " << point.size() << " points.\n" 
        << "x" << sep << "y" << sep << "bnd\n";
    
    for (MeshPoint p : point) {
        ofs << p.x << sep << p.y << sep;
        if (p.north == BND_INPUT || p.south == BND_INPUT || p.east == BND_INPUT || p.west == BND_INPUT) {
            ofs << "input";
        }
        else if (p.north == BND_OUTPUT || p.south == BND_OUTPUT || p.east == BND_OUTPUT || p.west == BND_OUTPUT) {
            ofs << "output";
        }
        else if (p.north == BND_OPEN || p.south == BND_OPEN || p.east == BND_OPEN || p.west == BND_OPEN) {
            ofs << "open";
        }
        else if (p.north == BND_MIRROR || p.south == BND_MIRROR || p.east == BND_MIRROR || p.west == BND_MIRROR) {
            ofs << "mirror";
        }
        else {
            ofs << "bulk";
        }
        ofs << "\n";
    }
    
    ofs.close();
}
