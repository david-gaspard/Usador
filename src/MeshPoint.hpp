/****
 * @date Created on 2025-07-03 at 13:52:26 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing the MeshPoint definition.
 ***/
#ifndef _MESHPOINT_H
#define _MESHPOINT_H
#include <string>

// Constants attached to the MeshPoint:
static const int BND_MIRROR  = -1;  // Index of "neigh" used for mirror boundary conditions (zero current condition).
static const int BND_OPEN    = -2;  // Index of "neigh" used for open boundary conditions (extrapolated-length boundary conditions).
static const int BND_INPUT   = -3;  // Index of "neigh" used for input boundary condition.
static const int BND_OUTPUT  = -4;  // Index of "neigh" used for input boundary condition.
//static const int BND_DEFAULT = BND_MIRROR; // This is the default boundary condition (zero current condition).

/**
 * Define the integer point type used in the square mesh.
 */
struct MeshPoint {
    
    int x;     // X coordinate of the point.
    int y;     // Y coordinate of the point.
    int north; // Index of the nearest neighboring point at the north (x, y+1). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int south; // Index of the nearest neighboring point at the south (x, y-1). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int east;  // Index of the nearest neighboring point at the east (x+1, y). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    int west;  // Index of the nearest neighboring point at the west (x-1, y). If not in the mesh, it must be BND_OPEN or BND_MIRROR.
    
    // Constructors:
    MeshPoint(const int x, const int y, const int bndtype);
    MeshPoint();
    
};

const std::string boundaryTypeString(const int bndtype);

#endif
