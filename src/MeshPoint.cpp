/****
 * @date Created on 2025-07-03 at 13:57:17 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing the MeshPoint methods.
 ***/
#include "MeshPoint.hpp"

/**
 * Constructor of a MeshPoint with default unset boundaries.
 */
MeshPoint::MeshPoint(const int x, const int y) {
    
    this->x = x;
    this->y = y;
    
    north = BND_DEFAULT;
    south = BND_DEFAULT;
    east  = BND_DEFAULT;
    west  = BND_DEFAULT;
}

/**
 * Default constructor of MeshPoint.
 */
MeshPoint::MeshPoint() : MeshPoint(0, 0) {}
