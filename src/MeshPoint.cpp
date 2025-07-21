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
MeshPoint::MeshPoint(const int x, const int y, const int bndtype) {
    
    this->x = x;
    this->y = y;
    
    north = bndtype;
    south = bndtype;
    east  = bndtype;
    west  = bndtype;
}

/**
 * Default constructor of MeshPoint.
 */
MeshPoint::MeshPoint() : MeshPoint(0, 0, BND_MIRROR) {}

/**
 * Returns a string version of the given boundary type.
 */
const std::string boundaryTypeString(const int bndtype) {
    switch (bndtype) {
        case BND_MIRROR:
            return "mirror";
        case BND_OPEN:
            return "open";
        case BND_INPUT:
            return "input";
        case BND_OUTPUT:
            return "output";
        default:
            return std::to_string(bndtype);
    }
}
