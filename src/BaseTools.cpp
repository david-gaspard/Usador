/****
 * @date Created on 2025-07-09 at 16:09:04 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing basic tools.
 ***/
#include "BaseTools.hpp"
#include "Constants.hpp"
#include <ctime>
#include <chrono>
#include <iomanip>

/**
 * Write a one-liner timestamp in the given stream with the given prefix (typically a comment character when writing to a file).
 */
void writeTimestamp(std::ofstream& ofs, const char* prefix) {
    std::time_t time_now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ofs << prefix << "Computed on " << std::put_time(std::localtime(&time_now), "%F at %T %z") << " by " << PROGRAM_COPYRIGHT << "\n";
}


