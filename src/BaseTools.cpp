/****
 * @date Created on 2025-07-09 at 16:09:04 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing basic tools for file manipulation.
 ***/
#include "BaseTools.hpp"
#include "Constants.hpp"
#include <ctime>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <iostream>

/**
 * Write a one-liner timestamp in the given stream with the given prefix (typically a comment character when writing to a file).
 */
void writeTimestamp(std::ofstream& ofs, const char* prefix) {
    std::time_t time_now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ofs << prefix << "Computed on " << std::put_time(std::localtime(&time_now), "%F at %T %z") << " by " << PROGRAM_COPYRIGHT << "\n";
}

/**
 * Create a unique filename at the given "path" and with the given "suffix". The unique filename found is written in the string "unique_filename".
 * The format is: <path><i><suffix> where <i> is an integer at least equal to 1 which is incremented until the filename is unique.
 * This function ensures data is not overwritten. If no unique filename is found, then prints a warning message.
 * This functions also creates the appropriate subdirectories if necessary, and warns the user if so.
 */
void uniqueFilename(const std::string& path, const std::string& suffix, std::string& unique_filename) {
    
    // 1. Ensure the output directory exists:
    std::string outputdir(path.substr(0, path.find_last_of("/")));
    //std::cout << TAG_INFO << "Output directory: " << outputdir << "\n";
    if (std::filesystem::create_directories(outputdir)) {// Create the directory if necessary.
        std::cout << TAG_INFO << "Created directory: '" << outputdir << "'.\n";
    }
    
    // 2. Create a unique filename:
    for (int i = 1; i <= MAX_FILENO; i++) {
        unique_filename = path + std::to_string(i) + suffix;
        if (!std::filesystem::exists(unique_filename)) return; // If the file does not exist, then exits the function.
    }
    
    std::cout << TAG_WARN << "Could not find unique filename: '" << unique_filename << "'. Data will be lost...\n";
}
