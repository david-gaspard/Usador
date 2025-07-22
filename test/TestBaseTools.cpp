/****
 * @date Created on 2025-07-22 at 13:35:40 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ code providing test for the BaseTools utilities.
 ***/
#include "BaseTools.hpp"
#include "Constants.hpp"
#include <thread>
#include <iostream>

/**
 * Test the progress bar.
 */
int testProgressBar() {
    
    using namespace std::chrono_literals;
    
    const int njob = 60;  // Number of jobs.
    const std::string msg("TestBaseTools");
    
    auto start = std::chrono::steady_clock::now();
    
    for (int i = 1; i <= njob; i++) {
        
        std::this_thread::sleep_for(0.1s);
        
        printProgressBar(i, njob, msg, start);
        
    }
    
    double elapsed_time = endProgressBar(start);
    std::cout << TAG_INFO << "Returned time = " << elapsed_time << " s.\n";
    
    return 0;
}

/**
 * Main function of this testing unit.
 */
int main(int argc, char** argv) {
    
    testProgressBar();
    
    return 0;
}
