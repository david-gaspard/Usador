/****
 * @date Created on 2025-07-01 at 15:20:52 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing global constants.
 ***/
#ifndef _CONSTANTS_H
#define _CONSTANTS_H
#include <complex>
#include <string>

// Define the complex type:
typedef std::complex<double> dcomplex;

// Printing tags:
static const std::string TAG_INFO = "[INFO] ";  // Information tag.
static const std::string TAG_WARN = "[WARN] ";  // Warning tag.
static const std::string TAG_ERROR = "[ERROR] ";  // Error tag.
static const std::string PROGRAM_NAME_SHORT = "Usador v0.1";  // Information tag.
static const std::string PROGRAM_NAME_FULL  = PROGRAM_NAME_SHORT + " - Usadel equation solver for arbitrary disordered regions";
static const std::string PROGRAM_COPYRIGHT  = PROGRAM_NAME_SHORT + " (c) 2025 David GASPARD (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>";

// Mathematical constants:
static const double   PI = 3.1415926535897932384626;  // The fundamental constant Pi = 3.1415926535897932384626...
static const dcomplex I  = dcomplex(0.0, 1.0);  // Define the imaginary unit.

// Physical constants:
static const int    DIMENSION  = 2;      // Number of spatial dimensions of the mesh, which is always equal to 2.
static const double EXTRAPOLEN = PI/4.;  // Diffusive extrapolation length in 2D. Exact: 0.818309886184. Std approx: pi/4 = 0.785398163397.

// Numerical constants:
static const double MEPS = 1.11e-16;     // Machine epsilon in double precision, MEPS = 2^(-53).
static const double SQRTEPS = 1.053e-8;  // Square root of the machine epsilon double precision, SQRTEPS = MEPS^(1/2) = 2^(-53/2), used for discrete derivatives with forward difference formula.
static const double CBRTEPS = 4.806e-6;  // Cubic root of the machine epsilon in double precision, CBRTEPS = MPS^(1/3) = 2^(-53/3), used for discrete derivatives with central difference formula.

#endif
