/****
 * @date Created on 2025-07-09 at 16:11:42 CEST
 * @author David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>
 * @copyright This program is distributed under the MIT License.
 * @file C++ header providing basic tools.
 ***/
#ifndef _BASE_TOOLS_H
#define _BASE_TOOLS_H
#include <fstream>

void writeTimestamp(std::ofstream& ofs, const char* prefix);

#endif
