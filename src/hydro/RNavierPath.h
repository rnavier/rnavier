/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RNAVIER_PATH_H
#define RNAVIER_PATH_H

#include <string>

// Retuns a string containing the path for data files. This string is 
// set by the environment variable RNAVIERDATA on first call.
const std::string &GetRNAVIERDATA() ;

// Retuns a string containing the nametag for the current run. 
void SetRNavierNametag(const std::string &nametag) ;
const std::string &GetRNavierNametag() ;

#endif
