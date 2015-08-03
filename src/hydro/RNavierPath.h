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
