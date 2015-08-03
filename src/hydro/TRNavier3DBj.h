/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RN_RNAVIER3DBJ_UTIL_H
#define RN_RNAVIER3DBJ_UTIL_H

#include<string>
#include<fstream>
#include<memory>
#include "THModel3D.h"

class TEOS ;

class TRNavier3DBj {

    private:

       std::string fEOSName;
       std::string fModelName;
       std::string fNametag ;
       std::ifstream fInputFile ;

       std::unique_ptr<TEOS> fEOS ;
       std::unique_ptr<THModel3D> fModel ;
       void read_rnavier3dbj(std::istream &in)  ;

    public:

       ~TRNavier3DBj() ;

       TRNavier3DBj(std::string nametag, std::unique_ptr<TEOS> (*rnavier3dbj_make_eos)(TRNavier3DBj *rn, const std::string &eosname))  ;

       std::string &getEOSName() { return fEOSName ; }
       std::string &getModelName() { return fModelName ; }
       std::string &getNametag() { return fNametag; }

       std::ifstream &getInputFile() {return fInputFile ;}

       TEOS *getEOS() {return fEOS.get() ; }
       THModel3D *getHModel3D() {return fModel.get() ; }

} ;


#endif
