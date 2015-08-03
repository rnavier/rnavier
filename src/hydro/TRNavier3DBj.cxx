#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"
#include "THModel3D.h"
#include "RNavierPath.h"


using namespace std ;

//! Read in the model
unique_ptr<THModel3D> read_model(TRNavier3DBj *rn, const string &eosname) ;

TRNavier3DBj::~TRNavier3DBj() { ; }

TRNavier3DBj::TRNavier3DBj(string nametag_ini, unique_ptr<TEOS> (*rnavier3dbj_make_eos)(TRNavier3DBj *rn, const string &eosname))  
{
    // strip off the .ini  to form the nametag
   int lastindex = nametag_ini.find_last_of("."); 
   fNametag = nametag_ini.substr(0, lastindex); 
    
   // Define the nametag for later use
   SetRNavierNametag(fNametag) ;

   // Open the input file
   fInputFile.open(nametag_ini) ;
   if (fInputFile.fail()) {
      cerr << "** TRNavier3DBj ** bad input filename " << nametag_ini << endl;
      exit(EXIT_FAILURE) ;
   }
   // Print out the name of the data inputfile
   cout << "# The data path is " << GetRNAVIERDATA() << endl;

   // Read the settings from the input file
   read_rnavier3dbj(fInputFile) ; 

   //fRNRandom = new TRNRandom(this) ;
   fEOS = rnavier3dbj_make_eos(this, getEOSName()) ;
   fModel = read_model(this, getModelName()) ;
}

void TRNavier3DBj::read_rnavier3dbj(istream &in)
{
   TIOBuffer buffer ;
   buffer.read_section(in, "TRNavier3DBj") ;
   buffer.getS("fEOSName",  fEOSName) ;
   buffer.getS("fModelName", fModelName) ;

   clog << "[TRNavier3DBj]" <<endl ;
   buffer.writeLineS(clog, fEOSName, "fEOSName : The EOS") ;
   buffer.writeLineS(clog, fModelName, "fModelName : The model") ;
   clog << "\n" <<endl ;
}

unique_ptr<THModel3D> read_model(TRNavier3DBj *rn, const string &name) 
{
   if (name == "TBRSSS") {
      return unique_ptr<THModel3D>(new TBRSSS(rn)) ;
   } else {
      clog << "Unrecognized type! \n" << name << endl ;
      abort() ;
      return unique_ptr<THModel3D>(nullptr) ;
   }
}

