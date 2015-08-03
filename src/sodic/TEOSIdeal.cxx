/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include <iostream>
#include "numeric.h"
#include "TRNavier3DBj.h"
#include "TEOSIdeal.h"
#include "TIOBuffer.h"

const double TEOSIdeal::fDim=3. ;
void TEOSIdeal::read(istream &in)
{

   TIOBuffer b ;
   b.read_section(in, "TEOSIdeal") ;
   b.getD("fGammaIndex", fGammaIndex) ;
   b.getD("fEtaOverS", fEtaOverS) ;
   b.getD("fSigmaOverS", fSigmaOverS) ;
   b.getD("fKappaTOverS", fKappaTOverS) ;
   b.getD("fT0", fT0) ;
   b.getD("fE0", fE0) ;

   b.getD("fTPi_EtaST", fTPi_EtaST) ;
   b.getD("fL1_EtaTPi", fL1_EtaTPi) ;
   b.getD("fL2_EtaTPi", fL2_EtaTPi) ;

   //Calculate  computed quantitites.
   fS0  = fE0 * fGammaIndex/fT0 ;
   //This code only works if fGammaIndex == 1.33333333 ;
   fDOF = fE0/(M_PI*M_PI/30.)/gsl_pow_4(fT0) ;

}

//! Write the Gamma Law equation of state to out
void TEOSIdeal::write(ostream &out) 
{
   TIOBuffer b ;
   out << "[TEOSIdeal]" << endl;
   b.writeLineD(out,fGammaIndex, "fGammaIndex : Adiabatic  Index") ;
   b.writeLineD(out,fEtaOverS,   "fEtaOverS   : shear viscosity /s ") ;
   b.writeLineD(out,fSigmaOverS, "fSigmaOverS : bulk  viscosity /s ") ;
   b.writeLineD(out,fKappaTOverS,"fKappaTOverS: conductivity / s  ") ;
   b.writeLineD(out,fT0, "fT0: T0 ") ;
   b.writeLineD(out,fE0, "fE0: E0 ") ;
   b.writeLineD(out,fTPi_EtaST, "fTPi_EtaST: TauPi/(eta/s T) ") ;
   b.writeLineD(out,fL1_EtaTPi, "fL1_EtaTPi: lambda_1 /(eta * tau_pi)") ;
   b.writeLineD(out,fL2_EtaTPi, "fL2_EtaTPi: lambda_2 /(eta * tau_pi)") ;

   out << "[END]" << endl;
}
double TEOSIdeal::eofs(const double &s, const double &n) 
{
   double ginv = 1./fGammaIndex ;
   double e =  fE0 * pow( s/fS0, fGammaIndex ) + (1.-ginv)*n;
   return e ;
}
   
void TEOSIdeal::eos(const double &e, const double &n, double &p, double &cs) 
{
   p  = (fGammaIndex-1.)*(e-n) ;
   cs = sqrt(fGammaIndex - 1.0) ;
}


void TEOSIdeal::stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) 
{
   double ginv = 1./fGammaIndex ;
   s  =  fS0 * pow(((e-(1-ginv)*n)/fE0),ginv) ;
   t  =  (fGammaIndex * e-n*(fGammaIndex-1)) / s ;
}

void TEOSIdeal::viscosity(const double &e, const double &n, 
         double &sigma_overs, double &kappaT_overs, double &eta_overs) 
{
   sigma_overs  = fSigmaOverS ;
   kappaT_overs = fKappaTOverS ;
   eta_overs    = fEtaOverS ;
}



std::unique_ptr<TEOS> make_eos_eosideal(TRNavier3DBj *rn, const std::string &icname) 
{
   TEOSIdeal *eos = new TEOSIdeal ;
   eos->read(rn->getInputFile()) ;
   eos->write(clog) ;
   return std::unique_ptr<TEOS>(eos) ;
}

