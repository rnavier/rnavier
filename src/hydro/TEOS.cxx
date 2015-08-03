#include <iostream>
#include "numeric.h"
#include "TEOS.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"

double TEOS::getBRSSSTauPi(const double & e, 
        const double &n, 
        const double &s,
        const double &temper, 
        const double &mu)
{
    double sigma_overs, kappaT_overs, eta_overs ;
    double tpi_etast, l1_ntpi, l2_ntpi ;
    viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
    getBRSSSParams(e, n, tpi_etast, l1_ntpi, l2_ntpi) ;
    return tpi_etast * eta_overs / temper ;
}

const double TGammaLaw::fDim=3. ;

TGammaLaw::TGammaLaw(TRNavier3DBj *rn) 
{
   read(rn->getInputFile()) ;
   write(clog) ;
}
void TGammaLaw::read(istream &in)
{

   TIOBuffer b ;
   b.read_section(in, "TGammaLaw") ;
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
void TGammaLaw::write(ostream &out) 
{
   TIOBuffer b ;
   out << "[TGammaLaw]" << endl;
   b.writeLineD(out,fGammaIndex, "fGammaIndex : Adiabatic  Index") ;
   b.writeLineD(out,fEtaOverS,   "fEtaOverS   : shear viscosity /s ") ;
   b.writeLineD(out,fSigmaOverS, "fSigmaOverS : bulk  viscosity /s ") ;
   b.writeLineD(out,fKappaTOverS,"fKappaTOverS: conductivity / s  ") ;
   b.writeLineD(out,fT0, "fT0: T0 ") ;
   b.writeLineD(out,fE0, "fE0: E0 ") ;
   b.writeLineD(out,fTPi_EtaST, "fTPi_EtaST: TauPi/(eta/s T) ") ;
   b.writeLineD(out,fL1_EtaTPi, "fL1_EtaTPi: lambda_1 /(eta * tau_pi)") ;
   b.writeLineD(out,fL2_EtaTPi, "fL2_EtaTPi: lambda_2 /(eta * tau_pi)") ;

   out << "\n" << endl;
}
std::unique_ptr<TEOS> make_eos_gammalaw(TRNavier3DBj *rn, const std::string &eosname) 
{
   return std::unique_ptr<TEOS>(new TGammaLaw(rn)) ;
}


double TGammaLaw::eofs(const double &s, const double &n) 
{
   double e =  fE0 * pow( s/fS0, fGammaIndex ) ;
   return e ;
}
   
void TGammaLaw::eos(const double &e, double &p, double &cs) 
{
   p  = (fGammaIndex-1.)*e ;
   cs = sqrt(fGammaIndex - 1.0) ;
}

void TGammaLaw::eos(const double &e, const double &n, double &p, double &cs) 
{
   eos(e, p, cs) ;
}

void TGammaLaw::st(const double &e, double &s, double &t) 
{
   double ginv = 1./fGammaIndex ;
   s  =  fS0 * pow((e/fE0),ginv) ;
   t  =  fGammaIndex * e / s ;
}

void TGammaLaw::stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) 
{
   mu = 0.0 ;
   st(e, s, t) ;
}

void TGammaLaw::viscosity(const double &e, const double &n, 
         double &sigma_overs, double &kappaT_overs, double &eta_overs) 
{
   sigma_overs  = fSigmaOverS ;
   kappaT_overs = fKappaTOverS ;
   eta_overs    = fEtaOverS ;
}




