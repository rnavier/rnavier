//! \addtogroup EOS Equations of State
//! @{
#ifndef RN_TEOSIdeal_h
#define RN_TEOSIdeal_h

#include <iostream>
#include <TEOS.h>


//! \brief This class describes  polytropic equation of state.
//!
//! The equation of state is :
//!
//!    p = (fGammaIndex -1) * (e-n)  and 
//!    s = fS0 * ((e-(1-1/fGammaIndex)*n)/fE0)^(1/fGammaIndex) 
//!
//! The formula for the entropy follows from the second law
//! thermodynamics. 
class TEOSIdeal : public virtual TEOS{

   private:

      static const double fDim ;
      double fGammaIndex ;  //!< The Gamma index 
      double fEtaOverS ;    //!< The Shear viscosity/entropy 
      double fSigmaOverS ;  //!< The Bulk viscosity/entropy
      double fKappaTOverS ; //!< The Conductivity/entropy

      double fT0 ;          //!< The temperature at energy density fE0
      double fE0 ;          //!< The reference energy density

      double fDOF ;          //!< The entropy density at fE0.
      double fS0 ;          //!< The entropy density at fE0.

      double fTPi_EtaST ;  //!< tau_pi / (eta/sT)
      double fL1_EtaTPi ;  //!< lambda_1 /(eta*tau_pi)
      double fL2_EtaTPi ;  //!< lambda_2 /(eta*tau_pi)
   
      void eos(const double &e, double &p, double &cs)  ;
      void st(const double &e, double &s, double &t)  ;

   public:

      double gamma() { return fGammaIndex ; }

      void read(std::istream &in)  ;
      void write(std::ostream &out) ;

      void eos(const double &e, const double &n, 
            double &p, double &cs) ;

      double eofs(const double &s, const double &n) ;

      void stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) ;

      void viscosity(const double &e, const double &n, 
         double &sigma_overs, double &kappaT_overs, double &eta_overs) ;

      virtual void getBRSSSParams(const double &e, const double &n, double &tpi_etast, double &l1_ntpi, double &l2_ntpi) 
      { tpi_etast =  fTPi_EtaST ; l1_ntpi = fL1_EtaTPi ; l2_ntpi = fL2_EtaTPi; }

};


std::unique_ptr<TEOS> make_eos_eosideal(TRNavier3DBj *rnavier, const std::string &icname) ;


#endif
//! @}
