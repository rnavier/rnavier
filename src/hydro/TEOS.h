/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
//! \addtogroup EOS Equations of State
//! @{
#ifndef RN_TGammaLaw_h
#define RN_TGammaLaw_h

#include <iostream>
#include <memory>

class TRNavier3DBj ;

//! \brief Base class for defining an Equation of State. All 
//! Hydro update rules will use only these functions
class TEOS {

   public:

      virtual ~TEOS() {} 

      //! Read the inputs for the EOS for a file
      virtual void read(std::istream &in)  = 0 ;
      //! Write the output of this EOS
      virtual void write(std::ostream &out) = 0 ;

      //! \brief Returns the pressure and the speed of sound for a
      //! given value of the energy density and number density
      virtual void eos(const double &e, const double &n, 
            double &p, double &cs) = 0 ;

      //! Returns the energy density from the entropy and number density
      virtual double eofs(const double &s, const double &n) = 0 ;

      //! \brief Returns the temperature and chemical potential for
      //! specified energy and number  densities
      virtual void stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) = 0 ;

      //! \brief Returns the viscosities relative to the entropy.
      //!
      //! \param sigma_overs The (bulk viscosity)/entropy
      //! \param kappaT_overs The (thermal conductivity*T)/entropy
      //! \param eta_overs The (shear viscosity)/entropy.
      virtual void viscosity(const double &e, const double &n, 
         double &sigma_overs, double &kappaT_overs, double &eta_overs) = 0 ;

      //! \brief Returns the second order transport coefficients.
      //!
      //! \param tpi_etast = tau_pi/(eta/(s*T) )
      //! \param l1_ntpi   = l1/(eta*tpi) 
      //! \param l2_ntpi   = l2/(eta*tpi) 
      virtual void getBRSSSParams(const double &e, const double &n, double &tpi_etast, double &l1_ntpi, double &l2_ntpi)   = 0 ;

      //! \brief Return the temperature,  baryon and strangeness
      //! chemical potentials, and the fraction of volume, \c x  ,which
      //! is in the hadron gas as a function  of energy density \c e
      //! and baryon density \c n.
      virtual void tmux(const double &e, const double &n,
            double &t, double &mub, double &mus, double &x) {;}
      
      //! Returns \f$\tau_pi\f$.
      //!
      //! This routine returns \f$\tau_\pi$ and is ment to be called
      //! after calling TEOS::stmu. The input parameters, s, temper, mu
      //! of this routine are the outputs of the TEOS::stmu function call.
      virtual double getBRSSSTauPi(const double & e, 
            const double &n, 
            const double &s,
            const double &temper, 
            const double &mu) ; 
} ;


//! \brief This class describes  polytropic equation of state.
//!
//! The equation of state is :
//!
//!    p = (fGammaIndex -1) * e  and 
//!    s = fS0 * (e/fE0)^(1/fGammaIndex) 
//!
//! The formula for the entropy follows from the second law
//! thermodynamics. 
class TGammaLaw : public virtual TEOS{

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

      TGammaLaw(TRNavier3DBj *rn) ;

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
      { tpi_etast =  fTPi_EtaST ; l1_ntpi = fL1_EtaTPi ; l2_ntpi = fL2_EtaTPi ; }

};

std::unique_ptr<TEOS> make_eos_gammalaw(TRNavier3DBj *rn, const std::string &eosname) ;

#endif
//! @}
