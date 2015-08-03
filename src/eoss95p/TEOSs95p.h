#ifndef RN_TEOSLats95p_h
#define RN_TEOSLats95p_h

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "TEOS.h"
#include "numeric.h"




// Lattice EoS s95p_v0 by Huovinen and Petreczky. Data files are stored in /s95p, correspondingly
// in *par*.dat (T, \mu_1, \mu_2) and *dens*.dat (e, P, s, n_B, plasma ratio), with natural units,
// e.g., T-GeV, \mu-GeV, e-GeV/fm^3. Each file is 
// particularly available in its range.
// file1: T~[0.0159,0.1801], e~[0.001,1.001]
// file2: T~[0.1801,0.2810], e~[1.001,11.00]
// file3: T~[0.2810,0.4243], e~[11.00,61.00]
// file4: T~[0.2423,0.6335], e~[61.00,311.0]
// There are also files avaiable for T<180MeV, with 28 non-zero chemical potential for baryon, 
// strangeness, etc. For the time being, and for simplicity we ignore these files. So basically
// we have \mu = 0 and n_B = 0. 

class TEOSs95p : public virtual TEOS {

   private:

      //minimum table values
      double fEmin;
      double fSmin;
      double fTmin;
      double fPmin;
	
      gsl_interp_accel *acc_eVSp;
      gsl_spline *eVSp;
      
      gsl_interp_accel *acc_tVSx;
      gsl_spline *tVSx;
      
      gsl_interp_accel *acc_eVSt;
      gsl_spline *eVSt;
      
      gsl_interp_accel *acc_eVSs;
      gsl_spline *eVSs;
      
      gsl_interp_accel *acc_sVSe;
      gsl_spline *sVSe;

      void eos(const double &e, double &p, double &cs)  ;
      void st(const double &e, double &s, double &t)  ; 
      
      double fEtaOverS ;    //!< The Shear viscosity/entropy 
      double fEtaOverSHad ; //!< The Shear viscosity/entropy 
      double fSigmaOverS ;  //!< The Shear viscosity/entropy 
      double fKappaTOverS ; //!< The Conductivity/entropy
      double fTPi_EtaST ;  //!< tau_pi / (eta/sT)
      double fL1_EtaTPi ;  //!< lambda_1 /(eta*tau_pi)
      double fL2_EtaTPi ;  //!< lambda_2 /(eta*tau_pi)
 

   public:

      TEOSs95p ();
      ~TEOSs95p ();

      virtual void read(std::istream &in)  ;

      virtual void write(std::ostream &out) ;

      virtual void eos(const double &e, const double &n, 
            double &p, double &cs) ;

      virtual double eofs(const double &s, const double &n) ;

      virtual void stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) ;

      virtual void viscosity(const double &e, const double &n, 
         double &sigma_overs, double &kappaT_overs, double &eta_overs) ;

      virtual void getBRSSSParams(const double &e, const double &n, double &tpi_etast, double &l1_ntpi, double &l2_ntpi) 
      { tpi_etast =  fTPi_EtaST ; l1_ntpi = fL1_EtaTPi ; l2_ntpi = fL2_EtaTPi ; }

      void tmux(const double &e, const double &n, double &t, double &mub, double &mus, double &x) ;

} ;

std::unique_ptr<TEOS> make_eos_eoss95p(TRNavier3DBj *rn, const std::string &icname) ;
#endif 
