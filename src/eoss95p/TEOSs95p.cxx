#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include "numeric.h"
#include "TEOSs95p.h"
#include "TIOBuffer.h"
#include "RNavierPath.h"
#include "TRNavier3DBj.h"


TEOSs95p::TEOSs95p()
{
   // Open the file
   clog << "# Reading in lattice data from file EOS_s95p_v0.dat" << endl;
   string s = GetRNAVIERDATA() + "EOS_s95p_v0.dat" ;
   FILE *LatData = fopen(s.c_str(), "r") ;
   if (!LatData) {
      clog << "** TEOSs95p ** Error in TEOSs95p ... FILE " << s << " cant be opened!" << endl;
      exit(EXIT_FAILURE) ;
      return;
   }

   // Set up storage array for the data
   static const int nC = 5 ; 
   static const int nR = 2000;

   double **LData ;
   int i,j;
   double v;
   LData = new double *[nC];
   for(i = 0; i < nC; i++) {
      LData[i] = new double [nR];
   }
   for(j = 0; j < nR; j++) {
      for(i = 0; i < nC; i++) {
         LData[i][j] = 0;
      }
   }

   // Read the data. Data stored in reverse order, due to technical
   // reason described in comine.py
   for(j = nR-1; j >=0 ; j--)
   {
      // printf("Row = %d ", j) ;
      for(i = 0; i < nC-1; i++)
      {
          fscanf(LatData, "%lf", &v);
          // printf("%f ", v) ;
          LData[i][j] = v ;
      }
      // printf("\n") ;
      fscanf(LatData, "%lf\n", &v);
      LData[i][j] = v;
   }
   fclose(LatData) ;

   //get minimum table values
   fEmin = LData[0][0];
   fPmin = LData[1][0];
   fSmin = LData[2][0];
   fTmin = LData[4][0];

   acc_eVSp = gsl_interp_accel_alloc ();
   eVSp = gsl_spline_alloc (gsl_interp_cspline, nR);

   acc_tVSx = gsl_interp_accel_alloc ();
   tVSx = gsl_spline_alloc (gsl_interp_cspline, nR);

   acc_eVSt = gsl_interp_accel_alloc ();
   eVSt = gsl_spline_alloc (gsl_interp_cspline, nR);

   acc_eVSs = gsl_interp_accel_alloc ();
   eVSs = gsl_spline_alloc (gsl_interp_cspline, nR);

   acc_sVSe = gsl_interp_accel_alloc ();
   sVSe = gsl_spline_alloc (gsl_interp_cspline, nR);

   gsl_spline_init (eVSp, LData[0], LData[1], nR);
   gsl_spline_init (tVSx, LData[4], LData[3], nR);
   gsl_spline_init (eVSt, LData[0], LData[4], nR);
   gsl_spline_init (eVSs, LData[0], LData[2], nR);
   gsl_spline_init (sVSe, LData[2], LData[0], nR);

   for (int i = 0 ; i < nC ; i++) {
      delete [] LData[i] ;
   }
   delete [] LData ;
   return ;
}

TEOSs95p::~TEOSs95p()
{
   gsl_spline_free (eVSp);
   gsl_interp_accel_free (acc_eVSp);
   gsl_spline_free (tVSx);
   gsl_interp_accel_free (acc_tVSx);
   gsl_spline_free (eVSt);
   gsl_interp_accel_free (acc_eVSt);
   gsl_spline_free (eVSs);
   gsl_interp_accel_free (acc_eVSs);
   gsl_spline_free (sVSe);
   gsl_interp_accel_free (acc_sVSe);
}

void TEOSs95p::read(istream &in)
{
   TIOBuffer b;
   
   b.read_section(in,"TEOSs95p") ;
   b.getD("fEtaOverS", fEtaOverS) ;
   b.getD("fSigmaOverS", fSigmaOverS) ;
   b.getD("fKappaTOverS", fKappaTOverS) ;
   b.getD("fTPi_EtaST", fTPi_EtaST) ;
   b.getD("fL1_EtaTPi", fL1_EtaTPi) ;
   b.getD("fL2_EtaTPi", fL2_EtaTPi) ;
}

//! Reads the input file.
//! 
//! \verbatim
//! [TEOSs95p]
//!   fEtaOverS            = 0.08         ;  shear viscosity /s 
//!   fSigmaOverS          = 1e-05        ;  bulk viscosity /s 
//!   fKappaTOverS         = 1e-05        ;  conductivity / s  
//!   fTPi_EtaST           = 2.61         ;  TauPi/(eta/s T) 
//!   fL1_EtaTPi           = 0.766        ;  lambda_1 /(eta * tau_pi)
//!   fL2_EtaTPi           = -1.06        ;  lambda_2 /(eta * tau_pi)
//! [END]
//!\endverbatim
void TEOSs95p::write(ostream &out) 
{
   TIOBuffer b;

   out << "[TEOSs95p]" << endl;
   b.writeLineD(out,fEtaOverS,   "fEtaOverS   : shear viscosity /s ") ;
   b.writeLineD(out,fSigmaOverS, "fSigmaOverS   : bulk viscosity /s ") ;
   b.writeLineD(out,fKappaTOverS,"fKappaTOverS: conductivity / s  ") ;
   b.writeLineD(out,fTPi_EtaST, "fTPi_EtaST: TauPi/(eta/s T) ") ;
   b.writeLineD(out,fL1_EtaTPi, "fL1_EtaTPi: lambda_1 /(eta * tau_pi)") ;
   b.writeLineD(out,fL2_EtaTPi, "fL2_EtaTPi: lambda_2 /(eta * tau_pi)") ;
   out << "\n" << endl;
}

// Return energy density e as a function of entropy density s and baryon density n_B=0
double TEOSs95p::eofs(const double &s0, const double &n) 
{
   double e;
   if (s0 < fSmin){
      double a = fPmin/fEmin ;
      e = fEmin * pow (s0/fSmin, 1. + a) ;
      return e;
   } 
   e = gsl_spline_eval(sVSe, s0, acc_sVSe); 
   return e;   
}


void TEOSs95p::eos(const double &e, const double &n, double &p, double &cs) 
{
   if (e < fEmin) {
       double a = fPmin/fEmin ;
       p  = a * e ;
       cs = sqrt( a ) ;
       return ;
   }
   p = gsl_spline_eval(eVSp, e, acc_eVSp) ;
   cs = sqrt( gsl_spline_eval_deriv(eVSp, e, acc_eVSp) ) ;
}


//! Returns the entropy density and temperature and chemical potential for energy density e and number density n.
void TEOSs95p::stmu(const double &e, const double &n, 
            double &s, double &t, double &mu) 
{
   if (e < fEmin)  {
      double a = fPmin/fEmin ;
      s = fSmin * pow( e/fEmin , 1./(1. + a) ) ;
      t = fTmin ;
      mu = 0. ;
      return ;
   }
   t = gsl_spline_eval(eVSt, e, acc_eVSt) ;
   s = gsl_spline_eval(eVSs, e, acc_eVSs) ;
   mu = 0.;
}

//! Returns the viscous coefficients divided by the entropy density for energy density e and number number density n.
void TEOSs95p::viscosity(const double &e, const double &n, double &sigma_overs, double &kappaT_overs, double &eta_overs) 
{
   kappaT_overs = fKappaTOverS ;
   eta_overs    = fEtaOverS ;   
   sigma_overs = fSigmaOverS ;
}


void TEOSs95p::tmux(const double &e, const double &n, double &t, double &mub, double &mus, double &x)
{
   if (e < fEmin) {
      double a = fPmin/fEmin ;
      double s = fSmin * pow( e/fEmin , 1./(1. + a) ) ;
      t = (1. + a)*e/s;
      mub= 0. ;
      mus=0. ;
      x=0. ;
   }
   t = gsl_spline_eval(eVSt, e, acc_eVSt) ;
   x = gsl_spline_eval(tVSx, t, acc_tVSx) ;
   mub = 0;
   mus = 0;
}


std::unique_ptr<TEOS> make_eos_eoss95p(TRNavier3DBj *rn, const std::string &icname) 
{
    TEOSs95p *s95p = new TEOSs95p ;
    s95p->read(rn->getInputFile()) ;
    s95p->write(clog) ;
    return std::unique_ptr<TEOS>(s95p) ;
}

