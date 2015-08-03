#include <memory>
#include <cstdio>
#include <boost/program_options.hpp>
#include "THYAnalysis.h"
#include "TICBjTest.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TIC3D.h"

using namespace std ;

//! Purpose:
//!
//! To simulate second order hydro for a boost invariant expansion, and 
//! to produce output that can be compared to the full hydro evolution.
//!
//! Usage
//!
//! test_icbj [--no-entropy-rewite] nametag
//!
//! Inputs:
//!
//! * The nametag. The program reads nametag.ini and generates the output Of
//! course the input file should be setup with the initial conditions
//! TICBjTest, so that the hydro solution will correspond to a (infinite in
//! transverse extent bjorken solution.
//!
//! * If the command line option is specified. Then the evolution is 
//! done without multiplying by the entropy and using ideal EOM.
//!
//! Outputs:
//!
//! On output the program produces a data file with the variables as a function of time. See the fprintf statement below for the variable list.
//!
//! The filename is nametage_icbj_srewrite.out if the entropy rewrite is used,
//! and nametag_icbj_nosrewrite.out if the entropy rewrite is not used.
//!
int main(int argc, char **argv) 
{
   bool use_entropy_rewrite = true ;
   std::string nametag;
   if (argc == 2) {
      nametag = argv[1] ;
   } else if (argc ==3) {
      use_entropy_rewrite = false ;
      std::string s = argv[1] ;
      if (s != "--no-entropy-rewrite") {
         cout << "Usage is test_icbj [--no-entropy-rewrite] nametag. " << endl;
      }
      nametag = argv[2] ;
   } else {
      cout << "Usage is test_icbj [--no-entropy-rewrite] nametag. " << endl;
   }

   TRNavier3DBj rn(nametag.c_str(), make_eos_gammalaw) ;
   THydro3DBj hy(&rn, make_ic_icbjtest) ;
   THModel3D *hm = rn.getHModel3D() ;

   // Extract the parameters from the class
   TICBjTest *ic = dynamic_cast<TICBjTest *>(hy.getIC()) ;
   double s0 = ic->getS0() ;
   double t0 = ic->getTau0() ;
   double novers = ic->getNOverS() ;
   double pioverpins = ic->getPiOverPiNS() ;
   double c1 = ic->getCutoffC1() ;

   // Set the time equal to the initial time
   double t = t0 ;

   // determine initial conditions for e and n and other thermo variables
   double s = s0/t;
   double n = novers*s ;
   double e = hm->getEOS()->eofs(s, n) ;
   double p, cs ;
   hm->getEOS()->eos(e, n, p, cs) ;
   double temper, mu ;
   hm->getEOS()->stmu(e, n, s, temper, mu) ;


   // Determine the viscosities needed
   double sigma_overs, kappaT_overs, eta_overs ; 
   hm->getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
   double eta = eta_overs*s ;
   double tpi_over_etast, l1, l2 ;
   //l1 and l2 are lambda_1/(eta*tau*pi) and lambda_2/(eta*tau_pi)
   hm->getEOS()->getBRSSSParams(e,n, tpi_over_etast, l1, l2)  ;
   double tpi = tpi_over_etast * eta/(s*temper) ;


   //Determine the initial conditions for pi.
   double pins = -4./3.*eta/t ; //navier stokes value
   double pi   = pioverpins * pins;
   double tr = 3./2.*pi*pi ;    //tr(piij*piij) 
   pi = pi/(1. + c1 * tr/(p*p)) ; // regulator used
   double spi = s * pi * t;

   // determine the nonlinear coupling
   double sig = 4./3./t ;
   tr = 3./2.*pi*sig ;
   double pi_sig = 0.5*pi*sig + 0.5*pi*sig - 1./3.*tr ; 

   // Set time stepping parameters
   double dt = 0.001 ;
   double tmax =  3. ;
   int nsteps = (tmax - hy.getTime())/dt ; 

   // TODO 
   // Determine the integration constant c0ns for the analytic
   // solution of first order viscous hydro. For speed of sound^2=1/3.
   // Also set the (analytic) expectations. Add the analytic
   // solution to the printout

   // Open the output file
   std::string name ;
   if (use_entropy_rewrite) {
      name = rn.getNametag() + "_srewrite" + "_icbj.out";
   } else {
      name = rn.getNametag() + "_nosrewrite" + "_icbj.out";
   }
   FILE *fp = fopen(name.c_str(), "w") ;

   // Forward euler
   for (int itime = 0 ;  itime < nsteps ; itime++) {
      // Print out the results and navier stokes expectation for ideal gas.
      fprintf(fp, 
              "%15.5e %15.5e %15.5e %15.5e "
              "%15.5e %15.5e %15.5e %15.5e \n", 
            t, e, n, p,
            s, temper, pins, pi) ;

      // Take a step
      e +=  -dt * (e + p + pi)/t ;

      n +=  -dt  * n / t ;

      pi += dt * (- 1./tpi*pi  + pins/tpi - 4./3.*pi/t - l1*pi_sig);

      spi +=  -dt * t * s * (1./tpi*pi  - pins/tpi + 4./3.*pi/t + l1*pi_sig);

      t += dt ;

      // determine parameters for next time step
      hm->getEOS()->eos(e, n, p, cs) ;
      hm->getEOS()->stmu(e, n, s, temper, mu) ;
      hm->getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
      eta = eta_overs*s ;
      hm->getEOS()->getBRSSSParams(e, n, tpi_over_etast, l1, l2) ; 
      tpi = tpi_over_etast * eta/(s*temper) ;
      pins = -4./3.*eta/t ; 
      sig = 4./3./t ;
      tr = 3./2.*pi*sig ;
      pi_sig = 0.5*pi*sig + 0.5*pi*sig - 1./3.*tr ; 

      // set spi or pi. Depending on which equations we want to solve
      if (use_entropy_rewrite) {
         pi = spi/(s*t);  // pi is determined by spi
      } else{
         spi = pi*s*t;  // spi is determined by pi
      }
   }
   fclose(fp) ;
   
}
