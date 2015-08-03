#include <memory>
#include "TBjFluctuation.h"
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TIC3D.h"

using namespace std ;

int main(int argc, char **argv) 
{
   if  (argc < 2) {
      cout << "Usage is: test_icbj_hydro.exe nametag" << endl;
      exit(1) ;
   }

   TRNavier3DBj rn(argv[1], make_eos_gammalaw) ;
   THydro3DBj hy(&rn, make_ic_icbjfluctuation) ;
   hy.nextEvent() ;

   auto icfluct = dynamic_cast<TICBjFluctuation *>(hy.getIC()) ;
   double tauR, tau_sound, tau_damp ;
   icfluct->get_timescales_cartesian(tauR, tau_sound, tau_damp) ;
       
   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;

   double dtout = 2.0*M_PI*tau_sound/16. ;
   //plots the x axis
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  dtout, "x", 0.0, 0.0, true) )) ;
   //plots the y axis
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  dtout, "y", 0.0, 0.0, true) )) ;
   //plots the z axis
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  dtout, "z", 0.0, 0.0, true) )) ;

   // Start the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   //

   int istep = 0  ;
   double tmax =  8.* 2.*M_PI * tau_sound;
   const int nsteps_max = 20000 ;
   while(hy.getTime() < tmax && (istep < nsteps_max))  {
      istep++ ;
      cout << "Time = " << hy.getTime() << endl;
      cout << "Time Step= " << hy.getDt() << endl;
      hy.step() ;
      icfluct->evolve_to_time(hy.getTime()) ;
      for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
         analyses[ia]->analyze()  ;
      }
      istep++ ;
   }
   
   // Stop the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->stop()  ;
   }

}

