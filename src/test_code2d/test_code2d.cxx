#include <memory>
#include <TRandom.h>

#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TIC3D.h"
#include "TEOSs95p.h"
#include "TICPhobosMC.h"

//#include <signal.h>
//#include <xmmintrin.h>
using namespace std ;

int main(int argc, char **argv) 
{
   if  (argc < 2) {
      cout << "Usage is: test_viscbj.exe nametag" << endl;
      exit(1) ;
   }
   // Floating point exceptions terminate
   //_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

   TRNavier3DBj rn(argv[1], make_eos_eoss95p) ;
   THydro3DBj hy(&rn, make_ic_icphobosmc) ;

   // This is the seed used in code2d.ini
   gRandom->SetSeed(2) ;
   hy.nextEvent() ;

   // The other 2d code takes a time step of 1/10 of dt to align
   // the staggered grid. We also use a CFL of 0.1 not 0.2, accounting
   // for the factor of 2.
   double dt = hy.getDt() ;
   hy.step(0.5*0.1*dt) ;

   // We use a CFL of 0.1. The other code (which is second order)
   // uses 0.2. The other code used 20 steps per output.
   const int nsteps_per_output = 10 ;
   cout << "------------------ test_code2d -------------------" << endl;
   cout << "Starting time = " << hy.getTime() << endl;
   cout << "Time step     = " << hy.getDt() << endl;
   cout << "The printing time step is = " << nsteps_per_output*hy.getDt() << endl;

   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0., "xy", 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0., "x", 0.0, 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0., "y", 0.0, 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_entropy_vs_time(&hy, 0.0) )) ;

   // Start the analysis
   cout << "Time Out = " << hy.getTime() << endl;
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   int istep = 0  ;
   double tmax =  14.5 ;
   const int nsteps_max = 10000. ;
   while(hy.getTime() < tmax && (istep < nsteps_max))  {
      hy.step() ;
      istep++ ;

      //Handle the printing
      if (istep % nsteps_per_output == 0) {
         cout << "Time Out = " << hy.getTime() << endl;
         for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
            analyses[ia]->analyze()  ;
         }
      }
   }
   
   // Stop the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->stop()  ;
   }

}
