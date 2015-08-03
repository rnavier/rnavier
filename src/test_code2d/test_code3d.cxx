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
   gRandom->SetSeed(45) ;
   hy.nextEvent() ;

   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0.4, "xy", 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0.4, "xz", 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0.4, "yz", 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0.4, "x", 0.0, 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0.4, "y", 0.0, 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0.4, "z", 0.0, 0.0) )) ;
   //analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_entropy_vs_time(&hy, 0.0) )) ;

   // Start the analysis
   cout << "Time Out = " << hy.getTime() << endl;
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   int istep = 0  ;
   double tmax =  3.0 ;
   const int nsteps_max = 10000. ;
   while(hy.getTime() < tmax && (istep < nsteps_max))  {
      hy.step() ;
      istep++ ;

         cout << "Time Out = " << hy.getTime() << endl;
         for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
            analyses[ia]->analyze()  ;
         }
      
   }
      // Stop the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->stop()  ;
   }

}
