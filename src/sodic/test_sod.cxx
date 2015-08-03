#include <memory>
#include <string>
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"

#include "TICSodTest.h"
#include "TEOSIdeal.h"

using namespace std ;

int main(int argc, char **argv) 
{
   if  (argc < 2) {
      cout << "Usage is: test_sod.exe nametag" << endl;
      exit(1) ;
   }

   TRNavier3DBj rn(argv[1], make_eos_eosideal) ;
   THydro3DBj hy(&rn, make_ic_icsodtest) ;
   hy.nextEvent() ;

   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;

   //plots the xy plane
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0.0, "xy", 0.0) )) ;
   //plots the x axis
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0.0, "x", 0.0, 0.0) )) ;

   // Start the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   int istep = 0  ;
   double tmax =  0.4 ;
   const int nsteps_max = 20000 ;
   while(hy.getTime() < (tmax - 1e-6) && (istep < nsteps_max))  {
      hy.step() ;
      istep++ ;
      cout << "Time = " << hy.getTime() << endl;
      cout << "Step= " << istep << endl;
   }
      for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
         analyses[ia]->analyze()  ;
      }
   
   // Stop the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->stop()  ;
   }
}
