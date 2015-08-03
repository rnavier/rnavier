#include <memory>
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TICBjTest.h"

using namespace std ;

int main(int argc, char **argv) 
{
   if  (argc < 2) {
      cout << "Usage is: test_icbj_hydro.exe nametag" << endl;
      exit(1) ;
   }


   TRNavier3DBj rn(argv[1], make_eos_gammalaw) ;
   THydro3DBj hy(&rn, make_ic_icbjtest) ;
   double dt = 0.001 ;
   double tmax =  3. ;
   int nsteps = (tmax - hy.getTime())/dt ; 
   hy.nextEvent() ;

   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;

   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_point(&hy,  0.0, 0.0, 0.0) )) ;

   // Start the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   for (int itime = 0 ;  itime < nsteps ; itime++) {
      hy.step(dt) ;
      cout << "# Time = " << hy.getGrid()->getTime() << endl;
      // Stop the analysis
      for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
         analyses[ia]->analyze()  ;
      }
   }
   
   // Stop the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->stop()  ;
   }

}
