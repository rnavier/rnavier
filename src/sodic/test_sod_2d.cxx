#include <memory>
#include <string>
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TIC3D.h"
#include "TICSodTest.h"
#include "TEOSIdeal.h"

using namespace std ;

int main(int argc, char **argv) 
{
   if  (argc < 2) {
      cout << "Usage is: test_gaus_2d.exe nametag" << endl;
   
   }

   string dir[3];
   string pln[3];
   dir[0]="x";
   dir[1]="x";
   dir[2]="y";
   pln[0]="xy";
   pln[1]="xz";
   pln[2]="yz";
   int i;
   string str1 ="sod2d_xy.ini"; 
   string str2 ="sod2d_xz.ini"; 
   string str3 ="sod2d_yz.ini"; 
   if (     argv[1]==str1) i = 0;
   else if (argv[1]==str2) i = 1;
   else if (argv[1]==str3) i = 2;
   else   {
      cout << "Wrong input file name" << argv[1] <<"!" << endl;
       
       exit(0);
   }

   TRNavier3DBj rn(argv[1], make_eos_eosideal) ;
   THydro3DBj hy(&rn, make_ic_icsodtest) ;
   hy.nextEvent() ;

   // Build up the analyses 
   vector<unique_ptr<THYAnalysis>> analyses ;

   //plots the xy plane
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_2dslice(&hy,  0.0, pln[i], 0.0) )) ;
   analyses.push_back(unique_ptr<THYAnalysis>(new  THYAnalysis_1dslice(&hy,  0.0, dir[i], 0.0, 0.0) )) ;

   // Start the analysis
   for( size_t ia=0 ;  ia < analyses.size() ; ++ia) {
      analyses[ia]->start()  ;
   }

   // Loop over time steps
   int istep = 0  ;
   double tmax =  .4 ;
   const int nsteps_max = 2000 ;
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
