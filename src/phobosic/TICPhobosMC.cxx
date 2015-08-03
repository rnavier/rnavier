#include "TICPhobosMC.h"
#include "TPhobosMC.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"
#include "THydro3DBj.h"
#include "TEOS.h"
#include "TBjRegulator.h"

TICPhobosMC::TICPhobosMC(TRNavier3DBj *rn)  : fRNavier(rn)
{
   fGlauberMC = std::unique_ptr<TPhobosMC>(new TPhobosMC(fRNavier->getInputFile())) ;
   fRegulator = std::unique_ptr<TBjRegulator>(new TBjRegulator(fRNavier)) ;
}

void TICPhobosMC::IC(const TICPoint3D &pnt, const TNDArray1D &auxi)   
{
   //Determine the entropy and baryon number from the glauber model
   double n = 0. ;
   double s = 0. ;
   fGlauberMC->entropy(pnt.t, pnt.x, pnt.y, pnt.z, s, n) ;
   double e = fRNavier->getEOS()->eofs(s, n) ;
   double u[3] ; 
   double pi[5] ;
   // Find the bjorken pi, u
   fRegulator->find_bjorken_state(pnt, e, n, u, pi) ;

   //Put the state into the aux
   fRNavier->getHModel3D()->setaux(e,n, u, pi, auxi) ;
   return ;
}

std::unique_ptr<TIC3D> make_ic_icphobosmc(THydro3DBj *hydro, const std::string &name)  
{
   clog << "# The IC is = " << name  << endl;
   if (name == "TICPhobosMC")  {
      return unique_ptr<TIC3D>(new TICPhobosMC(hydro->getRNavier())) ;
   } else {
      clog << "** read_ic *** Unrecognized IC named " << name << endl ;
      clog << " Check configuration options (by examining configure.h) for example!" << endl;
      abort();
      return unique_ptr<TIC3D>(nullptr) ;
   }
   return unique_ptr<TIC3D>(nullptr) ; 
}
