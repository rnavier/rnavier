#include "TBjRegulator.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"
#include "TIC3D.h"

TBjRegulator::TBjRegulator(TRNavier3DBj *rn) : fModel(rn->getHModel3D())
{
   TIOBuffer b;
   b.read_section(rn->getInputFile(), "TBjRegulator") ;
   b.getD("fCutoffC1",  fCutoffC1) ;
   b.getD("fEpsilon0",  fEpsilon0) ;
   clog << "[TBjRegulator]"  << endl;;
   b.writeLineD(clog, fCutoffC1, "fCutofffC1; The C1 parameter used to regulate the bjorken initial state") ;
   b.writeLineD(clog, fEpsilon0, "fEpsilon0; The minimum energy density used to regulate the initial state") ;
   clog << '\n' << endl;

}

void TBjRegulator::find_bjorken_state (const TICPoint3D &pnt, double &e, double &n, double u[3], double pi[5])
{
   e += fEpsilon0 ;
   // Determine the entropy
   double s,T, mu;
   fModel->getEOS()->stmu(e, n, s, T, mu) ;

   // Determine the viscosity
   double sigma_overs, kappaT_overs, eta_overs ; 
   fModel->getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
   double eta = eta_overs*s ;

   //Determine the initial conditions for pizz.
   double pizz = -4./3.*eta/pnt.t ; //navier stokes value
   u[0]=0.; u[1]=0.; u[2]=0.;
   pi[0]=-0.5*pizz; pi[1]=0.; pi[2]=0.; pi[3]=-0.5*pizz; pi[4]=0.;

   // Regulate the state
   regulate_pi(e, n, u, pi) ; 
} 

void TBjRegulator::regulate_pi (const double &e0, const  double &n0, const double u0[3], double pi0[5]) 
{
   //Create an aux view.
   double w1[THModel3D::Naux()]  ;
   TNDArray1D aux = THModel3D::view(w1) ;

   //Fill it up.
   fModel->setaux(e0, n0, u0, pi0, aux) ;

   //GetPi
   double pi[10] ;
   fModel->getPi(aux, pi) ;

   // Compute the trace
   double tr =  pi[0]*pi[0] 
               - 2*pi[1]*pi[1] - 2*pi[2]*pi[2] - 2*pi[3]*pi[3] 
               + pi[4]*pi[4] 
               + 2*pi[5]*pi[5] + 2*pi[6]*pi[6]
               + pi[7]*pi[7] 
               + 2*pi[8]*pi[8]
               + pi[9]*pi[9] ;
   
   double p = fModel->getP(aux) ;
   for (int i = 0 ; i < 5; i++) {
      pi0[i] = pi0[i]/(1. + fCutoffC1*tr/(p*p)) ;
   }
} 

