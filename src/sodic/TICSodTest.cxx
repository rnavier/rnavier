#include "TICSodTest.h"
#include "TGrid3D.h"
#include "THModel3D.h"
#include "TIOBuffer.h"
#include "THydro3DBj.h"
#include "TStep3D.h"

TICSodTest::TICSodTest(TRNavier3DBj *rn) : 
   fModel(rn->getHModel3D()) 
{
   read(rn->getInputFile()) ;
   write(std::clog) ;
   fRegulator = unique_ptr<TBjRegulator>(new TBjRegulator(rn)) ;
}

void TICSodTest::read(std::istream &in) 
{
   TIOBuffer b ;
   b.read_section(in,"TICSodTest") ;
   b.getD("fTau0", fTau0) ;
   b.getD("fS0", fS0) ;
   b.getD("fNOverS", fNOverS) ;
   b.getD("fPiOverPiNS", fPiOverPiNS) ;

   b.getD("fPhi", fPhi) ;
   b.getD("fLx", fLx) ;
   b.getD("fLy", fLy) ;
   b.getD("fLz", fLz) ;
}

void TICSodTest::write(std::ostream &out) 
{
   TIOBuffer b ;
   out << "# parameters for the initial condition test case " << endl;
   out << "[TICSodTest]" << endl;
   b.writeLineD(out, fTau0, "fTau0; Initial time" ) ;
   b.writeLineD(out, fS0, "fS0; Initial entropy per rapidity, s = S0/tau0") ;
   b.writeLineD(out, fNOverS, "fNOverS: Initial baryon to entropy density") ;
   b.writeLineD(out, fPiOverPiNS,"fPiOverS: Initial ratio of Pi to its navier stokes value") ;
   b.writeLineD(out, fPhi, "fPhi; angle of rotation around z axis" ) ;
   b.writeLineD(out, fLx, "fLx: box half-size in x" ) ;
   b.writeLineD(out, fLy, "fLy: box half-size in y" ) ;
   b.writeLineD(out, fLz, "fLz: box half-size in z" ) ;

   out << "[END]" << endl;
   out << '\n' << endl;
}

void TICSodTest::IC(const TICPoint3D &pnt, const TNDArray1D &auxi) 
{
   // Spatial gaussian * (t0/t)**(alpha) 
   //double gaus = exp( -0.5*pow(pnt.x*fSigmaxInverse, 2) 
   //                   - 0.5*pow(pnt.y*fSigmayInverse, 2)
   //                   - 0.5*pow(pnt.z*fSigmazInverse, 2)) 
   //              * fTau0/pnt.t ;
   double kc = cos(fPhi/180.*M_PI);
   double ks = sin(fPhi/180.*M_PI);
   double xr = (pnt.y*ks + kc*pnt.x);
   double yr = (pnt.y*kc - ks*pnt.x);
   double zr = pnt.z;
   double hatn = ((fabs(yr)<fLy and fabs(xr)<fLx and fabs(zr)<fLz) ? 10. : 1.);
   double hate = ((fabs(yr)<fLy and fabs(xr)<fLx and fabs(zr)<fLz) ? 30. : 1.+1e-6);


   //double hatn=1;
   //double hate=1+1e-6;
   //double u[3]   = {0., 0., 0. } ;
   //if (yr<0 and xr<0) {
   // hate = 2.0;
   // hatn = 0.5;
   //} else if (xr > 0 and yr<0){
   // hate = 1.6;
   // hatn = 0.1;
   // u[1] = 0.99/sqrt(1-0.99*0.99);
   //} else if (xr < 0 and yr > 0){
   // u[0] = 0.99/sqrt(1-0.99*0.99);
   // hate = 1.6;
   // hatn = 0.1;
   //} else if (xr >= 0 and yr >= 0){
   // hate = 2.762987e-3*3./2.;
   // hatn = 5.477875e-3;
   // hate +=hatn;
   //}


   // Determine the initial conditions at tau0
   double n = hatn ;
   double e = hate ;

   //Determine the initial conditions for pizz.
   double pi[5]  = {0, 0., 0., 0, 0.}  ;
   double u[3]   = {0., 0., 0. } ;

   // Regulate the initial state
   fRegulator->regulate_pi(e, n, u, pi) ; 

   // fill up the state auxi variables
   fModel->setaux(e, n, u, pi, auxi) ;
}

/////////////////////////////////////////////////////////////////////////

std::unique_ptr<TIC3D> make_ic_icsodtest(THydro3DBj *hy, const std::string &icname)  
{
   return unique_ptr<TIC3D>(new TICSodTest(hy->getRNavier())) ;
}
