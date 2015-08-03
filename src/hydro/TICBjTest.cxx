#ifndef RNAVIER_TICBJTEST_CXX
#define RNAVIER_TICBJTEST_CXX
#include "TICBjTest.h"
#include "TGrid3D.h"
#include "THModel3D.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"
#include "TStep3D.h"

TICBjTest::TICBjTest(TRNavier3DBj *rn) : 
   fModel(rn->getHModel3D()) 
{
   read(rn->getInputFile()) ;
   write(std::clog) ;
   fRegulator = unique_ptr<TBjRegulator>(new TBjRegulator(rn)) ;
}

void TICBjTest::read(std::istream &in) 
{
   TIOBuffer b ;
   b.read_section(in,"TICBjTest") ;
   b.getD("fTau0", fTau0) ;
   b.getD("fS0", fS0) ;
   b.getD("fNOverS", fNOverS) ;
   b.getD("fPiOverPiNS", fPiOverPiNS) ;

   b.getD("fSigmaxInverse", fSigmaxInverse) ;
   b.getD("fSigmayInverse", fSigmayInverse) ;
   b.getD("fSigmazInverse", fSigmazInverse) ;
}

void TICBjTest::write(std::ostream &out) 
{
   TIOBuffer b ;
   out << "# parameters for the initial condition test case " << endl;
   out << "[TICBjTest]" << endl;
   b.writeLineD(out, fTau0, "fTau0; Initial time" ) ;
   b.writeLineD(out, fS0, "fS0; Initial entropy per rapidity, s = S0/tau0") ;
   b.writeLineD(out, fNOverS, "fNOverS: Initial baryon to entropy density") ;
   b.writeLineD(out, fPiOverPiNS,"fPiOverS: Initial ratio of Pi to its navier stokes value") ;

   b.writeLineD(out, fSigmaxInverse, "fSigmaxInverse; the x-width of the gaussian which is initialized") ;
   b.writeLineD(out, fSigmayInverse, "fSigmayInverse; the y-width of the gaussian which is initialized" ) ;
   b.writeLineD(out, fSigmazInverse, "fSigmazInverse; the z-width of the gaussian which is initialized" ) ;
   out << "[END]" << endl;
   out << '\n' << endl;
}


void TICBjTest::IC(const TICPoint3D &pnt, const TNDArray1D &auxi) 
{
   // Spatial gaussian * (t0/t)**(alpha) 
   double gaus = exp( -0.5*pow(pnt.x*fSigmaxInverse, 2) 
                      - 0.5*pow(pnt.y*fSigmayInverse, 2)
                      - 0.5*pow(pnt.z*fSigmazInverse, 2)) 
                 * TBRSSS::g33(fTau0)/TBRSSS::g33(pnt.t) ;

   // Determine the initial conditions at tau0
   double s = fS0/TBRSSS::g33(fTau0)*gaus ;
   double n = fNOverS*s ;
   double e = fModel->getEOS()->eofs(s, n) ;

   // Determine the viscosity
   double sigma_overs, kappaT_overs, eta_overs ; 
   fModel->getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
   double eta = eta_overs*s ;
   //Determine the initial conditions for pizz.
   double pizz_ns = -4./3.*eta/TBRSSS::g33(fTau0) ; //navier stokes value
   double pizz   = TBRSSS::gamma33()*fPiOverPiNS * pizz_ns ;
   double u[3]   = {0., 0., 0. } ;
   double pi[5]  = {-0.5*pizz, 0., 0., -0.5*pizz, 0.}  ;

   // Regulate the initial state
   fRegulator->regulate_pi(e, n, u, pi) ; 

   // fill up the state auxi variables
   fModel->setaux(e, n, u, pi, auxi) ;
}

std::unique_ptr<TIC3D> make_ic_icbjtest(THydro3DBj *hy, const std:: string &icname)  
{
    return std::unique_ptr<TIC3D>(new TICBjTest(hy->getRNavier())) ;
}

/////////////////////////////////////////////////////////////////////////
TICStep3D::TICStep3D(TRNavier3DBj *rn, TIC3D *ic)  : fIC(ic)
{
    TIOBuffer b;
    std::ifstream &in = rn->getInputFile();
    b.read_section(in, "TICStep3D");
    b.getD("fCFL",  fCFL) ;
   // Write the inputs to log
    clog << "# Writing the stepper" << endl;
    clog << "[TICStep3D]" <<endl ;
    b.writeLineD(clog, fCFL, "fCFL: The default CFL condition used by stepper") ;
    clog << "\n" <<endl ;
    }

double TICStep3D::getDt(TGrid3D &gr)
{
    double tau = gr.getTime();
    double dx  = gr.getDx();
    double dy  = gr.getDy();
    double dz  = TBRSSS::g33(tau)*gr.getDz();
    return fCFL*GSL_MIN(dx,GSL_MIN(dy, dz));
}
int TICStep3D::step(TGrid3D &gr, const double &dt) 
{
   //Fill using initial conditions at the end of time step
   gr.updateTime(dt);
   fIC->fillGrid(&gr) ;
   //Restore time step
   gr.updateTime(-dt);
   return 1 ;
}


#endif
