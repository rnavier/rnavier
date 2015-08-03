/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include "TICGbTest.h"
#include "THydro3DBj.h"
#include "TGrid3D.h"
#include "THModel3D.h"
#include "TIOBuffer.h"
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"
#include "TStep3D.h"
#include "gsl/gsl_sf_hyperg.h"

double Htemp(double const &q, double const &tau,double const &x, double const &y);
TICGbTest::TICGbTest(TRNavier3DBj *rn) : 
   fModel(rn->getHModel3D()) 
{
   read(rn->getInputFile()) ;

   write(std::clog) ;
   fRegulator = unique_ptr<TBjRegulator>(new TBjRegulator(rn)) ;

}

void TICGbTest::read(std::istream &in) 
{
   TIOBuffer b ;
   b.read_section(in,"TICGbTest") ;
   b.getD("fTau0", fTau0) ;
   b.getD("fS0", fS0) ;
   b.getD("fNOverS", fNOverS) ;
   b.getD("fQTau0", fQTau0) ;
   b.getB("fIsIdealFluid", fIsIdealFluid) ;
   //Numerically determine H0
   if (fIsIdealFluid) {
       fH0 = 0;
   } else {
       double s,mu;
       fModel->getEOS()->stmu(1.0,0.0,s,fTfact,mu);
       // Determine the viscosity
       double sigma_overs, kappaT_overs, eta_overs ; 
       fModel->getEOS()->viscosity(1.0, 0.0, sigma_overs, kappaT_overs, eta_overs) ;
       fH0 = eta_overs*s ;
   }
      
   //Determine conformal temperature for given entropy density (fS0=s/tau) at initial time
   // we define e=(T_hat/tau)^4
   double T_0 = fTau0*pow(fModel->getEOS()->eofs(fS0*fTau0, 0.0),1./4.) ;
   //Determine \hat{T}_0*(2*q*tau0)^{2/3}
   fTqtau = pow(1.0 + fQTau0*fQTau0, 2./3.)*(T_0 - fH0*Htemp(fQTau0/fTau0, fTau0, 0.0, 0.0));
}

void TICGbTest::write(std::ostream &out) 
{
   TIOBuffer b ;
   out << "# parameters for the initial condition test case " << endl;
   out << "[TICGbTest]" << endl;
   b.writeLineD(out, fTau0, "fTau0; Initial time" ) ;
   b.writeLineD(out, fS0, "fS0; Initial entropy per rapidity, s = S0/tau0") ;
   b.writeLineD(out, fNOverS, "fNOverS: Initial baryon to entropy density") ;
   b.writeLineD(out, fQTau0, "fQTau0: dimensionless gubser parameter") ;
    b.writeLineB(clog, fIsIdealFluid, "fIsIdealFluid: The update step is assumed to be ideal") ;
   out << "[END]" << endl;
   out << '\n' << endl;
}

void TICGbTest::IC(const TICPoint3D &pnt, const TNDArray1D &auxi) 
{

   // Determine the initial conditions at tau0
   
   double q = fQTau0/fTau0;
   double tau = pnt.t;
   double x = pnt.x;
   double y = pnt.y;
   double r = sqrt(x*x+y*y);
   // cos(phi) an sin(phi) for azimuthal angle phi
   double cs = x/r;
   double ss = y/r;
   // sinh(kappa) and cosh(kappa) where tanh(kappa) = v_radial
   double kappa = atanh((2*q*q*tau*r)/(1.0+q*q*(tau*tau+r*r)));
   double shkappa = sinh(kappa);
   double chkappa = cosh(kappa);

   //Set gubser flow velocity in Bjorken coordinates
   double u[3]   = {shkappa*cs, shkappa*ss, 0. } ;
   double denom = 4*q*q*tau*tau + pow(1.0 - q*q*(tau*tau-r*r), 2);
   double T_hat = fTqtau*pow(tau/fTau0, 2./3.)/pow(denom, 1./3.) +
        fH0*Htemp(q,tau,x,y);
   double n = fNOverS*4./3.*pow(T_hat/tau, 3)/fTfact ;
   double e = pow(T_hat/tau, 4) ;

   //Determine the initial conditions for pizz.
   double pi[5]  = {0, 0., 0., 0, 0.}  ;
   if (fIsIdealFluid==false){
       double s,T,mu;
       fModel->getEOS()->stmu(e,n,s,T,mu);
       // Determine the viscosity
       double sigma_overs, kappaT_overs, eta_overs ; 
       fModel->getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs) ;
       double eta = eta_overs*s ;
       double thrho = -(1-q*q*(tau*tau-r*r))/sqrt(denom);
       //Determine the initial conditions for pizz.
       double pizz   =  TBRSSS::gamma33()*4./3.*eta*thrho/TBRSSS::g33(tau);
       double kxx = cs*cs*chkappa*chkappa+ss*ss;
       double kxy = cs*ss*chkappa*chkappa-ss*cs;
       double kyy = ss*ss*chkappa*chkappa+cs*cs;
       pi[0] = -0.5*kxx*pizz;
       pi[1] = -0.5*kxy*pizz;
       pi[2] =  0.;
       pi[3] = -0.5*kyy*pizz;
       pi[4] = 0. ;
   }
   // Regulate the initial state
   fRegulator->regulate_pi(e, n, u, pi) ; 

   // fill up the state auxi variables
   fModel->setaux(e, n, u, pi, auxi) ;
}


// Simple extension of hypergeometric function 2F1(0.5,1.0/6.0,1.5,-g^2) for g>1.
// It uses Taylor expansion around g=1.0 and g=Infinty calcuated with Mathematica 
// Absolute error < 10-9 (around g*1.8)
double Htemp(double const &q, double const &tau,double const &x, double const &y)
{
    if ((tau*q)==0) return -0.5;
    
    double gr = (1.0 - q*q*(tau*tau-x*x-y*y))/(2*q*tau);
    double g = fabs(gr);
    if (g < 1.0){
        return gr/sqrt(1.0+g*g)*(1.0 
                - pow(1.0 + g*g, 1./6.)*gsl_sf_hyperg_2F1(0.5,1.0/6.0,1.5,-g*g));
    } else if (g < 1.8){
        return gr/sqrt(1.0+g*g)*(1.0 
                - pow(1.0 + g*g, 1./6.)*(
         1.0000001718854175 + g*(-4.521000220209653e-6 + 
      g*(-0.0554984843446401 + 
         g*(-0.0004600184252026192 + 
            g*(0.022101393692005825 + 
               g*(-0.011695988075641858 + 
                  g*(0.030711261530574088 + 
                     g*(-0.11507649403981075 + 
                        g*(0.27402379577481983 + 
                           g*(-0.5191684539611702 + 
                           g*(0.8390308038802473 + 
                           g*(-1.150818608221189 + 
                           g*(1.3224583993675862 + 
                           g*(-1.2673876174647964 + 
                           g*(1.0136487225446893 + 
                           g*(-0.6776505361624287 + 
                           g*(0.3786827473381469 + 
                           g*(-0.17639972245638566 + 
                           g*(0.06806383987191858 + 
                           g*(-0.021519077516320054 + 
                           g*(0.005481943546668366 + 
                           g*(-0.0010971337268348008 + 
                           g*(0.00016592207347274068 + 
                           g*(-0.000017792769388513818 + 
                           (1.2023396143423515e-6 - 3.833164949340606e-8*g)*
                           g))))))))))))))))))))))) 
                           ));
    } else{
        return gr/sqrt(1.0+g*g)*(1.0 
                - pow(1.0 + g*g, 1./6.)*(
1.110126802162546e-21*pow(1/g,24.333333333333332)*
    (-8.694186270330961e17 + g*g*
       (1.021891296699348e18 + 
         g*g*(-1.220030762368018e18 + 
            g*g*(1.484512955608637e18 + 
               g*g*(-1.849384374156634e18 + 
                  g*g*(2.374093429149911e18 + 
                     g*g*
                      (-3.170490430025318e18 + 
                        g*g*
                         (4.470829822985011e18 + 
                           g*g*
                           (-6.828176456922563e18 + 
                           g*g*
                           (1.185946437254971e19 + 
                           g*g*
                           (-2.62732749176486e19 + 
                           g*g*
                           (1.125997496470654e20 + 
                           1.351196995764785e21*g*g)))))))))))) + 
   (sqrt(M_PI)*tgamma(-0.3333333333333333))/(2.*g*tgamma(0.16666666666666666))
    ));
    }

}


std::unique_ptr<TIC3D> make_ic_icgbtest(THydro3DBj *hy, const std::string &icname) 
{
   return unique_ptr<TIC3D>(new TICGbTest(hy->getRNavier())) ;
}

