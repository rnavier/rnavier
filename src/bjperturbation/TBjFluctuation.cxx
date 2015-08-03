#ifndef CMPOSER_TBjFluctuation_cxx
#define CMPOSER_TBjFluctuation_cxx
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "TRNavier3DBj.h"
#include "THydro3DBj.h"
#include "TBjFluctuation.h"
#include "TEOS.h"
#include "TIOBuffer.h"
#include "TIC3D.h"
#include "TGrid3D.h"
#include "TBRSSS.h"

using namespace std ;

TBjFluctuation::TBjFluctuation(THModel3D *hmodel) : fModel(hmodel),  fType(gsl_odeiv_step_rkf45)
{
    fStep = gsl_odeiv_step_alloc (fType, fSize);
    fControl = gsl_odeiv_control_y_new (fEpsAbs, fEpsRel);
    fEvolve = gsl_odeiv_evolve_alloc (fSize) ;

    gsl_odeiv_system sys={tbjfluctuation_func, 0, fSize, this} ;
    fSystem = sys ;
    
    // Default initial state corresponding
    // to a bjoken expansion (with no fluctuations) 
    // for tau0 =1, e0 = 43.8 * HBARC GeV/fm^3, 
    // No fluctuations wavenumber is zero, de, gx, gy, gz=  0
    setic(1., 
          0. , 0., 0.,
          43.8931,
          0. , 0., 0., 0.) ;
}

TBjFluctuation::~TBjFluctuation() 
{
    gsl_odeiv_evolve_free (fEvolve);
    gsl_odeiv_control_free (fControl);
    gsl_odeiv_step_free (fStep);
}

void TBjFluctuation::eos_params(const double &e0, 
       double &p0,
       double &cs0,
       double &s0,
       double &T0,
       double &eta0)
{
    double n=0.;
    //EOS calls   
    fModel->getEOS()->eos(e0, n, p0, cs0) ;
    double mu =0. ;
    fModel->getEOS()->stmu(e0, n, s0, T0, mu) ;

    // Determine the viscosities needed
    double sigma_overs, kappaT_overs, eta_overs ; 
    fModel->getEOS()->viscosity(e0, n, sigma_overs, kappaT_overs, eta_overs) ;
    eta0 = eta_overs*s0 ;
} 

void TBjFluctuation::setic(const double &t0, 
      const double &kx, const double &ky, const double &keta,
      const double &e0, 
      const std::complex<double> &de, 
      const std::complex<double> &gx, 
      const std::complex<double> &gy, 
      const std::complex<double> &gz) 
{
    gsl_odeiv_evolve_reset(fEvolve) ;

    fH = 1.e-4 ;
    fTau = t0 ;
    fKx = kx ;
    fKy = ky ;
    fKeta = keta ;

    setY(e0, de, gx, gy, gz, fY) ;
}

void TBjFluctuation::setY(const double &e0, 
      const std::complex<double> &de, 
      const std::complex<double> &gx, 
      const std::complex<double> &gy, 
      const std::complex<double> &gz, double y[]) 
{
    y[0] = e0 ;
    y[1] = real(de) ;
    y[2] = imag(de) ;

    y[3] = real(gx) ;
    y[4] = imag(gx) ;

    y[5] = real(gy) ;
    y[6] = imag(gy) ;

    y[7] = real(gz) ;
    y[8] = imag(gz) ;
}

void TBjFluctuation::getY(const double y[], double &e0, 
       std::complex<double> &de, 
       std::complex<double> &gx, 
       std::complex<double> &gy, 
       std::complex<double> &gz) 
{
    e0 = y[0] ;
    de = complex<double>(y[1],y[2]) ;
    gx = complex<double>(y[3],y[4]) ;
    gy = complex<double>(y[5],y[6]) ;
    gz = complex<double>(y[7],y[8]) ;
}

void TBjFluctuation::get_hydro_ic(const TICPoint3D &pnt, double &ehyd, double &nhyd, double *uI, double *piIJ_ns)
{
    const double &t = pnt.t;
    const double &x = pnt.x;
    const double &y = pnt.y;
    const double &eta = pnt.z; //pnt3d z means eta
    const double &z = TBRSSS::g33(pnt.t)*eta ;

    const double &kx = fKx;
    const double &ky = fKy;
    //const double &keta = fKeta;
    const double &kz = fKeta/TBRSSS::g33(pnt.t);

    double tinv = TBRSSS::gamma33()/t ;

    double e0 ;
    complex<double> de, gx, gy, gz ;
    getY(fY, e0, de, gx, gy, gz) ;

    const complex<double> I(0., 1.) ;
    complex<double> phase = exp( I * (kx * x + ky * y +  kz * z) ) ;

    double p0, cs0, s0, T0, eta0;
    eos_params(e0, p0, cs0, s0, T0, eta0) ;

    complex<double> ux = gx / (e0 + p0 + 2./3. * eta0 * tinv) ;
    complex<double> uy = gy / (e0 + p0 + 2./3. * eta0 * tinv) ;
    complex<double> uz = gz / (e0 + p0 - 4./3. * eta0 * tinv) ;

    ehyd  = real((e0 + phase * de)) ;
    nhyd  = 0. ;
    uI[0] = real( phase * ux) ;
    uI[1] = real( phase * uy) ;
    uI[2] = real( phase * uz) ;

    double p0_plus, cs0_plus, s0_plus, T0_plus, eta0_plus ;
    eos_params(ehyd, p0_plus, cs0_plus, s0_plus, T0_plus, eta0_plus) ;

    complex<double> tr = (kx * gx + ky * gy + kz * gz) ;

    piIJ_ns[0]  = real( -phase* eta0/(e0 + p0) * I *
            (kx * gx + kx * gx - 2./3. * tr)) + 2./3. * eta0_plus * tinv;

    piIJ_ns[1]  = real( -phase* eta0/(e0 + p0) * I *
            (kx * gy + ky * gx - 0.)) ;

    piIJ_ns[2]  = real( -phase* eta0/(e0 + p0) * I *
            (kx * gz + kz * gx - 0.)) ;

    piIJ_ns[3]  = real( -phase* eta0/(e0 + p0) * I *
            (ky * gy + ky * gy - 2./3. * tr )) + 2./3. * eta0_plus * tinv ;

    piIJ_ns[4]  = real( -phase* eta0/(e0 + p0) * I *
            (ky * gz + kz * gz -  0.)) ;

//    piIJ_ns[5]  = real( -phase* eta0/(e0 + p0) * I *
//            (kz * gz + kz * gz   - 2./3. * tr)) - 4./3 * eta0_plus * tinv ;

}

//! Returns the following parameters
//! 
//! tau_R = eta/[ c_s^2 (e + p)  ]
//!
//! tau_sound  =  1/(c_s k) = L/( c_s * 2 pi)
//!
//! tau_damp =  \tau_k^2/tau_R = 1/(eta/(sT) k^2) 
void TBjFluctuation::get_timescales_cartesian(double &tauR, double &tau_sound, double &tau_damp)  
{
    double e0, p0, cs0, s0, T0, eta0;
    e0 = fY[0] ;
    eos_params(e0, p0, cs0, s0, T0, eta0) ;

    tauR = eta0/(cs0*cs0*(e0 + p0)) ;

    double kz = fKeta*TBRSSS::g33(fTau) ;
    double k = sqrt(fKx*fKx + fKy*fKy + kz*kz) ;

    tau_sound = 1./GSL_MAX(cs0 * k, GSL_DBL_EPSILON) ;

    tau_damp = tau_sound * tau_sound /GSL_MAX(tauR, GSL_DBL_EPSILON) ;
    return ;
}

int TBjFluctuation::run(const double &t1, const bool &verbose) 
{
    while (fTau < t1) {
        int status = step(fTau, t1, fH, fY) ;
        if (status != GSL_SUCCESS) {
            return -1;
        }

        if (verbose) {
            double e0, p0, cs0, s0, T0, eta0;
            e0 = fY[0] ;
            eos_params(e0, p0, cs0, s0, T0, eta0) ;
            printf ("%12.4e "
                    "%15.5e %15.5e %15.5e %15.5e %15.5e" 
                    "%15.5e %15.5e %15.5e %15.5e " 
                    "%15.5e %15.5e %15.5e %15.5e \n", 
                    fTau, 
                    e0, cs0, T0, s0, eta0,
                    fY[1], fY[2], fY[3], fY[4],
                    fY[5], fY[6], fY[7], fY[8]
                    ) ;
        }
    }
    return 0 ;   
}

int tbjfluctuation_func (double t, const double y[], double f[],
        void *params)
{
    TBjFluctuation *cls = (TBjFluctuation *)params;

    double tinv = TBRSSS::gamma33()/t ;
    double g33 = TBRSSS::g33(t) ;

    double kx = cls->fKx ;
    double ky = cls->fKy ;
    double keta = cls->fKeta ; 
    double kz = keta / g33 ;

    double e0;
    complex<double> de, gx, gy, gz;
    cls->getY(y, e0, de, gx, gy, gz) ;

    double p0, cs0, s0, T0, eta0;
    cls->eos_params(e0, p0, cs0, s0, T0, eta0) ;

    double D = eta0/(e0 + p0) ; // Momentum diffusion coefficient in background

    f[0] =  - tinv * (e0 + p0 - 4./3. * eta0 * tinv ) ;

    using namespace std;
    complex<double> dp, deta ;


    dp = cs0*cs0 * de ;
    deta = eta0/s0 * 1./T0 * de ;  // This assumes eta/s constant

    const complex<double> I(0.,1.0) ;
    complex<double> dedot, gxdot, gydot, gzdot ;

    // delta Tzz = -4/3. deta *tinv - eta (dz uz + dz uz - 2/3 (d_a u^a) 
    dedot  = - tinv * (de + dp - 4./3. * deta * tinv 
             -  I * D * (2. * kz * gz - 2./3.* (kx * gx + ky * gy + kz * gz)))
             -  I * (kx * gx + ky * gy  + kz * gz) ;

    gxdot = - gx * tinv 
            - I * kx * ( dp  + 2./3. * deta * tinv )
            - D * ( kx*kx + ky*ky + kz*kz ) * gx
            - 1./3. * D * kx * ( kx * gx + ky * gy + kz * gz ) ;

    gydot = - gy * tinv
            - I * ky * ( dp  + 2./3. * deta * tinv )
            - D * ( kx*kx + ky*ky + kz*kz ) * gy
            - 1./3. * D * ky * ( kx * gx + ky * gy + kz * gz ) ;

    gzdot= - 2.* gz * tinv
             - I * kz * ( dp - 4./3 * deta * tinv )
             - D * ( kx*kx + ky*ky + kz*kz ) * gz
             - 1./3. * D * kz * ( kx * gx + ky * gy + kz * gz ) ;

    f[1] = real(dedot) ;
    f[2] = imag(dedot) ;

    f[3] = real(gxdot) ;
    f[4] = imag(gxdot) ;

    f[5] = real(gydot) ;
    f[6] = imag(gydot) ;

    f[7] = real(gzdot) ;
    f[8] = imag(gzdot) ;

    return GSL_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

TICBjFluctuation::TICBjFluctuation(TRNavier3DBj *rn, TGrid3D *gr)  : fModel(rn->getHModel3D())
{
    fTau0 = gr->getTStart() ;
    fLx = gr->Xmax() - gr->Xmin() ;
    fLy = gr->Ymax() - gr->Ymin() ;
    fLeta = gr->Zmax() - gr->Zmin() ;


    TIOBuffer b ;
    b.read_section(rn->getInputFile(),"TICBjFluctuation") ;
    b.getD("fInitialS0", fInitialS0) ;
    b.getD("fNKx", fNKx) ;
    b.getD("fNKy", fNKy) ;
    b.getD("fNKeta", fNKeta);
    b.getD("fDeOverE0", fDeOverE0);
    b.getD("fPhaseE", fPhaseE);

    // Log the inputs
    clog << "# parameters for the initial condition test case " << endl;
    clog << "[TBjFluctuation]" << endl;
    b.writeLineD(clog, fInitialS0, "fS0; Initial entropy per rapidity, s = S0/tau0") ;
    b.writeLineD(clog, fNKx, "fNKx; The wavenumber Kx= 2Pi nx/ L of the perturbation") ;
    b.writeLineD(clog, fNKy, "fNKy; The wavenumber Ky = 2Pi ny/L of the perturbation") ;
    b.writeLineD(clog, fNKeta, "fNKeta; The wavenumber Keta = 2Pi neta/L of the perturbation") ;
    b.writeLineD(clog, fDeOverE0, "fDeOverE0; Ratio of the initial perturbation"
            "to the initial backgroung") ;
    b.writeLineD(clog, fPhaseE, "fPhaseE; phase of the energy density") ;
    clog << "\n" << endl;

    double n = 0. ;
    fInitialE0 = fModel->getEOS()->eofs(fInitialS0, n) ;
    fFluctuation = unique_ptr<TBjFluctuation>(new TBjFluctuation(fModel)) ;

    nextEvent(0) ;
}

void TICBjFluctuation::nextEvent(const double &bgen)  
{
    double kx = 2.*M_PI*fNKx/fLx ;
    double ky = 2.*M_PI*fNKy/fLy ;
    double keta = 2.*M_PI*fNKeta/fLeta ;
    
    const complex<double> I(0.,1.) ;
    complex<double>de = fInitialE0 * fDeOverE0 * exp(I*fPhaseE) ;
    fFluctuation->setic(fTau0, 
            kx, ky, keta,
            fInitialE0,
            de,
            0., 0., 0.) ;
}

void TICBjFluctuation::IC(const TICPoint3D &pnt, const TNDArray1D &auxi) 
{
    using TBRSSSConsts::NPi33  ;

    double e, n;
    double uI[3] ;
    double piIJ[NPi33] ;

    fFluctuation->get_hydro_ic(pnt, e, n,  uI, piIJ) ;
    fModel->setaux(e, n, uI, piIJ, auxi) ;
    return;
}

std::unique_ptr<TIC3D> make_ic_icbjfluctuation(THydro3DBj *hy, const std::string &icname) 
{
    return unique_ptr<TIC3D>( new TICBjFluctuation(hy->getRNavier(), hy->getGrid()) ) ;
}
#endif
