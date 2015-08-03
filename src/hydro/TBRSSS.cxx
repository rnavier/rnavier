/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef rnavier_TBRSSS_CXX
#define rnavier_TBRSSS_CXX
#include "TRNavier3DBj.h"
#include "TEOS.h"
#include "numeric.h"
#include "TIOBuffer.h"
#include "TBRSSS.h"

bool TBRSSS::fIsCartesian = false ;

TBRSSS::TBRSSS(TRNavier3DBj *rn) : fEOS(rn->getEOS())
{
    TIOBuffer b;
    b.read_section(rn->getInputFile(), "TBRSSS") ;
    b.getB("fIsCartesian", fIsCartesian) ;
//    b.getD("fE0Regulator", fE0Regulator) ;
    
    clog << "# TBRSSS parameters " << endl;
    clog << "[TBRSSS]" << endl;
    b.writeLineB(clog, fIsCartesian, "fIsCartesian; Selects wether the model assumes cartesian coordinates") ;
    clog << "\n" << endl;

}

void TBRSSS::fillaux(const TNDArray1D &aux)
{
    aux[10] = evaluateP(aux);
    double ux = getU1(aux) ;
    double uy = getU2(aux) ;
    double uz = getU3(aux) ;
    aux[11] = sqrt(1. + ux*ux + uy*uy + uz*uz) ;
}
double TBRSSS::evaluateP(const TNDArray1D &aux)
{
    double e = getE(aux) ;
    double n = getN(aux) ;
    double p, cs ;
    fEOS->eos(e, n, p, cs) ;
    return p;
}
double TBRSSS::getCs(const TNDArray1D &aux)
{
    double e = getE(aux) ;
    double n = getN(aux) ;
    double p, cs ;
    fEOS->eos(e, n, p, cs) ;
    return cs ;
}
double TBRSSS::getS(const TNDArray1D &aux)
{
    double e = getE(aux) ;
    double n = getN(aux) ;
    double s, T, mu ;
    fEOS->stmu(e,n,s,T,mu) ;
    return s;
}

double TBRSSS::getT(const TNDArray1D &aux)
{
    double e = getE(aux) ;
    double n = getN(aux) ;
    double s, T, mu ;
    fEOS->stmu(e,n,s,T,mu) ;
    return T;
}
void TBRSSS::getIMVSC(const TNDArray1D &aux, double *vc)
{
    double s, T, mu ;
    double e = getE(aux) ;
    double n = getN(aux) ;
    getEOS()->stmu(e, n, s, T, mu) ;
    double taupi= getEOS()->getBRSSSTauPi(e, n, s, T, mu) ;
    vc[0] = 1./GSL_MAX(taupi, GSL_DBL_MIN) ;
}
void TBRSSS::getEXVSC(const TNDArray1D &aux, double *vc)
{
    double e = getE(aux)  ;
    double n = getN(aux)  ;
    double s, T, mu ;
    getEOS()->stmu(e, n, s, T, mu) ;
    double sigma_overs, kappaT_overs, eta_overs;
    getEOS()->viscosity(e, n, sigma_overs, kappaT_overs, eta_overs);
    double eta= s*eta_overs ;
    double taupi_over_etast, l1, l2 ;
    getEOS()->getBRSSSParams(e,n, taupi_over_etast, l1, l2)  ;
    double taupi = taupi_over_etast * eta/(s*T) ;

    // Regulator
    const double e0 = 0. ;
    double reg = e/sqrt(e*e + e0*e0) ;

    vc[0] = reg*eta/GSL_MAX(taupi, GSL_DBL_MIN) ;
    vc[1] = reg*4./3 ;
    vc[2] = reg*l1 ;
    vc[3] = reg*l2 ;
    vc[4] = reg/GSL_MAX(taupi, GSL_DBL_MIN);
    vc[5] = -reg*eta;
    return ;
}

void TBRSSS::getPi(const TNDArray1D &aux, double *pi) const
{
    const double ut = getU0(aux) ;
    const double ux = getU1(aux) ;
    const double uy = getU2(aux) ;
    const double uz = getU3(aux) ;

    pi[1] = getPi01(aux) ; //! pitx
    pi[2] = getPi02(aux) ; //! pity

    pi[4] = getPi11(aux) ; //! pixx
    pi[5] = getPi12(aux) ; //! pixy
    pi[6] = getPi13(aux) ; //! pixz
    pi[7] = getPi22(aux) ; //! piyy
    pi[8] = getPi23(aux) ; //! piyz
    pi[9] = getPi33(aux) ; //! pizz

    pi[0] = pi[4]+pi[7]+pi[9] ; //! pitt
    pi[3] = (ux*pi[6] + uy*pi[8] + uz*pi[9])/ut ; //! pitz
}

void TBRSSS::getU(const TNDArray1D &aux, double *u) const
{
    u[0] = getU0(aux) ;
    u[1] = getU1(aux) ;
    u[2] = getU2(aux) ;
    u[3] = getU3(aux) ;
}
void TBRSSS::setauxIdeal(const double &e0, const double &n0, const double uI[3], const TNDArray1D &aux)
{
    aux[0] = e0 ;
    aux[9] = n0 ;
    aux[1] = uI[0] ;
    aux[2] = uI[1] ;
    aux[3] = uI[2] ;
}
void TBRSSS::setaux(const double &e0, const double &n0, const double uI[3], const double piIJ[5], const TNDArray1D &aux)
{
    aux[0] = e0 ;
    aux[9] = n0 ;
    aux[1] = uI[0] ;
    aux[2] = uI[1] ;
    aux[3] = uI[2] ;
    aux[4] = piIJ[0] ;
    aux[5] = piIJ[1] ;
    aux[6] = piIJ[2] ;
    aux[7] = piIJ[3] ;
    aux[8] = piIJ[4] ;
    fillaux(aux) ;
}



//Checked DT1
void TBRSSS::fluxX(const double &tau, const TNDArray1D &aux, const TNDArray1D &fx)
{
    double e    = getE(aux) ;
    double n    = getN(aux) ;
    double p    = getP(aux) ;
    double s    = getS(aux) ;
    double ut   = getU0(aux) ;
    double ux   = getU1(aux);
    double uy   = getU2(aux) ;
    double uz   = getU3(aux) ;
    double pi10 = getPi01(aux) ;
    double pi11 = getPi11(aux) ;
    double pi12 = getPi12(aux) ;
    double pi13 = getPi13(aux) ;
    double pi22 = getPi22(aux) ;
    double pi23 = getPi23(aux) ;

    fx[0] = ((e + p) * ut * ux  +       pi10) * g33(tau) ;
    fx[1] = ((e + p) * ux * ux  +  p  + pi11) * g33(tau) ;
    fx[2] = ((e + p) * uy * ux  +       pi12) * g33(tau) ;
    fx[3] = ((e + p) * uz * ux  +       pi13) * g33(tau) ;
    fx[4] = ux * pi11 * s * g33(tau) ;
    fx[5] = ux * pi12 * s * g33(tau) ;
    fx[6] = ux * pi13 * s * g33(tau) ;
    fx[7] = ux * pi22 * s * g33(tau) ;
    fx[8] = ux * pi23 * s * g33(tau) ;
    fx[9] = ux * n * g33(tau) ;

}

// Checked DT1
void TBRSSS::fluxY(const double &tau, const TNDArray1D &aux, const TNDArray1D &fy)
{
    double e    = getE(aux) ;
    double n    = getN(aux) ;
    double p    = getP(aux) ;
    double s    = getS(aux) ;
    double ut   = getU0(aux) ;
    double ux   = getU1(aux);
    double uy   = getU2(aux) ;
    double uz   = getU3(aux) ;
    double pi20 = getPi02(aux) ;
    double pi11 = getPi11(aux) ;
    double pi12 = getPi12(aux) ;
    double pi13 = getPi13(aux) ;
    double pi22 = getPi22(aux) ;
    double pi23 = getPi23(aux) ;


    fy[0] = ((e + p) * ut * uy  +       pi20) * g33(tau) ;
    fy[1] = ((e + p) * ux * uy  +       pi12) * g33(tau) ;
    fy[2] = ((e + p) * uy * uy  +   p + pi22) * g33(tau) ;
    fy[3] = ((e + p) * uz * uy  +       pi23) * g33(tau) ;
    fy[4] = uy * pi11 * s * g33(tau) ;
    fy[5] = uy * pi12 * s * g33(tau) ;
    fy[6] = uy * pi13 * s * g33(tau) ;
    fy[7] = uy * pi22 * s * g33(tau) ;
    fy[8] = uy * pi23 * s * g33(tau) ;
    fy[9] = uy * n * g33(tau) ;

}

// Checked DT1
void TBRSSS::fluxZ(const double &tau, const TNDArray1D &aux, const TNDArray1D &fz)
{
    double e    = getE(aux) ;
    double n    = getN(aux) ;
    double p    = getP(aux) ;
    double s    = getS(aux) ;
    double ut   = getU0(aux) ;
    double ux   = getU1(aux);
    double uy   = getU2(aux) ;
    double uz   = getU3(aux) ;
    double pi30 = getPi03(aux) ;
    double pi11 = getPi11(aux) ;
    double pi12 = getPi12(aux) ;
    double pi13 = getPi13(aux) ;
    double pi22 = getPi22(aux) ;
    double pi23 = getPi23(aux) ;
    double pi33 = getPi33(aux) ;


    fz[0] = ((e + p) * ut * uz  +       pi30) * g33(tau) ;
    fz[1] = ((e + p) * ux * uz  +       pi13) * g33(tau) ;
    fz[2] = ((e + p) * uy * uz  +       pi23) * g33(tau) ;
    fz[3] = ((e + p) * uz * uz  +  p  + pi33) * g33(tau) ;
    fz[4] = uz * pi11 * s * g33(tau) ;
    fz[5] = uz * pi12 * s * g33(tau) ;
    fz[6] = uz * pi13 * s * g33(tau) ;
    fz[7] = uz * pi22 * s * g33(tau) ;
    fz[8] = uz * pi23 * s * g33(tau) ;
    fz[9] = uz * n * g33(tau) ;

}

//checked DT1
void TBRSSS::charges(const double &tau, const TNDArray1D &aux, const TNDArray1D &qi)
{
    double e    = getE(aux) ;
    double n    = getN(aux) ;
    double p    = getP(aux) ;
    double s    = getS(aux) ;
    double ut   = getU0(aux) ;
    double ux   = getU1(aux);
    double uy   = getU2(aux) ;
    double uz   = getU3(aux) ;
    double pi00 = getPi00(aux) ;
    double pi01 = getPi01(aux) ;
    double pi02 = getPi02(aux) ;
    double pi03 = getPi03(aux) ;
    double pi11 = getPi11(aux) ;
    double pi12 = getPi12(aux) ;
    double pi13 = getPi13(aux) ;
    double pi22 = getPi22(aux) ;
    double pi23 = getPi23(aux) ;


    qi[0] = ((e + p)*ut*ut  - p + pi00) * g33(tau) ;
    qi[1] = ((e + p)*ut*ux      + pi01) * g33(tau) ;
    qi[2] = ((e + p)*ut*uy      + pi02) * g33(tau) ;
    qi[3] = ((e + p)*ut*uz      + pi03) * g33(tau) ;
    qi[4] = ut * pi11 * s * g33(tau) ;
    qi[5] = ut * pi12 * s * g33(tau) ;
    qi[6] = ut * pi13 * s * g33(tau) ;
    qi[7] = ut * pi22 * s * g33(tau) ;
    qi[8] = ut * pi23 * s * g33(tau) ;
    qi[9] = ut * n * g33(tau) ;

}

void TBRSSS::getTtI(const double &tau, const TNDArray1D &aux, double TtI[3])
{
    double e    = getE(aux) ;
    double p    = getP(aux) ;
    double ut   = getU0(aux) ;
    double ux   = getU1(aux);
    double uy   = getU2(aux) ;
    double uz   = getU3(aux) ;
    double pi01 = getPi01(aux) ;
    double pi02 = getPi02(aux) ;
    double pi03 = getPi03(aux) ;

    TtI[0] = ((e + p)*ut*ux      + pi01) * g33(tau) ;
    TtI[1] = ((e + p)*ut*uy      + pi02) * g33(tau) ;
    TtI[2] = ((e + p)*ut*uz      + pi03) * g33(tau) ;

}

#endif

