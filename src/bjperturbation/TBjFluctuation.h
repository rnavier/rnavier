/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RNAVIER_TBjFluctuation_h
#define RNAVIER_TBjFluctuation_h

#include <gsl/gsl_odeiv.h>
#include <memory>
#include "THModel3D.h"
#include "TIC3D.h"

class THydro3DBj ;
class TRNavier3DBj ;

class TBjFluctuation{

   private:

      THModel3D *fModel ;

      double fTau ;
      double fH ;
      double fKx;
      double fKy;
      double fKeta;

      static constexpr double fEpsAbs = 1.e-10;
      static constexpr double fEpsRel = 1.e-4  ;

      static const size_t fSize= 9 ;
      double fY[fSize] ;


      const gsl_odeiv_step_type *fType ;
      gsl_odeiv_step      *fStep ;
      gsl_odeiv_control   *fControl ;
      gsl_odeiv_evolve    *fEvolve ;
      gsl_odeiv_system    fSystem  ;

      void setY(const double &e0, 
              const std::complex<double> &de, 
              const std::complex<double> &gx, 
              const std::complex<double> &gy, 
              const std::complex<double> &gz, double y[]) ;

      void getY(const double y[], double &e0, 
              std::complex<double> &de, 
              std::complex<double> &gx, 
              std::complex<double> &gy, 
              std::complex<double> &gz) ; 

      int step(double &t , const double &tf, double &h, double y[])  {
          return gsl_odeiv_evolve_apply(fEvolve, fControl, fStep,
                  &fSystem,  &t, tf, &h, y) ;
      }

      void eos_params(const double &e0, 
             double &p0,
             double &cs0,
             double &s0,
             double &T0,
             double &eta0) ;


   public:

      TBjFluctuation(THModel3D *hmodel) ;
      ~TBjFluctuation() ;

      void setic(const double &t0, 
              const double &kx, const double &ky, const double &keta,
              const double &e0, 
              const std::complex<double> &de, 
              const std::complex<double> &gx, 
              const std::complex<double> &gy, 
              const std::complex<double> &gz) ;


      void get_hydro_ic(const TICPoint3D &pnt, double &e, double &n, double *uI, double *piIJ_ns);


      int run(const double &t1, const bool &verbose = false) ;

      void get_timescales_cartesian(double &tauR, double &tau_sound, double &tau_damp) ;

      friend int tbjfluctuation_func (double t, const double y[], double f[],
         void *params) ;
} ;

int tbjfluctuation_func (double t, const double y[], double f[],
      void *params) ;

//////////////////////////////////////////////////////////////////////////////

class TICBjFluctuation : public TIC3D 
{
   private:

      THModel3D *fModel ;

      double fTau0 ;
      double fLx;
      double fLy ;
      double fLeta; 

      double fInitialS0;
      double fInitialE0;
      double fNKx;
      double fNKy;
      double fNKeta;
      double fDeOverE0;
      double fPhaseE;

      std::unique_ptr<TBjFluctuation> fFluctuation;

   public:

      TICBjFluctuation(TRNavier3DBj *rn, TGrid3D *gr) ;
      virtual void nextEvent(const double &bgen) override ;
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi) override ;
      using TIC3D::fillGrid ;

      void evolve_to_time(const double &t) 
      {
          fFluctuation->run(t) ;
      }
      void get_timescales_cartesian(double &tauR, double &tau_sound, double &tau_damp)  
      {
          fFluctuation->get_timescales_cartesian(tauR, tau_sound, tau_damp) ;
      }

};

std::unique_ptr<TIC3D> make_ic_icbjfluctuation(THydro3DBj *hy, const std::string &icname) ;

#endif //RNAVIER_TBjFluctuation_h




