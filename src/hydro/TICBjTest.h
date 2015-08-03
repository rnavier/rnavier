#ifndef RNAVIER_TICBJTEST_H
#define RNAVIER_TICBJTEST_H

#include <memory>
#include "TBjRegulator.h"
#include "THydro3DBj.h"
#include "TIC3D.h"
#include "TStep3D.h"
#include "THModel3D.h"

class TRNavier3DBj ;

//! Simple test intial conditions for Bjorken evoluation
//! in 0+1 dimensions.
//!
//! Optionally the energy density is multiplied by a gaussian.
//! This is no longer a solution. But it provides a simple initial
//! condition to monitor the code.
class TICBjTest : public TIC3D
{
   private:
      THModel3D *fModel ;
      std::unique_ptr<TBjRegulator> fRegulator ;

      //! The initial time.
      double fTau0 ;
      //! Initial energy density of background bjorken flow.
      double fS0 ;
      //! The initial number density per entropy
      double fNOverS ;
      //! Ratio between the initial Pi^zz to the navier stokes expectation
      double fPiOverPiNS ;

      double fSigmaxInverse ; //! The width of the gaussian in x
      double fSigmayInverse ; //! The width of the gaussian in y
      double fSigmazInverse ; //! The width of the gaussian in z

      void read(std::istream &in) ;
      void write(std::ostream &out) ;

   public:
      explicit TICBjTest(TRNavier3DBj *rn)  ;
      virtual ~TICBjTest() override {; } 
      virtual void nextEvent(const double &bgen) override {;}
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi) override ;

      double getS0()  {return fS0;}
      double getTau0() {return fTau0;}
      double getNOverS() {return fNOverS;}
      double getPiOverPiNS() {return fPiOverPiNS;} 
      double getCutoffC1() {return fRegulator->getCutoffC1();}

} ;

std::unique_ptr<TIC3D> make_ic_icbjtest(THydro3DBj *hy, const std:: string &icname)  ;


//! A very simple "stepper" for testing. The "stepper" just uses the known 
//! solution to fill the  grid at each moment in time.
class TICStep3D : public TStep3D 
{
   private:
      TIC3D *fIC ;
      double fCFL; 

   public:
      TICStep3D(TRNavier3DBj *rn, TIC3D *ic);
      virtual ~TICStep3D() override {; } 
      virtual int step(TGrid3D &gr, const double &dt) override;
      double getDt(TGrid3D &gr) override;
} ;

#endif


