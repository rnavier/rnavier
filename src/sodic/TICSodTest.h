/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RNAVIER_TICSODTEST_H
#define RNAVIER_TICSODTEST_H

#include <memory>
#include "TBjRegulator.h"
#include "TIC3D.h"
#include "TStep3D.h"
#include "THModel3D.h"

class THydro3DBj ;
class TRNavier3DBj ;

//! Simple test intial conditions for Bjorken evoluation
//! in 0+1 dimensions.
//!
//! Optionally the energy density is multiplied by a gaussian.
//! This is no longer a solution. But it provides a simple initial
//! condition to monitor the code.
class TICSodTest : public TIC3D
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

      double fPhi ; 
      double fLx ;
      double fLy ;
      double fLz ;

      void read(std::istream &in) ;
      void write(std::ostream &out) ;

   public:
      explicit TICSodTest(TRNavier3DBj *rn)  ;
      virtual ~TICSodTest() override {; } 
      virtual void nextEvent(const double &bgen) override {;}
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi) override ;

      double getS0()  {return fS0;}
      double getTau0() {return fTau0;}
      double getNOverS() {return fNOverS;}
      double getPiOverPiNS() {return fPiOverPiNS;} 
      double getCutoffC1() {return fRegulator->getCutoffC1();}

} ;

std::unique_ptr<TIC3D> make_ic_icsodtest(THydro3DBj *hy, const std::string &icname)  ;
#endif
