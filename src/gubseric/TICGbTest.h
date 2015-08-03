/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RNAVIER_TICGBTEST_H
#define RNAVIER_TICGBTEST_H

#include <memory>
#include "TBjRegulator.h"
#include "TIC3D.h"
#include "TStep3D.h"
#include "THModel3D.h"

class THydro3DBj ;
class TRNavier3DBj ;

//! Simple test intial conditions for invicid Gubser evoluation
//! in 2+1 dimensions.
//!
class TICGbTest : public TIC3D
{
   private:
      THModel3D *fModel ;
      std::unique_ptr<TBjRegulator> fRegulator ;
      //! The initial time.
      double fTau0 ;
      //! Initial energy density of background Gubser flow.
      double fS0 ;
      //! Initial number density 
      double fNOverS ;
      //! Conformal parameter
      double fQTau0 ;
      //! H0=eta/e^(3/4)
      double fH0 ;
      double fTqtau ;
      double fTfact ;
        bool fIsIdealFluid ;

      void read(std::istream &in) ;
      void write(std::ostream &out) ;

   public:
      explicit TICGbTest(TRNavier3DBj *rn)  ;
      virtual ~TICGbTest() override {; } 
      virtual void nextEvent(const double &bgen) override {;}
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi) override ;

      double getCutoffC1() {return fRegulator->getCutoffC1();}

} ;

std::unique_ptr<TIC3D> make_ic_icgbtest(THydro3DBj *hy, const std::string &icname) ;

#endif

