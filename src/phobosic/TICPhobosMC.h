/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef rnavier_TICPhobosMC_h
#define rnavier_TICPhobosMC_h
#include <memory>
#include "TNDArray.h"
#include "TPhobosMC.h"
#include "TIC3D.h"
#include "TBjRegulator.h"

class THydro3DBj ;
class TRNavier3DBj ;

//! This class uses the TPhobosMC object to determine the initial state
//!
//! It is simply provides the necessary interface code between the
//! (stand0alone) TPhobosMC object (which does the real work) and the hyro code
//! through the TIC3D interface.
class TICPhobosMC : public TIC3D {

   private:

      TRNavier3DBj *fRNavier ;
      //! Pointer to the Glauber Object
      std::unique_ptr<TPhobosMC> fGlauberMC ;
      
      //! Pointer to the regulator class. 
      //! This class takes the intial entropy density and produces the qi and
      //! aux arrays needed for the hydro
      std::unique_ptr<TBjRegulator> fRegulator ;

   public:

      //! Read input data  from section TICGlauberMC
      TICPhobosMC(TRNavier3DBj *rn) ;
      ~TICPhobosMC()  {;}

      //! pure virtual functions of TIC3D which must be overridden
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi) override   ;
      virtual void nextEvent(const double &bgen=-1) override {fGlauberMC->nextEvent(bgen);}

      //! Returns the phobos glauber object
      TPhobosMC *getGlauberMC() { return fGlauberMC.get(); } 

} ;

std::unique_ptr<TIC3D> make_ic_icphobosmc(THydro3DBj *hydro, const std::string &name)  ;
      
#endif


