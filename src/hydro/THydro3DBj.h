/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RN_THYDRO3DBJ_H
#define RN_THYDRO3DBJ_H
#include <memory>
#include "TStep3D.h" 
#include "TRNavier3DBj.h"
#include "TGrid3D.h"
#include "TIC3D.h"
#include "TGrid3D.h"
#include "THModel3D.h"

class TIC3D ;
class TEOS ;
class TGrid3D ;

class THydro3DBj  {

    private:
       TRNavier3DBj *fRNavier ;

       //! Name of the stepping algorithm
       std::string fStepperName ;
       //! The actual stepping algorithm.
       std::unique_ptr<TStep3D> fStepper ;

       //! A pointer to the local grid
       std::unique_ptr<TGrid3D> fGrid ;

       //! The name of intial conditions
       std::string fICName;
       //! The 3D initializer
       std::unique_ptr<TIC3D>   fIC ;

       //! The last time step
       double fLastDt ;

       //! read the inputs for the class from the stream
       void read(std::istream &in) ;
       //! write the inputs for the class to the stream
       void write(std::ostream &out) ;
       //! An integer counting the run
       int fRun ;

   public:

       THydro3DBj(TRNavier3DBj *rn, std::unique_ptr<TIC3D> (*hydro3dbj_make_ic)(THydro3DBj *hy, const std::string &icname)) ;

       TRNavier3DBj *getRNavier() { return fRNavier; } 

       THModel3D *getHModel3D() { return fRNavier->getHModel3D() ; }
       TEOS      *getEOS() { return fRNavier->getEOS() ; }
       TGrid3D   *getGrid() { return fGrid.get() ; }

       std::string getICName() { return fICName; }
       TIC3D     *getIC() { return fIC.get() ; }

       std::string getStepperName() { return fStepperName ; }
       TStep3D   *getStep3D () { return fStepper.get() ; }

       //! Returns the current time
       double getTime() { return fGrid->getTime() ; }
       double getLastDt() { return fLastDt ; }
       double getDt() { return fStepper->getDt(*fGrid) ; }

       //! Returns the run counter.
       int getRun() { return fRun; }
       //! initialize the hydro use the current initial conditions
       void nextEvent(const double &b=-1) ;
       //! Takes a step with the default step size set by the current CFL
       void step() { step(getDt()); }
       //! Takes a step of size dt
       void step(const double &dt) ;

} ;

class bad_event: public std::exception
{
   public:
   virtual const char * what() const throw()
   {
      return "The stepper is unable to evolve state." ;
   }
} ;


#endif
