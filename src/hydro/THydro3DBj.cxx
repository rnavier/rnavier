/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include "THydro3DBj.h"
#include "TRNavier3DBj.h"
#include "TIC3D.h"
#include "TIOBuffer.h"
#include "TGrid3D.h"
#include "TBRSSSStep.h"
#include "TICBjTest.h"

using namespace std ;

//! read in the stepper
unique_ptr<TStep3D> read_stepper(THydro3DBj *hy, const string &steppername)  ;

THydro3DBj::THydro3DBj(TRNavier3DBj *rn, unique_ptr<TIC3D> (*hydro3dbj_make_ic)(THydro3DBj *hy, const string &icname)) : fRNavier(rn) 
{
   read(rn->getInputFile())  ;
   write(clog) ;

   // Construct grid
   fGrid = unique_ptr<TGrid3D>(new TGrid3D(rn)) ;

   // Read the initial conditions
   fIC = hydro3dbj_make_ic(this, getICName()) ;

   // Make the stepper. 
   fStepper = read_stepper(this, getStepperName()) ;
   
   // Set the run counter this should be -1
   // calling nextEvent) will increment the run counter 
   // setting run 0 to the first run
   fRun = -1 ;

}

void THydro3DBj::read(istream &in) 
{
   TIOBuffer buffer ;
   buffer.read_section(in, "THydro3DBj") ;
   buffer.getS("fICName",  fICName) ;
   buffer.getS("fStepperName", fStepperName) ;
}

void THydro3DBj::write(ostream &out) 
{
   TIOBuffer b ;
   out << "[THydro3DBj]" << endl;
   b.writeLineS(out,fICName, "fICName : name of the hydro initializer") ;
   b.writeLineS(out,fStepperName, "fStepperName : name of the stepper") ;
   out << "\n"  << endl;
}

//! Initializes the hydro with the current value of the initial condition
void THydro3DBj::nextEvent(const double &b) 
{
   fRun++ ;
   clog << "#  Initialization of event " << fRun << endl;
   fIC->nextEvent(b) ;
   fIC->fillGrid(getGrid()) ;
}

void THydro3DBj::step(const double &dt)
{
   if (fRun < 0) {
       cout << "** THydro3DBj::step ** nextEvent() has never been called! "
               "The grid has never been initialized." << endl;
       abort() ;
   }
   fGrid->store() ;
   fStepper->step(*getGrid(),dt) ;
   fGrid->updateTime(dt);
   fLastDt = dt ;
}

unique_ptr<TStep3D> read_stepper(THydro3DBj *hy, const string &stepper_name)  
{
   if (stepper_name == "TICStep3D") {
      return unique_ptr<TICStep3D>(new TICStep3D(hy->getRNavier(),hy->getIC()));
   }
   else if (stepper_name == "TBRSSSStep") {
      return unique_ptr<TBRSSSStep>(new TBRSSSStep(hy->getRNavier(), hy->getGrid()));
   } else {
      clog << "** read_stepper *** Unrecognized stepper named " << stepper_name <<endl ;
      abort() ;
      return unique_ptr<TStep3D>(nullptr) ;
   }
   return unique_ptr<TStep3D>(nullptr) ;
}

