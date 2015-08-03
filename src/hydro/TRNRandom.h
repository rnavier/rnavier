#ifndef RNAVIER_TRNRandom_h
#define RNAVIER_TRNRandom_h

#include <string>
#include "gsl/gsl_rng.h"
class  TRNavier2DBj ;

//! An ultra thin  to set the seeds of the  gsl random number generator and  ROOT
//! random number generator
//!
//! [TRNRandom]
//! fSeed = 314159 
//! fSeedROOT = 628318 
//! fUseSeedFile = False 
//! 
//! On destruction this writes the state of the random number generator 
//! (gRandom or gGSL_RNG) to tape.
//!
//! nametag_gGSL_RNG.state
//! nametag_gRandom.state
//!
//! If eigher seed is negative, this initializes the seed with the default
//! seed.
class TRNRandom 
{
   private:
      std::string fNametag ;
   public: 
      TRNRandom(TRNavier2DBj *rn) ;
      ~TRNRandom() ;
} ;
extern gsl_rng *gGSL_RNG ;

#endif
