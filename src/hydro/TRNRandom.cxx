/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include "TRNRandom.h"
#ifdef RNAVIER_ROOT
#include  "TRandom3.h"
#endif 
#include "TIOBuffer.h"
#include "TRNavier2DBj.h"

#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>

unsigned long int get_random_seed() ;

gsl_rng *gGSL_RNG = 0 ;

TRNRandom::TRNRandom(TRNavier2DBj *rn) 
{
   TIOBuffer b ;
   b.read_section(rn->getInputFile(), "TRNRandom") ;
   int seed;
   int seed_root ;

   fNametag = rn->getNametag() ;
   b.getI("fSeed", seed) ;
   b.getI("fSeedROOT", seed_root) ;
   bool fUseSeedFile ;
   b.getB("fUseSeedFile", fUseSeedFile) ;

   clog << "[TRNRandom]" << endl;
   b.writeLineI(clog, seed, "fSeed ; the seed for GSL_RNG (use 0 for default)") ;
   b.writeLineI(clog, seed_root, "fSeedROOT ; the seed for gRandom (use 0 for default)") ;
   b.writeLineB(clog, fUseSeedFile, "fSeedFile ; read the random number state from nametag_gGSL_RNG.state and nametag_gRandom.state") ;

   if (!gGSL_RNG) {
      gGSL_RNG = gsl_rng_alloc(gsl_rng_default)   ;
   } else {
      gsl_rng_free(gGSL_RNG) ;
      gGSL_RNG = gsl_rng_alloc(gsl_rng_default) ;
   }

#ifdef RNAVIER_ROOT
   unsigned long int fSeedROOT = seed_root;
   if (fSeedROOT > 0 ) {
      gRandom->SetSeed(fSeedROOT)  ;
   }
   if (fUseSeedFile) {
      string  s = fNametag + "_gRandom.state" ;
      FILE *fp = fopen(s.c_str(), "r") ;
      // Use the default seed if cant find the file
      if (!fp) {
         gRandom->SetSeed() ;
         clog << "# Can't find seed file for gRandom. Using "<< gRandom->GetSeed() << endl;
      } else {
         clog << "# Reading Seed file " << s << endl;
         fclose(fp) ;
         gRandom->ReadRandom(s.c_str()) ;
      }
   }
#endif

   unsigned long int fSeed = seed;
   if (seed > 0) {
      gsl_rng_set(gGSL_RNG,fSeed) ;
   } 

   if (fUseSeedFile) {
      string  s= rn->getNametag() + "_gGSL_RNG.state" ;
      FILE *fp = fopen(s.c_str(), "r") ;
      if (!fp) {
         unsigned long int seedtry = get_random_seed() ;
         clog << "# Can't find seed file for gGSL_RNG: Using "<< seedtry << endl;
         gsl_rng_set(gGSL_RNG,seedtry);
      } else{
         clog << "# Reading Seed file " << s << endl;
         gsl_rng_fread(fp, gGSL_RNG) ;
         fclose(fp) ;
      }
   }
   clog << "# The GSL random number generator is " << gsl_rng_name(gGSL_RNG) <<endl ;
   clog << endl;
   clog << endl;

} 

TRNRandom::~TRNRandom() 
{
   {
      string  s= fNametag + "_gGSL_RNG.state" ;
      FILE *fp = fopen(s.c_str(), "w") ;
      gsl_rng_fwrite(fp, gGSL_RNG) ;
      fclose(fp) ;
   }

#ifdef RNAVIER_ROOT
   {
      string  s2 = fNametag + "_gRandom.state" ;
      gRandom->WriteRandom(s2.c_str()) ;
   }
#endif
}

unsigned long int get_random_seed()
{
   unsigned int seed;
   struct timeval tv;
   FILE *devrandom;

   if ((devrandom = fopen("/dev/random","r")) == NULL) {
      gettimeofday(&tv,0);
      seed = (tv.tv_sec + tv.tv_usec)*getpid();
   } else {
      fread(&seed,sizeof(seed),1,devrandom);
      fclose(devrandom);
   }
   return seed;
}

