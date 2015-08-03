/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include <iostream>
#include "gsl/gsl_math.h"
#include "TGlbEvent.h"
#include "TString.h"

using namespace std;

void  TGlbNucleon::Import(const TGlbNucleon &nucl)
{
   fX = nucl.fX ; fY= nucl.fY ; fZ=nucl.fZ ;  fWeight = nucl.fWeight ;
   fInNucleusA = nucl.fInNucleusA ; fNColl = nucl.fNColl ;
}

//------------------------------------------------------------------------------
ClassImp(TGlbEvent) 
//------------------------------------------------------------------------------
void TGlbEvent::fillGlbMoments(TGlbMoments &gm, const bool &shift) const 
{
   complex<double> zo(0.) ;
   double sumw=0;
   if (shift) {
      for(size_t in = 0 ; in < fN; in++) {
         if (fNucleons[in].IsSpectator()) continue ;

         double x= fNucleons[in].GetX() ;
         double y= fNucleons[in].GetY() ;
         double w= fNucleons[in].GetWeight() ;
         zo += w*complex<double>(x,y) ;
         sumw += w ;
      }
      if (sumw > 0) {
         zo /= sumw  ;
      }
   }
   complex<double> z, zc ;
   for (int i=0; i < GLB_MOMENT_DIM; i++) {
   for (int j=0; j < GLB_MOMENT_DIM; j++) {
      gm.M(i,j) = 0. ;
   }
   }
   for (size_t in = 0 ; in < fN ; in++) {
      if (fNucleons[in].IsSpectator()) continue ;

      double x= fNucleons[in].GetX() ;
      double y= fNucleons[in].GetY() ;
      double w= fNucleons[in].GetWeight() ;
      z = complex<double>(x,y)  - zo;
      zc = conj(z) ;
      for (int i=0; i < GLB_MOMENT_DIM; i++) {
      for (int j=0; j < GLB_MOMENT_DIM; j++) {
         if  (i + j < GLB_MOMENT_DIM) {
            gm.M(i,j) += w*pow(zc, i)*pow(z,j)  ;
         } else {
            gm.M(i,j) += 0.  ;
         }
      }
      }
   }
   gm.smearMoments(fSmearingSigma) ;
   gm.fillCumulants() ;
}

//------------------------------------------------------------------------------
ClassImp(TGlbMoments) 
//------------------------------------------------------------------------------

TGlbMoments::TGlbMoments() 
{
   for (int i = 0 ; i < GLB_MOMENT_DIM2 ; i++) {
      fMoments[i] = 0. ;
      fCumulants[i] = 0. ;
   }
   fMoments[0] = GSL_DBL_EPSILON ;
   fCumulants[0] = log(GSL_DBL_EPSILON) ;
   C(1,1) = GSL_DBL_EPSILON ;
}
void TGlbMoments::print(ostream &out)  const
{
   out << Form("*..............................................................................*") << endl;
   for (int i =0 ; i < 4 ; i++) {
      for(int j=0 ; j < 4 ; j++) {
         double re = Mij(i,j).real();
         double im = Mij(i,j).imag(); 
         out << Form("%+5.2e,", re)  ;
         out << Form("%+5.2e   ", im)  ;
      }
      out << endl;
   }

   out << Form("--------------------------------------------------------------------------------") << endl;
   for (int i =0 ; i < 4 ; i++) {
      for(int j=0 ; j < 4 ; j++) {
         double re = Cij(i,j).real();
         double im = Cij(i,j).imag();
         out << Form("%+5.2e,", re)  ;
         out << Form("%+5.2e   ", im)  ;
      }
      out << endl;
   }
   out << Form("*..............................................................................*") << endl;
}

void TGlbMoments::fillCumulants() 
{
#include "TGlbMoments_cumulants.cxx"
}

void TGlbMoments::smearMoments(const double &sigmar) 
{
   complex<double>Msm[GLB_MOMENT_DIM][GLB_MOMENT_DIM] ;
#include "TGlbMoments_smearMoments.cxx"
   for (int i=0; i < GLB_MOMENT_DIM ; i++) {
   for (int j=0; j < GLB_MOMENT_DIM ; j++) {
      M(i,j) = Msm[i][j];
   }
   }
   fillCumulants() ;
}

void TGlbMoments::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGlbMoments.
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      int R__i;
      for (R__i = 0; R__i < GLB_MOMENT_DIM2; R__i++) {
         double re, im ;
         R__b >> re ; 
         R__b >> im ;
         fMoments[R__i] = complex<double>(re, im)  ;
         R__b >> re ; 
         R__b >> im ;
         fCumulants[R__i] = complex<double>(re, im)  ;
      }
      R__b.CheckByteCount(R__s, R__c, TGlbMoments::IsA());
   } else {
      R__c = R__b.WriteVersion(TGlbMoments::IsA(), kTRUE);
      // Stream an object of class TGlbMoments.
      int R__i;
      for (R__i = 0; R__i < GLB_MOMENT_DIM2 ; R__i++) {
         R__b << fMoments[R__i].real() ;
         R__b << fMoments[R__i].imag() ;
         R__b << fCumulants[R__i].real() ;
         R__b << fCumulants[R__i].imag() ;
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

