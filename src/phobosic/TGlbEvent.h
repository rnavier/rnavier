//! \addtogroup GlbEvents Stored Glauber Events
//! @{
//!
//! This set of classes is used to store Glauber Events into 
//! root trees and to analyze the output.
#ifndef RNAVIER_TMCGLAUBER_H
#define RNAVIER_TMCGLAUBER_H

#include <complex>
#include <TObject.h>

//! Describes a Glauber nucleon. This is used by TPhobosMC and TGlbEvent.
//! 
//! This class describess a nucleon. During the collision process we assign 
//! a certain amount of entropy to  a given nucleon. The entropy that each
//! nucleon carries is stored as the nucleon weight. For example
//!
//!     double x = nucleon.GetX() 
//!     
//! yields the x cordinate of the participant
class TGlbNucleon : public TObject
{
   //private:
   public:
      Double32_t fX;            //Position of nucleon
      Double32_t fY;            //Position of nucleon
      Double32_t fZ;            //Position of nucleon
      Double32_t fWeight;       //Weight associated with nucleon
      Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
      Int_t      fNColl;        //Number of binary collisions

   public:
      TGlbNucleon() : fX(0), fY(0), fZ(0), fWeight(0.), fInNucleusA(0), fNColl(0) {}
      virtual   ~TGlbNucleon() {}

      void       Collide()            {fNColl++;}
      Int_t      GetNColl()     const {return fNColl;}
      Double_t   GetX()         const {return fX;}
      Double_t   GetY()         const {return fY;}
      Double_t   GetZ()         const {return fZ;}
      Double_t   GetWeight()    const {return fWeight;}
      Bool_t     IsInNucleusA() const {return fInNucleusA;}
      Bool_t     IsInNucleusB() const {return !fInNucleusA;}
      Bool_t     IsSpectator()  const {return !fNColl;}
      Bool_t     IsWounded()    const {return fNColl;}
      void       Reset()              {fX=0;fY=0;fZ=0;fWeight=0;fInNucleusA=0;fNColl=0;}
      void       SetInNucleusA()      {fInNucleusA=1;}
      void       SetInNucleusB()      {fInNucleusA=0;}
      void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
      void       SetWeight(Double_t weight)  {fWeight = weight;}
      void       Import(const TGlbNucleon &nucl) ;
      ClassDef(TGlbNucleon,1)
};


class TGlbMoments ;

//! The purpose of this class is to provide an ultra simple class to write
//! complete Phobos Events into root Trees.  
//!
//! Just the data required to reconstruct the event is recorded in this class.
//! (The filling of this table is done by a TPhobosMC). Once the data
//! is recorded it can be used to determine epsilon2 for example. For
//! example this snippet extracts the x and y components of the first
//! (i=0) nucleon using the TGlbNucleon data structure:
//!
//!      TGlbEvent &ge = //...
//!      // i runs from 0 ... ge.getNpart(). 
//!      int i = 0 
//!      x=ge[i].GetX() ;
//!      y=ge[i].GetY() ;
//!   
//! In fact you dont need to work too hard to compute epsilon2. e For filling
//! the moments, use fillGlbMoments, which uses this class to fill up the
//! TGlbMoments object which records all moments and cumulants for a given
//! initial condition (see below)
class TGlbEvent : public TObject {

   private:

      Double_t fSmearingSigma ; ///< The smearing width
      Double_t fB_MC ; ///< The current value of B
      Int_t fN ; ///< An integer recording the number of nucleons stored=npart
      static const int fNmax=500;  ///< The maximum number of nucleons
      TGlbNucleon fNucleons[fNmax]  ; ///< A fixed array of nucleons

   public:

      TGlbEvent() : fN(0) { clear();}
      virtual ~TGlbEvent()  {;}
     
      //! Returns a reference to the i-th nucleon with i = 0 ... getNpart()
      TGlbNucleon &operator[](const size_t &i)  { return fNucleons[i]; }
      //! Returns the impact parameter of the event
      double getB() { return fB_MC ;}
      //! Returns the number of participants
      Int_t  getNpart() {return fN; }
      
      //! Fills up the moments for the current set of nucleons.
      //!
      //! 1. If shift is true then the moments are around the center of mass.
      //!    The default is to NOT include the shift.
      //!
      //! 2. After filling up the moments (with pointlike nucleons) 
      //!    this routine calls gm.smearMoments, so the returned moments
      //!    already include the gaussian smearing of pointlike nucleons
      void fillGlbMoments(TGlbMoments &gm, const bool &shift=false) const ;

      //! Erases all the nucleons stored in the event
      void clear() {
        fN=0; fB_MC =0. ;
        for (size_t i = 0 ; i< fNmax; i++) fNucleons[i].Reset() ;
      }
      //! Inserts a nucleon into the event table
      void addNucleon(const TGlbNucleon &nucl) {
         fNucleons[fN].Import(nucl) ; fN++ ;
      }
      void setSmearingSigma(const double &sigma) {fSmearingSigma=sigma;}
      double getSmearingSigma() { return fSmearingSigma;}
      void setB(const double &bgen) { fB_MC = bgen ; }



      ClassDef(TGlbEvent,1) ;

} ;

#define GLB_MOMENT_DIM 7
#define GLB_MOMENT_DIM2 (GLB_MOMENT_DIM*GLB_MOMENT_DIM)
inline size_t GLB_AT(const size_t &i,const size_t &j)
{
   //return  i*GLB_MOMENT_DIM + j ;
  return (i*GLB_MOMENT_DIM + j < GLB_MOMENT_DIM2 ? i*GLB_MOMENT_DIM + j : GLB_MOMENT_DIM2-1) ;
}

//!  Contains enough data to evaluate all moments and cumulants for
//!  a given glauber model.
//! 
//!       TGlbMoments gm  ;
//!       glauberevent.fillGlbMoments(gm) ;
//! 
//! 
//!  Then access differemnt moments or cumulants one uses
//! 
//!       Mij = g.Mij(i, j)  // Returns <z^i z^j>   
//!       Cij = g.Cij(i, j)  // Returns the cumulants <z^i z^j> - subtractions
//! 
//!  The moments are Mij are un-nomrmalized -- cumulants are always normalized.
//!  So if you want the normalized moments use
//! 
//!       normMij = g.normMij(i,j)   // which returns,  Mij(i,j)/M(0,0).
//!
//!  If you wish to set the value of the cumulants yourself, this is how you do
//!  it:
//! 
//!       for (int i = 0 ; i < GLB_MOMENT_DIM ; i++) {
//!       for (int j = 0 ; j < GLB_MOMENT_DIM ; j++) {
//!           g.M(i, j) =  \int rho(x,y) \bar z^i z^j ;
//!       }
//!       }
//!       g.fillCumulants() ;
//!  
class TGlbMoments : public TObject
{
    private:

      std::complex<double> fMoments[GLB_MOMENT_DIM2] ;
      std::complex<double> fCumulants[GLB_MOMENT_DIM2] ;

    public:

      TGlbMoments() ;
      //TGlbMoments(std::complex<double>[GLB_MOMENT_DIM2])  ;

      void setMij(const int &i, const int &j, const std::complex<double> &v)
      { fMoments[GLB_AT(i,j)] =v; }

      //! Compute the cumulants based on the moments.
      //!
      //! This code is automatically generated by a mathematica notebook
      //! see makecode_cumulants
      void fillCumulants() ;

      //! Returns a reference to the i-th and j-th moment of the
      //! distribution:
      //!
      //! Examples:
      //!
      //!     TGlbMoments g;
      //!     // ... fill up g ...
      //!     g.M(0,0) = \int rho(x,y)  = the total entropy in the event. 
      //!     g.M(0,1) = \int rho(x,y) z =  the Center of Mass shift in the event times the total entropy
      //!
      //!     g.M(1,0) = conjugate of g(0,1)
      //!
      //!     g.M(0,2)/g.M(1,1) = Complex epsilon_2  = <z^2>/<r^2>
      std::complex<double> &M(const int &i, const int &j)
      { return fMoments[GLB_AT(i,j)]; }

      //! Same as the g.M but returns the cumulants. For examle
      std::complex<double> &C(const int &i, const int &j)
      { return fCumulants[GLB_AT(i,j)]; }

      //! Same as the M(i,j) == Mij(i,j) but returns an independent complex
      //! number rather than a reference
      std::complex<double> Mij(const int &i, const int &j) const
      { return fMoments[GLB_AT(i,j)]; }

      std::complex<double> normMij(const int &i, const int &j) const
      { return Mij(i,j)/Mij(0,0).real(); }

      //! Returns  <(z\bar z)^n z^m> unnormalized
      std::complex<double> Mnm(const int &n, const int &m) const
      { return Mij(n, n+m); }

      //! Returns  <(z\bar z)^n z^m>  unnormalized
      std::complex<double> normMnm(const int &n, const int &m) const
      { return Mij(n, n+m)/Mij(0,0).real(); }

      //! Same as the C(i,j) == Cij(i,j) but returns an independent complex
      //! number rather than a reference
      std::complex<double> Cij(const int &i, const int &j) const
      { return fCumulants[GLB_AT(i,j)]; }

      //! Returns  <(z\bar z)^n z^m> - subtractions 
      std::complex<double> Cnm(const int &n, const int &m) const
      { return Cij(n,n+m); }

      void print(std::ostream &out) const ;

      TGlbMoments &operator+=(const TGlbMoments &rhs)  ;
      TGlbMoments &operator/=(const double &d) ;
      TGlbMoments &operator=(const TGlbMoments &rhs) ;

      //! Changes the moments of the point-like distribution to the moments of
      //! the smeared distribution. 
      //!
      //! After calling this function, you should call ``fillCumulants`` to
      //! recalculate the cumulants for a given set of momoments.
      //!
      //! This code is automatically generated from a mathematica notebook.
      //! see ``makecode_cumulants``
      void smearMoments(const double &sigmar) ;
      
      //! Returns epsilon0_m = Cnm(0,m)/<r^2>**m/2
      std::complex<double> epsilonM(const int &m) const
      { return  epsilon0M(m)  ; }

      //! Returns epsilon0_m = Cnm(0,m)/<r^2>**m/2
      std::complex<double> epsilon0M(const int &m) const
      { return  Cnm(0,m)/pow(Cnm(1,0).real(), 0.5*m)  ; }

      //! Returns epsilon1_m = Cnm(1,m)/<r^2>**(m/2+1)
      std::complex<double> epsilon1M(const int &m) const
      { return  Cnm(1,m)/pow(Cnm(1,0).real(), 0.5*m+1)  ; }

      //! Returns epsilon2_m = Cnm(2,m)/<r^2>**(m/2+2)
      std::complex<double> epsilon2M(const int &m) const
      { return  Cnm(2,m)/pow(Cnm(1,0).real(), 0.5*m+2.)  ; }

      //! returns <(r^2)^n> around the current origin.
      double r2(const int &n=1) const
      { return  Mnm(n,0).real()/Mnm(0,0).real() ; }

      //! Same as r2(n)
      double r2_pown(const int &n) const
      { return  Mnm(n,0).real()/Mnm(0,0).real() ; }

      ClassDef(TGlbMoments, 1)
} ;

//! Add the (un-normalized) moments of one data structure to another. 
//! 
//! The cumulants are left unchanged during this operation.
//!
//! After adding  moments one should call fillCumulants to recompute
//! the cumulants
//! 
//!     for (i = 0 ; i < nsample ; i++)
//!          m += sample[i]
//!     m /= nsample
//!     m.fillCumulants()
inline TGlbMoments &TGlbMoments::operator+=(const TGlbMoments &rhs)
{
   for (size_t i=0 ; i < GLB_MOMENT_DIM2 ; i++) {
      fMoments[i] += rhs.fMoments[i] ;
   }
   return *this ;
}
//! Normalize the moments by a factor d. Call fill cumulants after this
inline TGlbMoments &TGlbMoments::operator/=(const double &d)
{
   for (size_t i=0 ; i < GLB_MOMENT_DIM2 ; i++) {
      fMoments[i] /= d;
   }
   return *this ;
}
inline TGlbMoments &TGlbMoments::operator=(const TGlbMoments &rhs)
{
   for (size_t i=0 ; i < GLB_MOMENT_DIM2 ; i++) {
      fMoments[i] =  rhs.fMoments[i];
      fCumulants[i] = rhs.fCumulants[i];
   }
   return *this ;
}

#endif
//! @}







