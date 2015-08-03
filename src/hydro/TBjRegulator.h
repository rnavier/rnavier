#ifndef RNAVIER_TBjRegulator_h
#define RNAVIER_TBjRegulator_h

#include "THModel3D.h"

struct TICPoint3D;

class TBjRegulator
{
   private:
      THModel3D *fModel;
      double fCutoffC1 ;
      double fEpsilon0 ;

   public:
      TBjRegulator(TRNavier3DBj *rn) ; 
      TBjRegulator(THModel3D *hm, const double &c1, const double &e0) : fModel(hm), fCutoffC1(c1), fEpsilon0(e0) { ; }

      //! Returns the primitives for a Bjorken expansion
      //!
      //! Given the energy density and number density and assuming no flow
      //! velocity and a bjorken longitudinal expansion fills up the primitives
      //! appropriately using (regulated) first order viscous hydrodynamics.
      //!
      //! Specifically, the (pre-regulated) bjorken estimate for pi^{munu}
      //! is  (pi_{xx} , p_yy , p_zz) = (+2/3 eta/t, +2/3 eta/t, -4/3 eta/t)
      void find_bjorken_state (const TICPoint3D &pnt, double &e0, double &n0, double u0[3], double pi0[5]) ;
      
      //! Applies the regulator to pi and returns a new regulated pi
      //!
      //! The regulator applies the map
      //!
      //! \pi^{\mu\nu} ->  \pi^{\mu\nu}/(1 + C1 * \pi^2 / p^2)
      void regulate_pi (const double &e0, const  double &n0, const double u[3], double pi0[5])  ;

      //! The cutoff parameter C1 used by the regulator (see regulate_pi)
      double getCutoffC1() {return fCutoffC1;}

} ;


#endif
