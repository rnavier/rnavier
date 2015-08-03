#ifndef RN_TIC_h
#define RN_TIC_h

#include "TNDArray.h"

class TGrid3D ;

//! Simple struct defining a space time point in 3d.
struct TICPoint3D {
   double t ;
   double x ;
   double y ;
   double z ;
} ;

//! Base class for setting the protocol of initial conditons.
class TIC3D {
    public :
      virtual ~TIC3D() { } 
      virtual void nextEvent(const double &bgen) = 0;
      virtual void IC(const TICPoint3D &pnt, const TNDArray1D &auxi)  = 0 ;
      virtual void fillGrid(TGrid3D *grid) ;
};


#endif

