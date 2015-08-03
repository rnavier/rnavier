/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
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

