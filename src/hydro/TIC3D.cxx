#include "TIC3D.h"
#include "TGrid3D.h"
#include "THModel3D.h"

void TIC3D::fillGrid(TGrid3D *gr)
{
   TICPoint3D pnt ;
   TNDArray1D auxi ;
   TNDArray4D Aux = gr->getAux() ;
   for (int i = gr->firstX() ; i < gr->lastX() ; i++) {
   for (int j = gr->firstY() ; j < gr->lastY() ; j++) {
   for (int k = gr->firstZ() ; k < gr->lastZ() ; k++) {
      double x, y, z ;
      gr->getXYZ(i, j, k, x, y, z) ;
      double t = gr->getTime() ;
      pnt.x = x ;
      pnt.y = y ;
      pnt.z = z ;
      pnt.t = t ;
      auxi = Aux[ndarray::view(i)(j)(k)()] ;
      IC(pnt, auxi) ;
      gr->getModel()->fillaux(auxi) ;
   }
   }
   }
   gr->fillbc() ;
}


