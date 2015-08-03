#include "TGrid3D.h"
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"

using namespace std ;

//! Allocates 4D arrays for current and past values of primitive variables
TGrid3D::TGrid3D(TRNavier3DBj *rn) 
   : fModel(rn->getHModel3D()),
     fAuxs(fNhistory),
     fTimes(fNhistory) 
{ 
   read(rn->getInputFile()) ; // Read the inputs from 
   write(clog) ;
   resetTime(getDomainTime()) ;
   for (int i = 0 ; i < Nhistory() ; i++) {
     fTimes[i] = getTime() ;
     fAuxs[i]  = makeGrid(fModel->Naux()) ;
     fAuxs[i].deep()  = 0. ;
   }
}

void TGrid3D::read(istream &in) 
{
   TDomain3D::read(in) ;
   TIOBuffer b;
   b.read_section(in, "TGrid3D")  ;
   b.getS("fBCName", fBCName) ;
}
void TGrid3D::write(ostream &out) const
{
   TDomain3D::write(out) ;

   out << "[TGrid3D]" << endl;
   TIOBuffer b ;
   b.writeLineS(out, fBCName, "fBCName;  acceptable choices are periodic/copy/periodicz") ;
   out << "\n" << endl;
}


void TGrid3D::store() 
{
   for (int i = Nhistory()-1 ; i > 0 ; i--) {
      fAuxs[i].deep() = fAuxs[i-1] ;
      fTimes[i] = fTimes[i-1] ;
   }
}

void TGrid3D::fillbc() 
{
   if (fBCName == "periodic") {
      fillbc_periodic();
   } else if (fBCName == "copy") {
      fillbc_copy();
   } else if (fBCName == "periodicz") {
      fillbc_periodicz() ;
   } else {
      cerr << "** TGrid3D::fillbc() ** Unrecognized boundary condition named " << fBCName << ".  Acceptable choices are \"periodic\", \"periodicz\", or \"copy\"."<< endl;
   }
}
void TGrid3D::fillbc_periodicz() 
{
   using namespace ndarray ;
   TNDArray4D Aux = getAux() ;
   for (int i=0 ; i < NbndryX()  ; i++) {
       int i1 = firstX() - 1 - i;
       int i2 = firstX();
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
       i1 = lastX() + i;
       i2 = lastX() - 1;
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
   }
   for (int j=0 ; j < NbndryY()  ; j++) {
       int j1 = firstY() - 1 - j;
       int j2 = firstY();
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
       j1 = lastY() + j;
       j2 = lastY() - 1;
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
   }
   for (int k=0 ; k < NbndryZ()  ; k++) {
       int k1 = firstZ() - 1 - k;
       int k2 = lastZ() - 1 - k;
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
       k1 = lastZ() + k;
       k2 = firstZ() + k;
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
   }
}

void TGrid3D::fillbc_periodic() 
{
   using namespace ndarray ;
   TNDArray4D Aux = getAux() ;
   for (int i=0 ; i < NbndryX()  ; i++) {
       int i1 = firstX() - 1 - i;
       int i2 = lastX() - 1 - i;
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
       i1 = lastX() + i;
       i2 = firstX() + i;
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
   }
   for (int j=0 ; j < NbndryY()  ; j++) {
       int j1 = firstY() - 1 - j;
       int j2 = lastY() - 1 - j;
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
       j1 = lastY() + j;
       j2 = firstY() + j;
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
   }
   for (int k=0 ; k < NbndryZ()  ; k++) {
       int k1 = firstZ() - 1 - k;
       int k2 = lastZ() - 1 - k;
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
       k1 = lastZ() + k;
       k2 = firstZ() + k;
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
   }
}
void TGrid3D::fillbc_copy() 
{
   using namespace ndarray ;
   TNDArray4D Aux = getAux() ;
   for (int i=0 ; i < NbndryX()  ; i++) {
       int i1 = firstX() - 1 - i;
       int i2 = firstX();
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
       i1 = lastX() + i;
       i2 = lastX() - 1;
       Aux[view(i1)()()()] =  Aux[view(i2)()()()];
   }
   for (int j=0 ; j < NbndryY()  ; j++) {
       int j1 = firstY() - 1 - j;
       int j2 = firstY();
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
       j1 = lastY() + j;
       j2 = lastY() - 1;
       Aux[view()(j1)()()] =  Aux[view()(j2)()()];
   }
   for (int k=0 ; k < NbndryZ()  ; k++) {
       int k1 = firstZ() - 1 - k;
       int k2 = firstZ();
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
       k1 = lastZ() + k;
       k2 = lastZ() - 1;
       Aux[view()()(k1)()] =  Aux[view()()(k2)()];
   }
}

TNDArray3D TGrid3D::makeGrid() 
{
  return ndarray::allocate(ndarray::makeVector(NgridX(), NgridY(), NgridZ()))  ;
}
TNDArray4D TGrid3D::makeGrid(const int &ncomponents) 
{ 
  return ndarray::allocate(ndarray::makeVector(NgridX(), NgridY(), NgridZ(), ncomponents))  ;
}


