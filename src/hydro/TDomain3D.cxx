#include "TDomain3D.h"
#include "TIOBuffer.h"
#include <iostream>

using namespace std ;


//! \brief Read in the 3D domain from an input file.
//!
//! The form of the input file is
//! \verbatim
//!#
//!# Set up the domain parameters
//!#
//![TDomain3D]
//!        60     = fNx: The number of grid points
//!    0.00000      = fXmin: The smallest value of x
//!    10.00000      = xmax: The largest value of x
//!        60      = fNy: The number of grid points
//!    0.00000      = fYmin: The smallest value of y
//!    10.0000      = ymax: The largest value of y
//!        60      = fNz: The number of grid points
//!    0.00000      = fZmin: The smallest value of z
//!    10.0000      = zmax: The largest value of z
//!    1.00000      = fTStart: The starting time
//![END]
//! 
//! \endverbatim
//!
void TDomain3D::read(istream &in) 
{
   TIOBuffer b;
   double xmin, xmax ;
   b.read_section(in, "TDomain3D") ;
   b.getI("fNx",fNx) ;
   b.getD("fXmin",xmin) ;
   b.getD("xmax",xmax) ;

   fXmin = xmin ;
   fDx   = (xmax - xmin)/(double)fNx ;

   double ymin, ymax ;
   b.getI("fNy",fNy) ;
   b.getD("fYmin",ymin) ;
   b.getD("ymax",ymax) ;

   fYmin = ymin ;
   fDy   = (ymax - ymin)/(double)fNy ;

   double zmin, zmax ;
   b.getI("fNz",fNz) ;
   b.getD("fZmin",zmin) ;
   b.getD("zmax",zmax) ;

   fZmin = zmin ;
   fDz   = (zmax - zmin)/(double)fNz ;

   b.getD("fTime",fTime) ;
}

//! Write the domain to an output file.
void TDomain3D::write(ostream &out) const
{
   TIOBuffer b;
   out << "[TDomain3D]" <<endl; 
   b.writeLineI(out, fNx,   "fNx: The number of grid points") ;
   b.writeLineD(out, fXmin, "fXmin: The smallest value of x") ;

   double xmax = fXmin + (double)fNx*fDx ;
   b.writeLineD(out, xmax,  "xmax: The largest value of x ") ;

   b.writeLineI(out, fNy,   "fNy: The number of grid points") ;
   b.writeLineD(out, fYmin, "fYmin: The smallest value of y") ;
   double ymax = fYmin + (double)fNy*fDy ;
   b.writeLineD(out, ymax,  "ymax: The largest value of y ") ;

   b.writeLineI(out, fNz,   "fNz: The number of grid points") ;
   b.writeLineD(out, fZmin, "fZmin: The smallest value of z") ;
   double zmax = fZmin + (double)fNz*fDz ;
   b.writeLineD(out, zmax,  "zmax: The largest value of z ") ;

   b.writeLineD(out, fTime,   "fTime: The starting time") ;
   out << "[END]" << '\n' <<endl; 
}

void TDomain3D::getIJK(double x, double y, double z, int &ix, int &jy, int &kz)
{
   ix = getI(x) ;
   jy = getJ(y) ;
   kz = getK(z) ;
}

int TDomain3D::getI(double x)
{
   int ixp = firstX() + std::round((x - Xfirst())/fDx) ;
   return GSL_MIN(GSL_MAX(ixp, 0), NgridX()-1) ;
}

int TDomain3D::getJ(double y)
{
   int iyp = firstY() + std::round((y- Yfirst())/fDy) ;
   return GSL_MIN(GSL_MAX(iyp, 0), NgridY()-1) ;
}
int TDomain3D::getK(double z)
{
   int izp = firstZ() + std::round((z- Zfirst())/fDz);
   return GSL_MIN(GSL_MAX(izp, 0), NgridZ()-1) ;
}

