#ifndef RN_TDomain3D_h
#define RN_TDomain3D_h

#include "gsl/gsl_math.h"
#include "TNDArray.h"
#include<iostream>

//! \brief Describes the lattice in three dimensions.
//!
//!\par Implementation Notes:
//!
//! - This cls uses the default constructor, and in general the default
//!   constructor needs to work,
//!
//! - Futher the asignment operator \c=. This operator needs to work
class TDomain3D {

   private:

      static const int fNbndryX=2 ; //!< The number of boundary points 2
      static const int fNbndryY=2 ; //!< The number of boundary points 2
      static const int fNbndryZ=2 ; //!< The number of boundary points 2

      int fNx ;        //!<  The number of X grid points
      int fNy ;        //!<  The number of Y grid points 
      int fNz ;        //!<  The number of Z grid points 

      double fDx     ; //!<  The grid spacing dx
      double fXmin   ; //!<  The minimum value of x
      double fDy     ; //!<  The grid spacing dy 
      double fYmin   ; //!<  The minimum value of y
      double fDz     ; //!<  The grid spacing dy 
      double fZmin   ; //!<  The minimum value of y

      double fTime   ;    //!<  Gives the curren time of the grid

   public:

      //! TDomain3D(const TDomain3D &) using default
      //! TDomain3D& operator= (const TDomain3D &) using default
      
      //! Read in the domain parameters. 
      void   read(std::istream &in)  ;
      void   write(std::ostream &ostream) const ;
      
      //! set the Current time
      void setDomainTime(double &t)  { fTime =t; }
      //! defines an equally spaced grid with cell edges between xmin and xmax
      void setXgrid(const double &xmin, const double &xmax, const int &nx) 
      {
         fXmin = xmin; fNx = nx; fDx = (xmax - xmin)/fNx;
      }
      //! defines an equally spaced grid with cell edges between ymin and ymax.
      void setYgrid(const double &ymin, const double &ymax, const int &ny) 
      {
         fYmin = ymin; fNy = ny; fDy = (ymax - ymin)/fNy;
      }
      //! defines an equally spaced grid with cell edges between zmin and zmax.
      void setZgrid(const double &zmin, const double &zmax, const int &nz) 
      {
         fZmin = zmin; fNz = nz; fDz = (zmax - zmin)/fNz;
      }

      //! Returns the current times
      double getDomainTime()   const { return fTime; }

      //! Returns delta x
      double getDx() const {return fDx ;}
      //! Returns delta y
      double getDy() const {return fDy ; }
      //! Returns delta z
      double getDz() const {return fDz ; }

      //! Returns the value of x at the center of lattice cite i 
      double getX(const int &i) const 
         { return fXmin + (i-fNbndryX) * fDx + fDx*0.5 ; } 
      //! Returns the value of y at the center of lattice cite j 
      double getY(const int &j) const 
         { return fYmin + (j-fNbndryY) * fDy + fDy*0.5 ; } 
      //! Returns the value of z at the center of lattice cite k 
      double getZ(const int &k) const 
         { return fZmin + (k-fNbndryZ) * fDz + fDz*0.5 ; } 

      //! Returns x y, and z at the center of lattice cite i & j & k
      void getXYZ(const int &i, const int &j, const int &k, double &x, double &y, double &z) 
         { x = fXmin + (i - fNbndryX)*fDx +  0.5*fDx ; 
           y = fYmin + (j - fNbndryY)*fDy +  0.5*fDy ; 
           z = fZmin + (k - fNbndryZ)*fDz +  0.5*fDz ; 
           return ;
         }


      //! The number of ghost cells in the X direction
      int NbndryX()  const {return fNbndryX ;}
      //! The number of ghost cells in the Y direction
      int NbndryY()  const {return fNbndryY ;}
      //! The number of ghost cells in the Z direction
      int NbndryZ()  const {return fNbndryZ ;}
      
      //! The number of cells in the x direction not including the boundary cells
      int Nx()      const {return fNx ;}
      //! The number of cells in the y direction not including the boundary cells
      int Ny()      const {return fNy ;}
      //! The number of cells in the z direction not including the boundary cells
      int Nz()      const {return fNz ;}
      
      //! The total number of lattice sites in the X-direction
      int NgridX()   const { return 2*fNbndryX+fNx; }
      //! The total number of lattice sites in the Y-direction
      int NgridY()   const { return 2*fNbndryY+fNy; }
      //! The total number of lattice sites in the Z-direction
      int NgridZ()   const { return 2*fNbndryZ+fNz; }

      //! The index of the first lattice site
      int firstX()   const { return fNbndryX ;}
      //! The index of the first lattice site
      int firstY()   const { return fNbndryY ;}
      //! The index of the first lattice site
      int firstZ()   const { return fNbndryZ ;}

      //! \brief The index of the last lattice site
      //! in the computational domain. 
      int lastX()    const { return fNbndryX + fNx  ; }
      //! \brief The index of the last lattice site
      //! in the computational domain. 
      int lastY()    const { return fNbndryY + fNy  ; }
      //! \brief The index of the last lattice site
      //! in the computational domain. 
      int lastZ()    const { return fNbndryZ + fNz  ; }
      

      //! \brief Returns the nearest gridpoint (i,j,k) to physical point (x,y,z) 
      void getIJK(double x, double y, double z, int &i, int &j, int &k)  ;
      //! \brief Returns the nearest gridpoint (i) to physical point (x) 
      int getI(double x) ;
      //! \brief Returns the nearest gridpoint (j) to physical point (y) 
      int getJ(double y) ;
      //! \brief Returns the nearest gridpoint (k) to physical point (z) 
      int getK(double z) ;

      //! Returns the min coordinate of the domain
      double Xmin() { return getX(firstX()) - 0.5*fDx ; }
      double Ymin() { return getY(firstY()) - 0.5*fDy ; }
      double Zmin() { return getZ(firstZ()) - 0.5*fDz ; }
      //! Returns the max coordinate of the domain
      double Xmax() { return getX(lastX()) - 0.5*fDx ; }
      double Ymax() { return getY(lastY()) - 0.5*fDy ; }
      double Zmax() { return getZ(lastZ()) - 0.5*fDz ; }
      
      //! Returns the center of the first cell on the lattice
      double Xfirst() { return getX(firstX()); }
      double Yfirst() { return getY(firstY()); }
      double Zfirst() { return getZ(firstZ()); }
      //! Returns the center of the last cell on the lattice
      double Xlast() { return getX(lastX()-1); }
      double Ylast() { return getY(lastY()-1); }
      double Zlast() { return getZ(lastZ()-1); }

} ;

#endif
