/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#ifndef RN_TGrid3D_h
#define RN_TGrid3D_h

#include <vector>
#include "TNDArray.h"
#include "THModel3D.h" 
#include "TDomain3D.h"

class TRNavier3DBj ;

//#include "TBC3D.h"

//! \brief This class holds the grid in three dimensions. 
//!
//! It must initialized by passing a THModel3D object(See the THModel3D
//! documenation for a description of this). The THModel3D object
//! contains information on the number of conserved fields and the
//! number of auxiliary fields.  
//!
//! For instance, \c fAux is  a vector of 4D arrays. fAux[1][i,j,k, n] is the 
//! value of the n-th primitive field at x-position  x=gr->GetX(i),
//! y-position y =gr->GetY(j) and z-position z = gr->GetZ(k) and time slice
//! time = gr.getTimeSlice(1);
//!
//! The number of stored time slices is \c fNhistory=2, the most
//! recent is at  ffAux[0]. 
//! Generally  a stepping algorithm will take the data in fAux[1]  
//! and place the updated data in fAux[0].  Thus in order to have acess to
//! more than one time, one should typically call 
//!
//! \code
//! //TGrid *gr ;
//! gr->store() 
//! \endcode
//!
//! before the stepper. 
//! 
class TGrid3D: public TDomain3D {

   private: 
      
      THModel3D   *fModel ;

      //! The size of the history of past grids stored
      static const int fNhistory=2 ;
      //! \brief A vector containing past history of the grids
      //! fAuxs[0] is the current time slice, fAuxs[1] is the previous slice
      //! etc
      std::vector<TNDArray4D> fAuxs ;
      //! \brief A vector containing the times of the previous stored slices
      //! fTimes[0]  is the same as getTime() ;
      std::vector<double> fTimes ;

      double fTStart ;

      //! Read in the domain parameters. 
      void   read(std::istream &in)  ;
      void   write(std::ostream &ostream) const ;

      std::string fBCName ;
      void   fillbc_periodic() ;
      void   fillbc_copy() ;
      void   fillbc_periodicz() ;

   public:

      TGrid3D(TRNavier3DBj *rn) ; 
      THModel3D *getModel() { return fModel; }
      int Nhistory() const {return fNhistory ;} 

      //! Returns the starting time
      double getTStart() const { return fTStart; }

      //! Sets the starting time. 
      void   setTStart(const double &to) { fTStart = to ; }
      
      //! increments the time
      void updateTime(const double &dt) 
      { fTimes[0]+=dt ; setDomainTime(fTimes[0]); }

      //! Sets the time
      void setTime(double &t) { fTimes[0]= t; setDomainTime(t); } 
      
      //! Get the most recent data slice
      double getTime() { return fTimes[0] ; }

      //! Restarts the clock at TStart and nothing else
      void resetTime() { fTimes[0] = fTStart ; setDomainTime(fTimes[0]); }

      //! Restarts the clock and TStart at to
      void resetTime(const double &to) { setTStart(to) ; resetTime() ; }
      
      //! Get the most recent data slice
      TNDArray4D &getAux() { return fAuxs[0] ; }
     
      //! Get a particular data slice (islice=0 is the most recent)
      TNDArray4D &getAuxSlice(const int &islice) { return  fAuxs[islice] ; }
      //! Get time at a particular data slice (islice=0 is the most recent)
      double getTimeSlice(const int &islice) { return  fTimes[islice] ; }

      //! Returns vector of all past arrays
      std::vector<TNDArray4D> &getAuxs() { return fAuxs; }
      //! Returns vector of all past time slices
      std::vector<double> &getTimes() { return fTimes; }

      //! Pushes data arrays by one past slice poistion.
      void store() ;
      
      //! fills up the boundary conditions of Auxs[0]
      void fillbc() ;

      //! Allocates a 3D grid for a single scalar field
      TNDArray3D makeGrid() ;
      //! Allocates a 4D grid for a n-component field
      TNDArray4D makeGrid(const int &ncomponents) ;
};



#endif
