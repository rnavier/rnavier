#ifndef RN_TFOMaker_h
#define RN_TFOMaker_h

#include<string>
#include "../hydro/TNDArray.h"
class TBRSSS ;
class TGrid3D;
class TRNavier3DBj ;

//! This class defines temperature/radius  surfaces
class TFOMaker{
    private:

        TGrid3D *fGrid ;
        TBRSSS *fModel ;

        std::string fThetaName; //! Freezout surface type
        double fThetaCut; //! Threshold for determening surface neighborhood
        double fS0; //! The position of the surface
        double fWidth; //! Width of the surface

     //! theta is a  scalar potential used to determinate the freezout surface.
     //! It is a smooth function interpolating between 0 and 1. Its argument
     //! is givent by getS(i,j,k) function. Position and width of the surface
     //! is determined by fS0 and fWidth
     double theta(const int &i, const int &j, const int &k, const int &tSlice=0);
        //! Takes the average of theta over neighbouring cells (including one
        //! past time slice)
        double thetaAverage( const int &i, const int &j, const int &k);
   public:

     ~TFOMaker() {; } 
     TFOMaker(TRNavier3DBj *rn, TGrid3D *gr, const double s0, const double width);
     //! Returns the threshold for surface determination
     double getThetaCut() const {return fThetaCut;} ;
     double getTheta(const int &i, const int &j, const int &k) {return theta(i,j,k);} ;
     double getThetaAverage(const int &i, const int &j, const int &k) {return thetaAverage(i,j,k);} ;
     //! Returns scalar argument used to determine freezout surface (e.g. temperature)
     virtual double getScl(const int &i, const int &j, const int &k, const int &tSlice=0);

     //! Returns true if a given cell is near the surface
     bool nearSurface(const int &i, const int &j, const int &k);
     bool isSurface(const int &i, const int &j, const int &k);
     //! Returns true if a given cell is inside the surface
     bool isInside(const int &i, const int &j, const int &k)
     { return (theta(i,j,k,0) >= 1.0-fThetaCut) ? true : false; }; 
     //! Returns true if a given cell is outside the surface
     bool isOutside(const int &i, const int &j, const int &k)
     { return (theta(i,j,k,0) <= fThetaCut) ? true : false; }; 
     //! \brief Returns \partial_\mu theta, where theta is the scalar potential determening
     //! the freezout surface. It also returns primitives at the center of the hypercube.
     void make(const int &i, const int &j, const int &k, const TNDArray1D &dSmu) ;
} ;


#endif //define RN_TFOMaker_h
