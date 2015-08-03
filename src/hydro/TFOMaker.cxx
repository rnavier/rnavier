/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include "TFOMaker.h"
#include "TLimiter.h"
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"

#include "TGrid3D.h"
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
TFOMaker::TFOMaker(TRNavier3DBj *rn, TGrid3D *gr, const double s0, const double width):
    fGrid(gr), fS0(s0), fWidth(width)
{
    fModel =  dynamic_cast<TBRSSS*>(rn->getHModel3D());
    
    // Read inputs
    TIOBuffer b;
    std::ifstream &in = rn->getInputFile();
    b.read_section(in, "TFOMaker");
    b.getS("fThetaName", fThetaName) ;
    b.getD("fThetaCut", fThetaCut) ;

    // Write the inputs to log
    clog << "# Writing the FO maker" << endl;
    clog << "[TFOMaker]" <<endl ;
    b.writeLineS(clog, fThetaName, "fThetaName: Type of freezout surface") ;
    b.writeLineD(clog, fS0, "fS0: The freezout location") ;
    b.writeLineD(clog, fWidth, "fWidth: The width of freezout surface") ;
    clog << "\n" <<endl ;
};

double TFOMaker::getScl(const int &i,const int &j, const int &k, const int &tSlice)
{
    if (fThetaName == "temperature"){
        auto pn = fGrid->getAuxSlice(tSlice);
        return fModel->getT(pn[i][j][k]);
           }
    else if (fThetaName == "radius"){
        double x = fGrid->getX(i);
        double y = fGrid->getY(j);
        double z = fGrid->getZ(k)*TBRSSS::g33(fGrid->getTime());
        return sqrt(x*x+y*y+z*z);
    }
        else {

                clog << "Unrecognized theta function! \n" << fThetaName << endl ;
        abort();
    }
};

double TFOMaker::theta( const int &i, const int &j, const int  &k, const int &tSlice)
{
        double scl = getScl(i,j,k, tSlice);
        return 0.5*erfc(-(scl-fS0)/fWidth);
};

//! Take the average of theta values and check against the threshold fThetaCut
bool TFOMaker::nearSurface( const int &i, const int &j, const int &k)
{
    double theta_average = thetaAverage(i,j,k);
    return (theta_average  > fThetaCut and 
            theta_average < 1.0-fThetaCut) ? true : false; 
};
bool TFOMaker::isSurface( const int &i, const int &j, const int &k)
{
    double th = theta(i,j,k);
    return (th > fThetaCut and 
            th < 1.0-fThetaCut) ? true : false; 
};
double TFOMaker::thetaAverage( const int &i, const int &j, const int &k) 
{
    double theta_sum = 0;
    int ntotal = 0;
    for (int di = -1; di <= 1; di++) {
    for (int dj = -1; dj <= 1; dj++) {
    for (int dk = -1; dk <= 1; dk++) {
        theta_sum += theta(i+di, j+dj, k+dk, 0);
        ntotal++;
        theta_sum += theta(i+di, j+dj, k+dk, 1);
        ntotal++;
    }
    }
    }
    theta_sum /= ((double) ntotal);
    return  theta_sum;
};

void TFOMaker::make( const int &i, const int &j, const int &k, const TNDArray1D &dSmu)
{
    using ndarray::view;
    double dt = fGrid->getTimeSlice(0) - fGrid->getTimeSlice(1);
    double tau = 0.5*(fGrid->getTimeSlice(0) + fGrid->getTimeSlice(1));

    double dx = fGrid->getDx();
    double dy = fGrid->getDy();
    double dz = fGrid->getDz()*TBRSSS::g33(tau);
    dSmu[0] = -(theta(i,j,k,0)-theta(i,j,k,1))/dt;
    dSmu[1] = -(theta(i+1,j,k)- theta(i-1,j,k))/(2*dx);
    dSmu[2] = -(theta(i,j+1,k)- theta(i,j-1,k))/(2*dy);
    dSmu[3] = -(theta(i,j,k+1)- theta(i,j,k-1))/(2*dz);

}; 

