#include <memory>
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TDomain3D.h"
#include "TGrid3D.h"
#include "THydro3DBj.h"
#include "TICGbTest.h"
#include "TBRSSS.h"
#include "TNDArray.h"
#include <fenv.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include "TIOBuffer.h"
#include "TFOMaker.h"

using namespace std ;
//! This program evolves gubser flow for given time and maps the constant
//! temperature freezout surface.
//!
//! Usage:
//!  ./test_fo.exe nametag time T Twidth
//!  ./test_fo.exe gubserfo.ini 3.0 1.0 0.01   
int main(int argc, char **argv) 
{
    if  (argc < 5) {
        cout << "Usage is: test_fo.exe nametag time T Twidth" << endl;
        exit(1) ;
    }

    TRNavier3DBj rn(argv[1], make_eos_gammalaw) ;
    THydro3DBj hy(&rn, make_ic_icgbtest) ;
    hy.nextEvent() ;
    double T ;
    double width ;
    sscanf(argv[3],"%lf",&T);
    sscanf(argv[4],"%lf",&width);
    TGrid3D *grid = hy.getGrid();
    unique_ptr<TFOMaker> maker = unique_ptr<TFOMaker>(new TFOMaker(&rn, grid, T, width)) ;

    //Set up for FO surface
    double EntropyFlux = 0;
    double InitialEntropy = 0;
    double FinalEntropy = 0;
    TNDArray1D v = ndarray_alloc(4);
    // Loop over time steps
    int istep = 0  ;
    double tmax;
    sscanf(argv[2],"%lf",&tmax);

    const int nsteps_max = 20000 ;
//Take first step;
    {
        cout << "Time = " << hy.getTime() << endl;
        cout << "Time Step= " << hy.getDt() << endl;
        hy.step() ;

        double taun = grid->getTimeSlice(1);
        double taun1 = grid->getTimeSlice(0);
        double tau = 0.5*(taun+taun1);
        double dt = taun1-taun;
        double dx = grid->getDx();
        double dy = grid->getDy();
        double dz = grid->getDz()*TBRSSS::g33(tau);
        auto pn = grid->getAuxSlice(0);
        double u[4];
        double dEntropyFlux = 0;
        for (int i=grid->firstX(); i<grid->lastX(); i++) {
        for (int j=grid->firstY(); j<grid->lastY(); j++) {
        for (int k=grid->firstZ(); k<grid->lastZ(); k++) {

        auto pn1 = pn[i][j][k];
            if (maker->isInside(i,j,k)) {
                maker->make(i,j,k,v);
                double s = hy.getHModel3D()->getS(pn1);
                double ut = hy.getHModel3D()->getU0(pn1);
                InitialEntropy += s*ut*dx*dy*dz;
            }
            if (maker->nearSurface(i,j,k)) {
                maker->make(i,j,k, v);
                double s = hy.getHModel3D()->getS(pn1);
                hy.getHModel3D()->getU(pn1,u);
                dEntropyFlux += s*(v[0]*u[0]
                        +v[1]*u[1]
                        +v[2]*u[2]
                        +v[3]*u[3]
                        )*dx*dy*dz*dt;

            }
        }
        }
        }

        EntropyFlux += dEntropyFlux;
    }
    std::ofstream Out;
    string str = argv[1];
    str.erase(str.end()-4, str.end());
    Out.open(str+".out") ;
    TIOBuffer b;
    while(hy.getTime() < tmax && (istep < nsteps_max))  {
        cout << "Time = " << hy.getTime() << endl;
        cout << "Time Step= " << hy.getDt() << endl;
        hy.step() ;

        double taun = grid->getTimeSlice(1);
        double taun1 = grid->getTimeSlice(0);
        double tau = 0.5*(taun+taun1);
        double dt = taun1-taun;
        double dx = grid->getDx();
        double dy = grid->getDy();
        double dz = grid->getDz()*TBRSSS::g33(tau);
        auto pn = grid->getAuxSlice(0);

        int totcount = 0;
        int count = 0;
        double u[4];
        double rave = 0;
        double dEntropyFlux = 0;

        for (int i=grid->firstX(); i<grid->lastX(); i++) {
        for (int j=grid->firstY(); j<grid->lastY(); j++) {
        for (int k=grid->firstZ(); k<grid->lastZ(); k++) {

        auto pn1 = pn[i][j][k];
            totcount++;
                double x = grid->getX(i);
                double y = grid->getY(j);
                double z = grid->getZ(k)*TBRSSS::g33(tau);
                double tau = 0.5*(grid->getTimeSlice(0) + grid->getTimeSlice(1));
            if (maker->nearSurface(i,j,k)) {
                count++;
                maker->make(i,j,k,v);
                double s = hy.getHModel3D()->getS(pn1);
                double T = hy.getHModel3D()->getT(pn1);
                rave +=sqrt(x*x+y*y);
                //cout << sqrt(x*x+y*y) << " " << T  << " " << maker->getThetaAverage(i,j,k)<<endl;
                hy.getHModel3D()->getU(pn1,u);
                dEntropyFlux += s*(v[0]*u[0]
                        +v[1]*u[1]
                        +v[2]*u[2]
                        +v[3]*u[3]
                        )*dx*dy*dz*dt;
    Out << b.format(
            "%15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n",
            tau, x, y, z, maker->getTheta(i,j,k),
            v[0], v[1], v[2], v[3]);
            }


        }
        }
        }
        //if  (dEntropyFlux ==0.) {
        //   break ;
        //} 
        Out << "\n\n";
    
        EntropyFlux += dEntropyFlux;

        cout << "Average transverse radius = " <<rave/ ((double) count +1e-6) << endl;
        cout << "EntropyFlux: " <<dEntropyFlux<< endl;
        printf("Counts: %d \t out of %d\n", count, totcount);

        istep++ ;
    }

    Out.close();

 {
        cout << "Time = " << hy.getTime() << endl;
        cout << "Time Step= " << hy.getDt() << endl;
        hy.step() ;

        double taun = grid->getTimeSlice(1);
        double taun1 = grid->getTimeSlice(0);
        double tau = 0.5*(taun+taun1);
        double dt = taun1-taun;
        double dx = grid->getDx();
        double dy = grid->getDy();
        double dz = grid->getDz()*TBRSSS::g33(tau);
        auto pn = grid->getAuxSlice(0);

        double u[4];
        double dEntropyFlux = 0;
        for (int i=grid->firstX(); i<grid->lastX(); i++) {
        for (int j=grid->firstY(); j<grid->lastY(); j++) {
        for (int k=grid->firstZ(); k<grid->lastZ(); k++) {
        auto pn1 = pn[i][j][k];
            if (maker->isInside(i,j,k)) {
                maker->make(i,j,k, v);
                double s = hy.getHModel3D()->getS(pn1);
                double ut = hy.getHModel3D()->getU0(pn1);
                FinalEntropy += s*ut*dx*dy*dz;
            }
            if (maker->nearSurface(i,j,k)) {
                maker->make(i,j,k,v);
                double s = hy.getHModel3D()->getS(pn1);
                hy.getHModel3D()->getU(pn1,u);
                dEntropyFlux += s*(v[0]*u[0]
                        +v[1]*u[1]
                        +v[2]*u[2]
                        +v[3]*u[3]
                        )*dx*dy*dz*dt;
            }
        }
        }
        }

        EntropyFlux += dEntropyFlux;
 }

    cout << "Total Entropy Flux" << endl;
    cout << EntropyFlux << endl;
    cout << "Initial Entropy" << endl;
    cout << InitialEntropy<< endl;
    cout << "Final Entropy" << endl;
    cout << FinalEntropy<< endl;
    cout << "Entropy (non)-conservation" << endl;
    cout << (InitialEntropy-FinalEntropy-EntropyFlux) << endl;
    cout << "Entropy (non)-conservation/Entropy Flux " << endl;
    cout << (InitialEntropy-FinalEntropy-EntropyFlux)/EntropyFlux << endl;
}
