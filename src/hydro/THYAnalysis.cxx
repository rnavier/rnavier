/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
#include "THydro3DBj.h"
#include "TBRSSS.h"
#include "THYAnalysis.h"
#include "TRNavier3DBj.h"
#include "TIOBuffer.h"
#include "TGrid3D.h"
#include "THModel3D.h"
#include "TNDArray.h"
#include "numeric.h"
#include "TEOS.h"

//! Write the run tag to the output stream.
void write_run_tag(std::ostream &out, const int &run, const double &time);  

//! Write a specified point to the stream
void write_point(std::ostream &out, 
      const double &t,
      const double &x,
      const double &y,
      const double &z,
      THModel3D *hm, 
      const TNDArray1D &auxi) ;

void write_point(std::ostream &out,
      THydro3DBj *hy,
      const int &i,
      const int &j,
      const int &k, bool plot_solution) ;

////////////////////////////////////////////////////////////////////////////////

THYAnalysis_2dslice::THYAnalysis_2dslice(THydro3DBj *hy, const double &time_out, const string &plane, const double &slice, const bool &plot_solution):
    fTOut(time_out), fPlane(plane), fHydro(hy), fPlotSolution(plot_solution)  
{
    fTNextPrint = fHydro->getTime() ; 
    if (fPlane=="xy") {
       fLocation=hy->getGrid()->getK(slice) ;
    } else if (fPlane=="xz") {
       fLocation=hy->getGrid()->getJ(slice) ;
    } else if (fPlane=="yz") {
       fLocation=hy->getGrid()->getI(slice) ;
    } else {
       cerr << "** THYAnalysis_2dslice::THYAnalysis_2dslice ** Bad selection in switch. Selected plane " << plane << ". Choices are xy, xz, or yz. " << endl;
    }
    // Open filename nametag nametag_2dslice_xz_0035.out
    TIOBuffer b ;
    string s = hy->getRNavier()->getNametag() + "_2dslice" + "_" + fPlane + b.format("_%04d.out",static_cast<int>(10*slice)) ;
    fOut.open(s.c_str()) ;
} 
void THYAnalysis_2dslice::start() 
{
    write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
    write_2dslice() ;
    fTNextPrint = fHydro->getTime() + fTOut ;
}
void THYAnalysis_2dslice::analyze() 
{
    double time = fHydro->getTime() ;
    if (time > fTNextPrint) {
        write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
        write_2dslice() ;
        fTNextPrint += fTOut ; 
    }
}
void THYAnalysis_2dslice::stop() 
{
    fOut << "\n\n" ;
}

void THYAnalysis_2dslice::write_2dslice() 
{
   if (fPlane == "yz")  {
      int i = fLocation ;
      for (int j=0 ; j < fHydro->getGrid()->NgridY() ; j++) {
      for (int k=0 ; k < fHydro->getGrid()->NgridZ() ; k++) {
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut << '\n' ;
      }
      fOut << '\n' ; 
      }
   } else if (fPlane == "xz") {
      int j = fLocation ;
      for (int i=0 ; i < fHydro->getGrid()->NgridX() ; i++) {
      for (int k=0 ; k < fHydro->getGrid()->NgridZ() ; k++) {
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut << '\n' ; 
      }
      fOut << '\n' ; 
      }
   } else if (fPlane == "xy") {
      int k = fLocation ;
      for (int i=0 ; i < fHydro->getGrid()->NgridX() ; i++) {
      for (int j=0 ; j < fHydro->getGrid()->NgridY() ; j++) {
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut << '\n' ;
      }
      fOut << '\n' ; 
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

THYAnalysis_1dslice::THYAnalysis_1dslice(THydro3DBj *hy, double time_out, const string &line, const double &yslice, const double &zslice, const bool &plot_solution):
    fTOut(time_out), fLine(line), fHydro(hy), fPlotSolution(plot_solution) 
{
    fTNextPrint = fHydro->getTime() ; 
    if (fLine=="x") {
       fLocation2=hy->getGrid()->getJ(yslice) ;
       fLocation3=hy->getGrid()->getK(zslice) ;
    } else if (fLine=="y") {
       fLocation2=hy->getGrid()->getI(yslice) ;
       fLocation3=hy->getGrid()->getK(zslice) ;
    } else if (fLine=="z") {
       fLocation2=hy->getGrid()->getI(yslice) ;
       fLocation3=hy->getGrid()->getJ(zslice) ;
    } else {
       cerr << "** THYAnalysis_2dslice::THYAnalysis_2dslice ** Bad selection in switch. Selected line " << line << ". Choices are x, y, or z. " << endl;
    }

    // Open filename nametag nametag_2dslice_xz_0035.out
    TIOBuffer b ;
    string s = hy->getRNavier()->getNametag() + "_1dslice" + "_" + fLine + b.format("_%04d_%04d.out",(int)(10*yslice),(int)(10*zslice)) ;
    fOut.open(s.c_str()) ;
} 

//! Write the run tag and the line of interest, Update the time for next
//! printing.
void THYAnalysis_1dslice::start() 
{
    write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
    write_1dslice() ;
    fTNextPrint = fHydro->getTime() + fTOut ;
}

//! If it is time for printing, write the run tag, the line of interest, and update the time for next printing.
void THYAnalysis_1dslice::analyze() 
{
    double time = fHydro->getTime() ;
    if (time > fTNextPrint) {
        write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
        write_1dslice() ;
        fTNextPrint += fTOut ; 
    }
}
void THYAnalysis_1dslice::stop() 
{
    fOut << "\n\n" ;
}

void THYAnalysis_1dslice::write_1dslice() 
{
   if (fLine == "x")  {
      for (int i=0 ; i < fHydro->getGrid()->NgridX() ; i++) {
         int j = fLocation2 ;
         int k = fLocation3 ;
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut<< '\n' ;
      }
   } else if (fLine == "y") {
      for (int j=0 ; j < fHydro->getGrid()->NgridY() ; j++) {
         int i = fLocation2 ;
         int k = fLocation3 ;
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut<< '\n' ;
      }
   } else if (fLine == "z") {
      for (int k=0 ; k < fHydro->getGrid()->NgridZ() ; k++) {
         int i = fLocation2 ;
         int j = fLocation3 ;
         write_point(fOut, fHydro, i, j, k, fPlotSolution) ;
         fOut<< '\n' ;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
THYAnalysis_point::THYAnalysis_point(THydro3DBj *hy, 
            const double &x, 
            const double &y, 
            const double &z, 
            const bool &plot_solution) : 
    fHydro(hy), fPlotSolution(plot_solution)  
{
    fLocation1=hy->getGrid()->getI(x) ;
    fLocation2=hy->getGrid()->getJ(y) ;
    fLocation3=hy->getGrid()->getK(z) ;

    TIOBuffer b ;
    string s = hy->getRNavier()->getNametag() + b.format("_point_%04d_%04d_%04d.out",
          (int)(10*x),
          (int)(10*y),
          (int)(10*z)) ;
    fOut.open(s.c_str()) ;
} 
void THYAnalysis_point::start() 
{
    write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
    write_point(fOut, fHydro, fLocation1, fLocation2, fLocation3, fPlotSolution) ;
    fOut << '\n' ;
}
void THYAnalysis_point::analyze() 
{
    write_point(fOut, fHydro, fLocation1, fLocation2, fLocation3, fPlotSolution) ;
    fOut << '\n' ;
}
void THYAnalysis_point::stop() 
{
    fOut << "\n\n" ;
}

////////////////////////////////////////////////////////////////////////////////


//! Write the gnuplot ready header at the start of a run.
void write_run_tag(ostream &out, const int &run, const double &time)  
{
    out << endl;
    out << endl;
    out << "# Run  " << run << endl;
    out << "#"<<endl;
    out << "# Time = " << time << endl;
}

void write_point(std::ostream &out,
      THydro3DBj *hy,
      const int &i,
      const int &j,
      const int &k, bool plot_solution)
{
    double t = hy->getTime() ;
    double x, y, z ;
    hy->getGrid()->getXYZ(i, j, k, x, y, z) ;
    TNDArray1D aux = hy->getGrid()->getAux()[ndarray::view(i)(j)(k)()]; 
    THModel3D *hm  = hy->getHModel3D() ;
    write_point(out, t, x, y, z, hm, aux) ;
    if (plot_solution) {
       double ax[THModel3D::Naux()] ;
       auto auxsol = THModel3D::view(ax) ;
       TICPoint3D pnt = {t, x, y, z};
       hy->getIC()->IC(pnt, auxsol) ;
       write_point(out, t, x, y, z, hm, auxsol) ;
    }
}

//! Write a specified point to the stream
void write_point(std::ostream &out, 
      const double &t,
      const double &x,
      const double &y,
      const double &z,
      THModel3D *hm, 
      const TNDArray1D &aux) 
{
    TIOBuffer b;
    double e = hm->getE(aux) ;
    double n = hm->getN(aux) ;
    double p = hm->getP(aux) ;
    TEOS *eos = hm->getEOS() ;
    double s, temper, mu ;
    eos->stmu(e, n, s, temper, mu) ;

    double u0 = hm->getU0(aux) ;
    double u1 = hm->getU1(aux) ;
    double u2 = hm->getU2(aux) ;
    double u3 = hm->getU3(aux) ;

    double pi11 = hm->getPi11(aux) ;
    double pi12 = hm->getPi12(aux) ;
    double pi13 = hm->getPi13(aux) ;
    double pi22 = hm->getPi22(aux) ;
    double pi23 = hm->getPi23(aux) ;
    double pi33 = hm->getPi33(aux) ; 

    out << b.format(
            "%15.5e %15.5e %15.5e %15.5e "
            "%15.5e %15.5e %15.5e %15.5e "
            "%15.5e %15.5e "
            "%15.5e %15.5e %15.5e %15.5e "
            "%15.5e %15.5e %15.5e %15.5e "
            "%15.5e %15.5e ",
            t, x, y, z,
            e, n, p, s, 
            temper, mu,
            u0, u1, u2, u3, 
            pi11, pi12, pi13, pi22, 
            pi23, pi33) ;
}

////////////////////////////////////////////////////////////////////////////////
THYAnalysis_entropy_vs_time::THYAnalysis_entropy_vs_time(THydro3DBj *hy, const double &z) : fHydro(hy) 
{
    fLocation3=hy->getGrid()->getK(z) ;
    TIOBuffer b ;
    string s = hy->getRNavier()->getNametag() + b.format("_entropy_vs_time_z%04d.out",
          (int)(10*z)) ;
    fOut.open(s.c_str()) ;
} 
void THYAnalysis_entropy_vs_time::start() 
{
    write_run_tag(fOut, fHydro->getRun(), fHydro->getTime()) ;
    write_entropy_vs_time(fHydro->getTime()) ;
    fOut << '\n' ;
}
void THYAnalysis_entropy_vs_time::analyze() 
{
    write_entropy_vs_time(fHydro->getTime()) ;
    fOut << '\n' ;
}
void THYAnalysis_entropy_vs_time::stop() 
{
    fOut << "\n\n" ;
}
void THYAnalysis_entropy_vs_time::write_entropy_vs_time(const double &tau) 
{
    TGrid3D *gr = fHydro->getGrid() ;
    THModel3D *hm = fHydro->getHModel3D() ;
    TNDArray4D aux = gr->getAux() ;

    double stot = 0. ;
    int k = fLocation3 ;
    for (int i=gr->firstX(); i < gr->lastX(); i++) {
    for (int j=gr->firstY(); j < gr->lastY(); j++) {
        double dx = gr->getDx() ;
        double dy = gr->getDy() ;
        double s = hm->getS(aux[i][j][k]) ;
        double u0 = hm->getU0(aux[i][j][k]) ;

        stot += dx*dy*TBRSSS::g33(tau)*s*u0 ;
    }
    }
    TIOBuffer b;
    fOut << b.format("%15.5e %15.5e ", tau, stot) ;
}


