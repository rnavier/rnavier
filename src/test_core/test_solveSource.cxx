#include <iostream>
#include "TRNavier3DBj.h"
#include "TBRSSS.h"
#include "TEOS.h"
#include "TBRSSSSolverS2.h"
#include "TNDArray.h"

using namespace std ;
using namespace ndarray;

//!  Run hydro and do some analysis
int main(int argc, char **argv) 
{
    if  (argc < 2) {
        cout << "Usage is: test_solveSource.exe solveSource.ini" << endl;
        exit(1) ;
    }
    TRNavier3DBj rn(argv[1], make_eos_gammalaw) ;
    TBRSSS *model =  dynamic_cast<TBRSSS*>(rn.getHModel3D()) ;
    TBRSSSSolverS2 solver2(model);

    cout << " ---- Inversion test ----" << endl;  
    TNDArray1D qs, dq;
    TNDArray1D p0,p1,p2, p3;
    qs =ndarray_alloc(TBRSSS::Ncharge());
    dq =ndarray_alloc(TBRSSS::Ncharge());
    p0 =ndarray_alloc(TBRSSS::Naux());
    p1 =ndarray_alloc(TBRSSS::Naux());
    p2 =ndarray_alloc(TBRSSS::Naux());
    p3 =ndarray_alloc(TBRSSS::Naux());
    p0.deep() = 0.;
    // Primitives
    double tau = 0.6;
    double e = 9.43;
    double ux = 0.11;
    double uy = 1.0;
    double uz =-0.5;
    double ut = sqrt(1.+ux*ux+uy*uy+uz*uz);
    double eta = 4.5;
    double theta = ut/tau;
    double pixx = -eta*2./3.*theta;
    double pixy = 0.5;
    double pixz = 0.3;
    double piyy = -eta*2./3.*theta;
    double piyz = -0.1;
    double n = 0.93;
    double u[3] = {ux, uy, uz} ;
    double pi[5] = {pixx, pixy, pixz, piyy, piyz} ;
    model->setaux(e, n, u, pi, p0) ;


    p1.deep() = p0;
    p1[0] +=0.1;
    p1[1] +=0.2;
    p1[2] +=0.2;
    p1[3] +=0.2;
    p1[4] +=0.1;
    p1[5] +=0.1;
    p1[6] +=0.1;
    p1[7] +=0.1;
    p1[8] +=0.1;
    p1[9] +=0.1;
    model->fillaux(p1) ;

    cout << "Input primitives" << endl;
    cout << p1 << endl;
    //p2.deep()=p1;
    p3.deep()=p1;
    model->charges(tau, p0, qs);
    cout << "Charges" << endl;
    cout << qs << endl;
    int iter2 = solver2.invert(tau, qs, p3);
    cout << "Solver 2 did " << iter2 << " iterations"  << endl;

    printf(" Input\tSolver2\tAnswer\n");

    for (int m =0; m < TBRSSS::Ncharge(); m++){
        printf("%6.3f\t%6.3f\t%6.3f\n", p1[m],p3[m], p0[m]);
    }
} 







