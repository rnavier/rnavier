#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include "gsl/gsl_math.h"
#include "numeric.h"
#include "TBRSSSStep.h"
#include "TBRSSSSolverS2.h"
#include "TIOBuffer.h"
#include "TGrid3D.h"
#include "TRNavier3DBj.h"
#include "TLimiter.h"

#include "THydro3DBj.h"

std::ofstream TBRSSSStep::fInversionLog;

TBRSSSStep::TBRSSSStep(TRNavier3DBj *rn, TGrid3D *gr) : fRNavier(rn), fGrid(gr)
{
    fModel =  dynamic_cast<TBRSSS*>(rn->getHModel3D());

    // Read inputs
    TIOBuffer b;
    std::ifstream &in = rn->getInputFile();
    b.read_section(in, "TBRSSSStep");
    b.getS("fStepType",  fStepType);
    b.getD("fCFL",  fCFL) ;
    b.getB("fIsIdealFluid", fIsIdealFluid) ;
    int icase ;
    b.getI("fAccelerationFlag", icase) ;
    fAccelerationFlag = static_cast<EAcceleration>(icase) ;
    b.getI("fPropagationFlag", icase) ;
    fPropagationFlag = static_cast<EPropagation>(icase) ;

    // Write the inputs to log
    clog << "# Writing the stepper" << endl;
    clog << "[TBRSSStep]" <<endl ;
    b.writeLineS(clog, fStepType, "fStepType : The name of time stepping algorithm") ;
    b.writeLineD(clog, fCFL, "fCFL: The default CFL condition used by stepper") ;
    b.writeLineB(clog, fIsIdealFluid, "fIsIdealFluid: The update step is assumed to be ideal") ;

    b.writeLineI(clog, static_cast<int>(fAccelerationFlag), "fAccelerationFlag: 0/1, The acceleration is computed using ideal hydro estimate/viscous hydro estimate") ;
    b.writeLineI(clog, static_cast<int>(fPropagationFlag), "fPropagationFlag: 0/1, The signal velocities are taken with the ideal charachteristics/speed of light") ;
    clog << "\n" <<endl ;

    //Allocate the solver
    fSolver = unique_ptr<TBRSSSSolver>(new TBRSSSSolverS2(fModel, fIsIdealFluid));

    // Allocate workspace for various stages
    if (fStepType == "Euler_IMEX") {
        fTau_rk.resize(2);
        fTau_rk[0] = fGrid->getTime();
        fTau_rk[1] = fGrid->getTime();
        fP_rk.resize(2);
        fP_rk[0] = fGrid->getAuxSlice(1);
        fP_rk[1] = fGrid->getAuxSlice(0);
    } else if (fStepType == "ARS222") {
        fTau_rk.resize(3);
        fTau_rk[0] = fGrid->getTime();
        fTau_rk[1] = fGrid->getTime();
        fTau_rk[2] = fGrid->getTime();
        fP_rk.resize(3);
        fP_rk[0] = fGrid->getAuxSlice(1);
        fP_rk[1] = fGrid->makeGrid(TBRSSS::Naux());
        fP_rk[2] = fGrid->getAuxSlice(0);
    } else {
        clog << "Unrecognized alogrithm! \n" << fStepType << endl ;
        abort() ;
    }
    // Allocate common workspace
    fQstar = fGrid->makeGrid(TBRSSS::Ncharge());
    fDaI   = fGrid->makeGrid(3);

    // Open the inversion log
    fBadStates = 0 ;
    fRealBadStates = 0 ;
    if (!fInversionLog.is_open()) {
        std::string s = fRNavier->getNametag() +  "_bad_points.out" ;
        fInversionLog.open(s.c_str()) ;
    }
}

double TBRSSSStep::getDt(TGrid3D &gr)
{
    if (&gr != fGrid) {
        cout << "** TBRSSSStep::getDt ** The grid passed as argument does not match the initialization grid!" << endl;
        return 0.;
    }
    double tau = fGrid->getTime();
    double dx  = fGrid->getDx();
    double dy  = fGrid->getDy();
    double dz  = TBRSSS::g33(tau)*fGrid->getDz();
    return fCFL*GSL_MIN(dx,GSL_MIN(dy, dz));
}

int TBRSSSStep::step(TGrid3D &gr, const double &dt)
{
    if (&gr != fGrid) {
        cout << "** TBRSSSStep::getDt ** The grid passed as argument does not match the initialization grid!" << endl;
        return -1;
    }


    if (fStepType == "Euler_IMEX") {
        eulerUpdateGrid(dt);
    } else if (fStepType == "ARS222") {
        ars222UpdateGrid(dt);
    } else {
        clog << "Unrecognized algorithm! \n" << fStepType << endl ;
        abort() ;
    }
    // Put a line break in the inversion log to separate steps from each
    // other.
    if (fBadStates > 0 or fRealBadStates> 0) {
        fInversionLog << "\n" << endl;
    }


    return 1 ;
}

void TBRSSSStep::accelerationCorrection(const int &kstage)
{
    using ndarray::view;
    TNDArray4D pn   = fP_rk[kstage];
    TNDArray4D pn1  = fP_rk[kstage+1];

    for (int i=fGrid->firstX(); i<fGrid->lastX(); i++) {
    for (int j=fGrid->firstY(); j<fGrid->lastY(); j++) {
    for (int k=fGrid->firstZ(); k<fGrid->lastZ(); k++) {
        double ut = fModel->getU0(pn[view(i)(j)(k)]);
        for (int I=0 ; I<3; I++) {
            double dtuI = Dtp(kstage,  pn1(i,j,k,1+I), i, j, k, 1+I);
            //Acceleration correction -ut dtuI_estimate was
            //already placed in fDaI. We hadd + ut dtuI_exact.
            //The acceleration correction  (delta a)_I, from
            // this time step.
            fDaI(i,j,k,I)  += ut*dtuI;
        }
    }
    }
    }
}

//! \brief Euler update q = q_0 - dt*(Flux(p^i)+ S_ex(p^i)) - dt*S_im(p^i))
//!
//! Butcher tableau for Euler IMEX
//! \verbatim
//!    _1_|_1_0   _1_|_0_1__
//!       | 1 0      | 0 1
//! \endverbatim
void TBRSSSStep::eulerUpdateGrid(const double &dt)
{
    double aex = 1.0;
    double aim = 1.0;
    TNDArray2D aex_rk = ndarray_view(&aex, 1, 1) ;
    TNDArray2D aim_rk= ndarray_view(&aim, 1, 1) ;

    fTau_rk[0] = fGrid->getTime();
    fTau_rk[1] = fGrid->getTime() + dt; //ex and im times must match

    TNDArray4D pn, pn1,qs;

    int kstage = 0;
    pn  = fP_rk[kstage];
    pn1 = fP_rk[kstage+1];
    qs  = fQstar;

    fillqs(kstage, dt, aex_rk, aim_rk, qs) ;
    stepStage(kstage, dt, fTau_rk[kstage], fTau_rk[kstage+1], pn, qs, aex, aim, pn1);
    accelerationCorrection(kstage);
    fGrid->fillbc();
}

//! \brief ARS222  diagonally implicit 2nd order rk update
//! U.M. Ascher et al. / Applied Numerical Mathematics 25 (1997) 151-167
//! See section 2.6
//!
//! update q = q_0 - dt*(Flux(p^i)+ S_ex(p^i)) - dt*S_im(p^i))
//!
//! Butcher tableau for ARS222
//! \verbatim
//! g = 1-1/sqrt(2)
//! d = 1-1/(2*g)
//!     g |  g  0    g | 0  g  0
//!    _1_|_1-d_d   _1_|_0_1-g_g_
//!       | 1-d d      | 0 1-g g
//! \endverbatim
//!
void TBRSSSStep::ars222UpdateGrid(const double &dt)
{
    const double g = 1.0-1.0/sqrt(2.0);
    const double d = 1.0-1.0/(2*g);
    double aex[2][2] ={{ g,   0 },
                       { d,  1-d }};
    double aim[2][2] ={{ g,   0 },
                       { 1-g, g}};
    TNDArray2D aex_rk = ndarray_view(&aex[0][0], 2, 2) ;
    TNDArray2D aim_rk = ndarray_view(&aim[0][0], 2, 2) ;

    fTau_rk[0] = fGrid->getTime();
    fTau_rk[1] = fGrid->getTime() + g*dt; //ex and im times must match
    fTau_rk[2] = fGrid->getTime() +   dt; //ex and im times must match

    fP_rk[1].deep() = fP_rk[0];

    TNDArray4D pn, pn1, qs;

    // Take the first step
    int ks = 0;
    pn  = fP_rk[ks];
    pn1 = fP_rk[ks+1];
    qs  = fQstar;

    fillqs(ks, dt, aex_rk, aim_rk, qs) ;
    stepStage(ks, dt, fTau_rk[ks], fTau_rk[ks+1], pn, qs, aex[ks][ks], aim[ks][ks], pn1);
    accelerationCorrection(ks);
    fillBoundary(pn1);

    // The second stage
    ks = 1;
    pn  = fP_rk[ks];
    pn1 = fP_rk[ks+1];

    fillqs(ks, dt, aex_rk, aim_rk, qs) ;

    stepStage(ks, dt, fTau_rk[ks], fTau_rk[ks+1], pn, qs, aex[ks][ks], aim[ks][ks], pn1);

    accelerationCorrection(ks);
    fillBoundary(pn1) ;
}

//! Fill up the boundary conditions for a specified stage.
//! Overwrite
void TBRSSSStep::fillBoundary(const TNDArray4D &pn1)
{
    fGrid->getAux().deep() = pn1 ;
    fGrid->fillbc() ;
    pn1.deep() = fGrid->getAux() ;
}

void TBRSSSStep::fillqs(const int kstage, 
        const double &dt,
        const TNDArray2D &aex_rk, 
        const TNDArray2D &aim_rk,
        const TNDArray4D &qs)
{
    using ndarray::view ;
    TNDArray2D qs_rk = ndarray_alloc(kstage+1, TBRSSS::Ncharge()) ;
    TNDArray1D pn = ndarray_alloc(TBRSSS::Naux()) ;
    TNDArray1D pn1 = ndarray_alloc(TBRSSS::Naux()) ;
    TNDArray1D qn1 = ndarray_alloc(TBRSSS::Ncharge()) ;
    TNDArray1D src_im = ndarray_alloc(TBRSSS::Ncharge()) ;
    TNDArray1D dq_ex = ndarray_alloc(TBRSSS::Ncharge()) ;

    for (int i=fGrid->firstX(); i<fGrid->lastX(); i++) {
    for (int j=fGrid->firstY(); j<fGrid->lastY(); j++) {
    for (int k=fGrid->firstZ(); k<fGrid->lastZ(); k++) {

        // Extract the original charges and copy it to all stages
        pn.deep() = fP_rk[0][view(i)(j)(k)()] ;
        double taun = fTau_rk[0] ;
        fModel->charges(taun, pn, qs_rk[ view(0)() ]) ;

        for (int istg = 1 ; istg <= kstage; istg++) {
            qs_rk[ view(istg)() ] = qs_rk[ view(0)() ] ;
        } 

        // Loop over stages up to kstage and fill up qs_rk.
        // istg is an index for the contibution of istg-1, 
        // to the q* of  all higher stages up to nstag
        for (int istg = 1 ; istg <= kstage ; istg++) {
            // extract the src_im for the previous stage (istg-1)
            double taun1 = fTau_rk[istg] ;
            pn1.deep() = fP_rk[istg][ view(i)(j)(k)() ] ;
            fSolver->source(taun1, pn1, src_im) ;
           
            // Determine the dq_ex for the previous stage (istg-1)
            // Using qn1 = qst + aex*dq_ex + aim *dt * src
            fModel->charges(taun1, pn1, qn1) ;

            double aex = aex_rk(istg-1, istg-1) ;
            double aim = aim_rk(istg-1, istg-1) ;

            dq_ex.deep() = (qn1 - qs_rk[view(istg-1)()] - aim* dt*src_im)/aex ;

            //Fill up the qs_rk  for all higher stages with the dq_ex, and src
            // from the (i-1)-th stage.
            for (int jstg=istg; jstg <= kstage; jstg++) {
                aex = aex_rk(jstg, istg-1)  ;
                aim = aim_rk(jstg, istg-1)  ; 
                qs_rk[ view(jstg)() ] += aex * dq_ex + aim * dt * src_im ;
            }
        }
        qs[view(i)(j)(k)()] = qs_rk[view(kstage)()] ;
    }
    }
    }
}
//! loops over the grid of pn and updates pn1
void TBRSSSStep::stepStage(const int &kstage,
        const double &dt,
        const double &taun0,
        const double &taun1,
        const TNDArray4D &pn,
        const TNDArray4D &qs,
        const double &aex,
        const double &aim,
        const TNDArray4D &pn1)
{
    using ndarray::view;
    //used for x,y,z loop
    //CFL condition variables
    char direction;
    double dxI[3];
    dxI[0] = fGrid->getDx() ;
    dxI[1] = fGrid->getDy() ;
    dxI[2] = fGrid->getDz()*TBRSSS::g33(taun0) ;
    double lamx = dt/dxI[0] ;
    double lamy = dt/dxI[1] ;
    double lamz = dt/dxI[2] ;

    TNDArray2D pnXslice =  ndarray_alloc(fGrid->NgridX(), TBRSSS::Naux());
    TNDArray2D pnYslice =  ndarray_alloc(fGrid->NgridY(), TBRSSS::Naux());
    TNDArray2D pnZslice =  ndarray_alloc(fGrid->NgridZ(), TBRSSS::Naux());

    TNDArray2D dqnXslice = ndarray_alloc(fGrid->NgridX(), TBRSSS::Ncharge());
    TNDArray2D dqnYslice = ndarray_alloc(fGrid->NgridY(), TBRSSS::Ncharge());
    TNDArray2D dqnZslice = ndarray_alloc(fGrid->NgridZ(), TBRSSS::Ncharge());

    TNDArray2D plXslice =  ndarray_alloc(fGrid->NgridX(), TBRSSS::Naux());
    TNDArray2D plYslice =  ndarray_alloc(fGrid->NgridY(), TBRSSS::Naux());
    TNDArray2D plZslice =  ndarray_alloc(fGrid->NgridZ(), TBRSSS::Naux());

    TNDArray2D prXslice =  ndarray_alloc(fGrid->NgridX(), TBRSSS::Naux());
    TNDArray2D prYslice =  ndarray_alloc(fGrid->NgridY(), TBRSSS::Naux());
    TNDArray2D prZslice =  ndarray_alloc(fGrid->NgridZ(), TBRSSS::Naux());

    TNDArray1D dqnsrc =  ndarray_alloc(TBRSSS::Ncharge());


    // Compute the flux in difference in the x direction
    direction = 'x';
    for (int j=fGrid->firstY(); j<fGrid->lastY(); j++) {
    for (int k=fGrid->firstZ(); k<fGrid->lastZ(); k++) {
        pnXslice.deep() = pn[view()(j)(k)()];
        //AM checked numerically
        conservedFlux(direction, taun0, pnXslice,
                fGrid->firstX(), fGrid->lastX(),
                dqnXslice, plXslice, prXslice);
        qs[view()(j)(k)()].deep() += aex*lamx*dqnXslice ;
    }
    }

    // Compute the flux in difference in the y direction
    direction = 'y';
    for (int i=fGrid->firstX(); i<fGrid->lastX(); i++) {
    for (int k=fGrid->firstZ(); k<fGrid->lastZ(); k++) {
        pnYslice.deep() = pn[view(i)()(k)()];
        //AM checked numerically
        conservedFlux(direction,taun0, pnYslice,
                fGrid->firstY(), fGrid->lastY(),
                dqnYslice, plYslice, prYslice);
        qs[view(i)()(k)()].deep() += aex*lamy*dqnYslice ;
    }
    }

    // Compute the flux in difference in the z direction
    direction = 'z';
    for (int i=fGrid->firstX(); i<fGrid->lastX(); i++) {
    for (int j=fGrid->firstY(); j<fGrid->lastY(); j++) {
        pnZslice.deep() = pn[view(i)(j)()()];
        conservedFlux(direction,taun0, pnZslice,
                fGrid->firstZ(), fGrid->lastZ(),
                dqnZslice, plZslice, prZslice);
        qs[view(i)(j)()()].deep() += aex*lamz*dqnZslice ;
    }
    }


    // Add the source terms
    int ok;
    for (int i=fGrid->firstX(); i<fGrid->lastX(); i++) {
    for (int j=fGrid->firstY(); j<fGrid->lastY(); j++) {
    for (int k=fGrid->firstZ(); k<fGrid->lastZ(); k++) {
        const auto qs_ijk = qs[view(i)(j)(k)()] ;
        const auto pn1_ijk = pn1[view(i)(j)(k)()] ;


        // Evaluate the geometric source terms
        // AM numerically checked
        findGeometricSource(taun0, pn, i, j, k, dqnsrc);

        qs_ijk.deep() += aex*dt*dqnsrc;

        ok = fSolver->invert(taun1, qs_ijk, pn1_ijk) ;

        if (!ok) {
            // the aimdt=zero signifies this is inversion not implicit solve
            handleBadState(kstage, taun0, pn, i, j, k, 0., qs, taun1, pn1) ;
        }

        // For the viscous code we are finished
        if (fIsIdealFluid) {
            continue ;
        }

        // Evaluate the viscous source terms
        findViscousSource(kstage, taun0, pn, taun1, pn1, i, j, k, dqnsrc);

        qs_ijk.deep() += aex*dt*dqnsrc;;

        ok = fSolver->solveSource(aim*dt, taun1, qs_ijk, pn1_ijk ) ;

        if (!ok) {
            handleBadState(kstage, taun0, pn, i, j, k, aim*dt, qs, taun1, pn1) ;
        }
    }
    }
    }

}

//------------------ stepper stage functions ----------------------------------
//------------------ KT scheme conserved fluxes -------------------------------
//! q^{n+1} = q^* + a_ex*dt*(F(p^n)+S_ex(p^n)) + a_im*dt*S_im(p^n)
//! dqn = (F(p^n)+S_ex(p^n))
void TBRSSSStep::conservedFlux(const char &direction,
        const double &tau,
        const TNDArray2D &p,
        const int  &ilo,
        const int  &ihi,
        const TNDArray2D &dH,
        const TNDArray2D &pl,
        const TNDArray2D &pr)
{
    double ap, am;
    TNDArray1D Fp = ndarray_alloc(TBRSSS::Ncharge());
    TNDArray1D Fm = ndarray_alloc(TBRSSS::Ncharge());
    TNDArray1D qp = ndarray_alloc(TBRSSS::Ncharge());
    TNDArray1D qm = ndarray_alloc(TBRSSS::Ncharge());

    // reconstruct to the cell interfaces.
    reconstruct(ilo-1, ihi+1, p, pl, pr);

    // After this loop dH(i,:) contains the flux from
    // the interface between the (i-1)-th and i-th cell.
    for (int i=ilo; i<ihi+1; i++){
        propagationVelocity(direction, pl[i], pr[i-1], ap, am);
        flux(direction, tau, pr[i-1], Fp);
        flux(direction, tau, pl[i],   Fm);
        fModel->charges(tau,pr[i-1], qp);
        fModel->charges(tau,pl[i],   qm);
        for (size_t m=0; m<TBRSSS::Ncharge(); m++){
            dH(i,m) = ( ap*Fp[m] - am*Fm[m] )/(ap-am) + ap*am/(ap-am)*( qm[m] - qp[m] );
        }
    }

    // Now subtract compute the flux difference  - (Flux(i+1)  - Flux(i))
    // and store this in  dH(i).
    for (int i=ilo; i<ihi; i++) {
        for (size_t m=0 ;  m<TBRSSS::Ncharge() ; m++) {
            dH(i,m) -= dH(i+1,m);
        }
    }
    // Prevent poorly evaluated flux difference from being added to
    // the boundary points
    dH[ihi].deep() = 0. ;

}

void TBRSSSStep::propagationVelocity(const char &direction, const TNDArray1D &pl, const TNDArray1D &pr, double &ap, double &am)
{
    double apl, aml;
    double apr, amr;
    switch(fPropagationFlag) {
    case EPropagation::kIdealHydro:
        idealPropagationVelocity(direction, pl, apl,aml);
        idealPropagationVelocity(direction, pr, apr,amr);
        ap = GSL_MAX(GSL_MAX(apl,apr), 0.);
        am = GSL_MIN(GSL_MIN(aml,amr), 0.);
        if (abs(ap)>1.0 or abs(am)>1.0){
        cout << "**TBRSSSStep::propagationVelocity*** superluminal velocity!" << endl;
            std::cout << ap << " " << am << std::endl;
            abort() ;
        }
        break;
    case EPropagation::kSpeedOfLight:
        ap =1.01 ;
        am =-1.01 ;
        break;
    default:
        cout << "**TBRSSSStep::propagationVelocity*** bad selection in switch, fPropagationType is specified!" << endl;
        abort() ;
   }
}

void TBRSSSStep::idealPropagationVelocity(const char &direction, const TNDArray1D &aux, double &ap, double &am)
{
    double ut = fModel->getU0(aux)  ;
    double uk;
    switch (direction){
        case 'x':
            uk = fModel->getU1(aux)  ;
            break;
        case 'y':
            uk = fModel->getU2(aux)  ;
            break;
        case 'z':
            uk = fModel->getU3(aux)  ;
            break;
        default:
            clog << "** TBRSSSStep::idealPropagationVelocity ** Unrecognized direction \n" << direction << endl ;
            abort() ;
    }
    double cs = fModel->getCs(aux) ;
    const double cs2 = cs*cs;
    const double A = ut*uk*(1.-cs2);
    const double B = (ut*ut-uk*uk-(ut*ut-uk*uk-1.)*cs2)*cs2;
    const double D = ut*ut*(1.-cs2)+cs2;
    ap = (A+sqrt(B))/D;
    am = (A-sqrt(B))/D;
}


void TBRSSSStep::reconstruct(const int  &ilo,
        const int  &ihi,
        const TNDArray2D &p,
        const TNDArray2D &pl,
        const TNDArray2D &pr)
{
    TLimiter lim;
    using ndarray::view;
    double dphalf;

    // Extrapolate to the edges
    for (size_t m=0; m<TBRSSS::Ncharge(); m++){
        auto pv = p[view()(m)];
        for (int i=ilo; i<ihi; i++){
            dphalf = 0.5*lim.D(pv, i);
            pl(i,m) = p(i,m) - dphalf;
            pr(i,m) = p(i,m) + dphalf;
        }
    }
    // Reset the state
    for (int i=ilo; i<ihi; i++){
        fModel->fillaux(pl[i]);
        fModel->fillaux(pr[i]);
    }
}

void TBRSSSStep::flux(const char &direction, const double &tau, const TNDArray1D &auxi, const TNDArray1D &f)
{
    switch(direction){
        case 'x':
            fModel->fluxX(tau,auxi,f);
            break;
        case 'y':
            fModel->fluxY(tau,auxi,f);
            break;
        case 'z':
            fModel->fluxZ(tau,auxi,f);
            break;
        default:
            clog << "** TBRSSSStep::flux ** Unrecognized direction \n" << direction << endl ;
            abort() ;
    }
}

//------------------ stepper stage functions ----------------------------------
//------------------------- sources -------------------------------------------
void TBRSSSStep::findGeometricSource(const double &tau,
        const TNDArray4D &aux4D,
        const int &i, const int &j, const int &k,
        const TNDArray1D &src)
{
    using namespace TBRSSSConsts;

    TNDArray1D aux = aux4D[i][j][k];
    double e = fModel->getE(aux)  ;
    double p = fModel->getP(aux)  ;
    double s = fModel->getS(aux) ;
    double u[4];
    fModel->getU(aux, u) ;
    double pi[NPi44];
    fModel->getPi(aux, pi) ;

    //Source terms for conservation equations
    src[0] = - ((e+p)*u[3]*u[3] + p + pi[9]) *TBRSSS::gamma33();
    src[1] = 0. ;
    src[2] = 0. ;
    src[3] = - ((e+p)*u[0]*u[3]     + pi[3]) *TBRSSS::gamma33();

    //Geometric source terms for the pi equations
    src[kPiIJ+0] = 0. ;
    src[kPiIJ+1] = 0. ;
    src[kPiIJ+2] = -s*u[3]*pi[1]*TBRSSS::gamma33();
    src[kPiIJ+3] = 0. ;
    src[kPiIJ+4] = -s*u[3]*pi[2]*TBRSSS::gamma33();

    // Geometric source terms for the N equations
    src[kNb] = 0. ;
}

void TBRSSSStep::findViscousSource(const int &kstage,
        const double &taun0,
        const TNDArray4D &pn,
        const double &taun1,
        const TNDArray4D &pn1,
        const int &i, const int &j, const int &k,
        const TNDArray1D &src)
{
    using namespace TBRSSSConsts ;

    double dxI[3];
    dxI[0] = fGrid->getDx() ;
    dxI[1] = fGrid->getDy() ;
    dxI[2] = fGrid->getDz()*TBRSSS::g33(taun0) ;

    TNDArray1D aux = pn[i][j][k];

    double u[4];
    fModel->getU(aux, u) ;

    // ideal source is zero
    src[0] = 0  ;
    src[1] = 0. ;
    src[2] = 0. ;
    src[3] = 0. ;
    src[kNb] = 0. ;

    //AM Checked numericall (if no limiter)
    //Find local derivatives for accelleration and sigma/omega
    double DIuJ[3][3];
    for (int J=0; J<3; J++){
        for (int I=0; I<3; I++){
            DIuJ[I][J] = DIp(I, pn, i, j, k, kUI + J)/dxI[I] ;
        }
    }
    DIuJ[2][2] += TBRSSS::gamma33() * u[0]/TBRSSS::g33(taun0);


    // Extract u, DIuJ, Find acceleration, theta
    double a[4];
    double theta;
    //AM Checked numerically
    acceleration(kstage, pn, pn1, i, j, k, dxI, u, DIuJ, fDaI, a, theta);

    double pi[NPi44];
    fModel->getPi(aux, pi) ;
    
    double pia[4];
    //pia[0]=-pi[ 0 ]*a[0]+pi[ 1 ]*a[1]+pi[ 2 ]*a[2]+pi[ 3 ]*a[3];
    pia[1] = -pi[ 1 ]*a[0]+pi[ 4 ]*a[1]+pi[ 5 ]*a[2]+pi[ 6 ]*a[3];
    pia[2] = -pi[ 2 ]*a[0]+pi[ 5 ]*a[1]+pi[ 7 ]*a[2]+pi[ 8 ]*a[3];
    pia[3] = -pi[ 3 ]*a[0]+pi[ 6 ]*a[1]+pi[ 8 ]*a[2]+pi[ 9 ]*a[3];

    // upia^{IJ} = (u^I \pi^{J \mu} + u^J \pi^{I \mu})a_\mu
    double upia[NPi33];
    upia[0] = u[1]*pia[1]*2;
    upia[1] = u[1]*pia[2]+u[2]*pia[1];
    upia[2] = u[1]*pia[3]+u[3]*pia[1];
    upia[3] = u[2]*pia[2]*2;
    upia[4] = u[2]*pia[3]+u[3]*pia[2];

    double sg[NPi44];
    double om[NOmega44];
#include "TBRSSSStep_sigmaomega_util.inc"

    // Evaluate the transport parameters.
    double vc[20] ;
    fModel->getEXVSC(aux,vc) ; //Explicit Viscous Source Coefficients

    double pisigma[NPi33];
    double piomega[NPi33];
#include "TBRSSSStep_pisigma_util.inc"

    double s = fModel->getS(aux) ;
    double g33 = TBRSSS::g33(taun0) ;


    src[kPiIJ + 0] =  -g33*s*(vc[0]*sg[4] + vc[1]*pi[4]*theta + 
            vc[2]*pisigma[0] + vc[3]*piomega[0] - upia[0]);
    src[kPiIJ + 1] =  -g33*s*(vc[0]*sg[5] + vc[1]*pi[5]*theta + 
            vc[2]*pisigma[1] + vc[3]*piomega[1] - upia[1]);
    src[kPiIJ + 2] =  -g33*s*(vc[0]*sg[6] + vc[1]*pi[6]*theta + 
            vc[2]*pisigma[2] + vc[3]*piomega[2] - upia[2]);
    src[kPiIJ + 3] =  -g33*s*(vc[0]*sg[7] + vc[1]*pi[7]*theta + 
            vc[2]*pisigma[3] + vc[3]*piomega[3] - upia[3]);
    src[kPiIJ + 4] =  -g33*s*(vc[0]*sg[8] + vc[1]*pi[8]*theta + 
            vc[2]*pisigma[4] + vc[3]*piomega[4] - upia[4]);
}

double TBRSSSStep::DIp(const int &direction, const TNDArray4D &p,
      const int &i,
      const int &j,
      const int &k,
      const int &m)
{
    using ndarray::view;
    ndarray::Array<double,1> pv;
    TLimiter lim(TLimiter::kNolim);
    switch(direction){
        //x direction
        case 0:
            pv = p[view()(j)(k)(m)];
            return lim.D(pv,i);
            break;
        //y direction
        case 1:
            pv = p[view(i)()(k)(m)];
            return lim.D(pv,j);
            break;
        //z direction
        case 2:
            pv = p[view(i)(j)()(m)];
            return lim.D(pv,k);
            break;
        default:
            clog << "** TBRSSSStep::DIp ** Unrecognized direction \n" << direction << endl ;
           abort();
    }
    return 0;
}
double TBRSSSStep::Dtp(const int &kstage,
      const double &pm,
      const int &i,
      const int &j,
      const int &k,
      const int &m)
{
    switch(kstage) {
         // 2nd order time derivatve for 1st stage
        case 1: 
            {
            double dt10 = fTau_rk[kstage]-fTau_rk[kstage-1];
            double dt21 = fTau_rk[kstage+1]-fTau_rk[kstage];
            double dy10 = fP_rk[kstage](i,j,k,m) - fP_rk[kstage-1](i,j,k,m);
            double dy21 = pm - fP_rk[kstage](i,j,k,m);
            return (dy21*(dt10/dt21) + dy10*(dt21/dt10))/(dt10+dt21); 
            break;
            }
            //default 1st order time derivatve
        default:
            {
            double dt = fTau_rk[kstage+1]-fTau_rk[kstage];
            return (pm - fP_rk[kstage](i,j,k,m))/dt;
            }

    }
    return 0;
}

//! Purpose:  given the covariant derivatives  DIuJ determine provide an
//! estimated for the acceleration and the expansion scalar
//!
//! Inputs:  pn (the completed primitives of the nth stage) pn1 (an estimate
//! for the completed primitives of the n+1th stage) the coordinate points.
//!
//! Outputs: a[4] the four acceleration, theta the expansion scalar
//!
//! InOuts:  DaI, 
//!
//! On input this contains the correction (delta aI) from the pr,  on output
//! this containts an estimate for,  -u0 D_t u^I.  Later, after the step is
//! completed u0 D_t u^I will be added to this quanity (by
//! accelerationCorrection), to determine (delta aI) from the stage, kstage.
void TBRSSSStep::acceleration(const int &kstage,
      const TNDArray4D &pn,
      const TNDArray4D &pn1,
      const int &i, const int &j, const int &k,
      const double (&dxI)[3],
      const double (&u)[4],
      const double (&DIuJ)[3][3],
      const TNDArray4D &DaI,
      double (&a)[4],
      double &theta)
{
    // Just selects the functions which do the real work.
    switch (fAccelerationFlag) {
    case EAcceleration::kIdealAcceleration:
        accelerationIdeal(kstage, pn, pn1, i, j, k, dxI, u, DIuJ, DaI, a, theta);
        break;
    case EAcceleration::kNonIdealAcceleration:
        accelerationNonIdeal(kstage, pn, pn1, i, j, k, dxI, u, DIuJ, DaI, a, theta);
        break;
    default:
        clog << "*** TBRSSSS:::acceleration *** Bad selction in switch" << endl;
        exit(-1) ;
    }
}

//! Implements the acceleration (see above) using the explicit stage.
void TBRSSSStep::accelerationNonIdeal(const int &kstage,
      const TNDArray4D &pn,
      const TNDArray4D &pn1,
      const int &i, const int &j, const int &k,
      const double (&dxI)[3],
      const double (&u)[4],
      const double (&DIuJ)[3][3],
      const TNDArray4D &DaI,
      double (&a)[4],
      double &theta)
{
    // Evaluate DIuI, uIDIuJ
    double DIuI = 0.; 
    double uIDIuJ[3] = {0., 0., 0.} ;
    for (int J=0; J<3; J++) {
        for (int I=0; I<3; I++) {
            uIDIuJ[J] += u[1+I]*DIuJ[I][J];
        }
        DIuI += DIuJ[J][J];
    }

    // Extract the time derivaives:
    // covariant derivatives with  with respect to time.
    // Note (tau D_t u^eta) = partial_{\tau} (tau u^\eta)
    double DtuJ[3];
    for (int J = 0; J < 3; J++) {
       DtuJ[J] = Dtp(kstage, pn1(i,j,k,1+J), i, j, k, 1+J) ;
    }

    // Compute a^mu,  including the acorrection stored in DaI
    for (int J = 0; J < 3; J++) {
       a[1+J] =  u[0]*DtuJ[J] +  uIDIuJ[J]  + DaI(i,j,k,J);
    }
    a[0] = (a[1]*u[1]+a[2]*u[2]+a[3]*u[3])/u[0];

    // Store the a correction based on this approximate time derivative.
    // updateAccelerationCorrection will modify this value after the stage
    // to determine the difference between the actual time derivaive  of
    // the step and the approximate derivative
    for (int J = 0; J < 3; J++) {
        DaI(i,j,k,J) = -u[0]*DtuJ[J];
    }

    // Compute theta, and note D_{\tau} u^tau = partial_{\tau} u^\tau = uI
    // partial_t uI / u[0]
    theta = DIuI ;
    for (int I = 0 ; I< 3; I++) {
       theta += u[1+I]*DtuJ[I]/u[0] ;
    }
}

//! Implements acceleration (see above) using the ideal equations of motion
void TBRSSSStep::accelerationIdeal(const int &kstage,
      const TNDArray4D &pn,
      const TNDArray4D &pn1,
      const int &i, const int &j, const int &k,
      const double (&dxI)[3],
      const double (&u)[4],
      const double (&DIuJ)[3][3],
      const TNDArray4D &DaI,
      double (&a)[4],
      double &theta)
{
    // Evaluate DIuI, uIuJDIuJ and uIDIuJ
    double DIuI = 0; 
    double uIuJDIuJ = 0;
    double uIDIuJ[3] = {0., 0., 0.} ;
    for (int J=0; J<3; J++) {
        for (int I=0; I<3; I++) {
            uIDIuJ[J] += u[1+I]*DIuJ[I][J];
        }
        uIuJDIuJ += u[1+J]*uIDIuJ[J] ;
        DIuI += DIuJ[J][J];
    }

    // For extracting the time derivatives using the ideal EOM
    // we need spatial derivatives of energy density and the speed of sound
    double DIe[3];
    for (int I=0; I<3; I++){
        DIe[I] = DIp(I, pn, i, j, k, 0 ) /dxI[I];
    }
    TNDArray1D aux = pn[ndarray::view(i)(j)(k)] ;
    double e = fModel->getE(aux) ;
    double p = fModel->getP(aux) ;
    double cs = fModel->getCs(aux) ;
    double cs2 = cs*cs;

    // Construct the acceleration using the ideal EOM
    const double cs2ep = cs2/(e+p);
    const double uIDIe = (u[1]*DIe[0]+u[2]*DIe[1]+u[3]*DIe[2]);
    const double term2 = cs2/((1.-cs2)*u[0]*u[0]+cs2)*
        (cs2ep * uIDIe + uIuJDIuJ - u[0]*u[0]*DIuI);

   // We now can determine the acceleration,
    a[1] = -cs2ep*DIe[0] - u[1]*term2 + DaI(i,j,k,0);
    a[2] = -cs2ep*DIe[1] - u[2]*term2 + DaI(i,j,k,1);
    a[3] = -cs2ep*DIe[2] - u[3]*term2 + DaI(i,j,k,2);
    a[0] = (a[1]*u[1]+a[2]*u[2]+a[3]*u[3])/u[0];

    // Store the a correction based on this approximate time derivative.
    // updateAccelerationCorrection will modify this value after the stage
    // to determine the difference between the actual time derivaive  of
    // the step and the approximate derivative
    DaI(i,j,k,0) = -(-cs2ep*DIe[0] - u[1]*term2 - uIDIuJ[0]);
    DaI(i,j,k,1) = -(-cs2ep*DIe[1] - u[2]*term2 - uIDIuJ[1]);
    DaI(i,j,k,2) = -(-cs2ep*DIe[2] - u[3]*term2 - uIDIuJ[2]);

    // Extract theta using D_{\tau} u^\tau, using a[0]
    theta = a[0]/u[0]-uIuJDIuJ/(u[0]*u[0])+DIuI;
}

//--------------------- handling bad states -----------------------------------
//! Try to recover from a failures of solveSource
//!
//! First, we  a different initial guess. If this doesn't work, the bad state is
//! recorded and the output of kstage, is replaced with an average of the input
//! data of that stage.
//!
//! Outputs:
//!
//! pn1 at the best possible states we can find
//!
//! Inputs:
//!
//! The remaining arguements
void TBRSSSStep::handleBadState(const int &kstage,
        const double &taun0,
        const TNDArray4D &pn,
        const int &i, const int &j, const int &k,
        const double &aimdt,
        const TNDArray4D &qs,
        const double &taun1,
        const TNDArray4D &pn1)
{
    // Try again with a different initial guess
    find_averaged_state(i, j, k, pn, pn1[i][j][k]) ;
    int ok = fSolver->solveSource(aimdt, taun1, qs[i][j][k], pn1[i][j][k]) ;

    // It worked!
    if (ok != 0) {
       clog << " Recovered from inversion failure at " << "("<<i<<","<<j<<","<<k<<")" << endl;
    }

    // It didn't work, replace the state  with an average of the the input
    // state and log the result
    find_averaged_state(i, j, k, pn, pn1[i][j][k]) ;
    log_bad_state(kstage, aimdt, i, j, k, pn[i][j][k]) ;
}

//! Computes an average starting guess for calls to solveSource (an iterative
//! newton solver)
//!
//! Experience has shown that the average can succeed when the original
//! guess fails.
//!
void TBRSSSStep::find_averaged_state(const int &i, const int &j, const int
        &k, const TNDArray4D &pn, const TNDArray1D &pn1i)
{
    pn1i.deep() = 0;
    int ntotal = 0;
    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
                pn1i.deep() += pn[i+di][j+dj][k+dk];
                ntotal++;
            }
        }
    }
    pn1i.deep() /= ((double) ntotal);
    fModel->fillaux(pn1i) ;
}

//! This is called when all attempts and solving for the primitives
//! have failed.
//!
//! - The failure is recorded here (on screen and in a file)
//!
//! - If the number of failures excedes a certain threshold then
//!   than an exception bad_event() is thrown.
//!
//! So a typical usage is
//!
//! try {
//!    hy.step()
//! } catch (exception &bad_event) {
//!    goto the next event.
//! }
void TBRSSSStep::log_bad_state(const int &kstage, const double &aimdt, const int &i, const int &j, const int &k, const TNDArray1D &pni)
{
    // Record the location
    double tau,x,y,eta;
    tau = fGrid->getTime();
    fGrid->getXYZ(i,j,k,x,y,eta);

    // Compute quantities of interest like the temperature in GeV
    double e = fModel->getE(pni) ;
    double n = fModel->getN(pni) ;

    double s, T, mu ;
    fModel->getEOS()->stmu(e, n, s, T, mu) ;
    double T_in_gev = M_HBARC*T ;

    // Log the information to file and screen
    TIOBuffer b ;
    fInversionLog << b.format("%d %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n",
            kstage,
            aimdt, tau, x, y,
            eta, T_in_gev, e) ;

    clog << b.format("Bad inversion: kstage=%d, implict dt=%4f, at tau=%4f [i,j,k]=(%d,%d,%d) and (%3f,%3f,%3f) and T (GeV) =%3f,e (GeV/fm**3)=%3f \n",
            kstage, aimdt,
            tau,
            i,j, k,
            x, y, eta,
            T_in_gev, e*M_HBARC) ;

    // If temperature excedes 90\, MeV  then this is a real bad state
    // and we record it
    const double Tmax_not_so_bad_in_gev =0.09;
    if (T_in_gev < Tmax_not_so_bad_in_gev) {
        fBadStates++ ;
    } else {
        fRealBadStates++ ;
    }

    // Allow a certain number of bad states, before terminating the event.
    const int bad_states_max = 50000 ;
    if (fBadStates> bad_states_max) {
        cerr << "The number of bad states excedes "<< bad_states_max <<". Terminating event."<< endl;
        throw bad_event() ;
    }

    // Allow a certain number of real bad states, before terminating the event.
    const int real_bad_states_max = 4 ;
    if (fRealBadStates > real_bad_states_max) {
        cerr << "The number of real bad states excedes "<< real_bad_states_max <<". Terminating event."<< endl;
        throw bad_event();
    }
}
