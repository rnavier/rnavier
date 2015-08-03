#ifndef RN_TBRSSSStep_h
#define RN_TBRSSSStep_h

#include<string>
#include<memory>
#include<vector>
#include "TStep3D.h"
#include "TGrid3D.h"

#include "TBRSSS.h"
#include "TNDArray.h"
class TBRSSS ;
class TRNavier3DBj ;


//! \brief The purpose of TBRSSSSolver is to provide a protocol for recovering
//! the primitives from the charges, and to solve for the primitives when
//! implicit source terms are included.
//!
//! First consider the simpler case of invert.
//!
//! TBRSSSSolver   bsolver(model)  ;
//! ..
//! bsolver.invert(taun1, qns, pn1)
//!
//! Given some charges, such as the stress tensor, and charges associated
//! with the viscous strain, invert determines the primitives, by numerically
//! solving the system of equations:
//!
//!    f(p) = q* - q(p) = 0
//!
//! The BRSSSSolver is also used by the stepper algorithim to accomplish
//! the same task but with implicit source terms. In this case
//! a call to 
//! 
//! bsolver.solveSouce(aimdt, taun1, qns, pn1)
//!
//! will solve the following generalized inversion problem:
//!
//! f = q(p^{n+1}) - q^* + a_im * dt * S_im(p^{n+1}) = 0 
//!
//! Clearly calling 
//!
//! bsolver.solveSource(0, taun1, qns, pn1) is the same as calling invert.
class TBRSSSSolver 
{
    public:
    virtual ~TBRSSSSolver() {;}
    virtual int solveSource(const double &aimdt, 
                    const double &taun1, 
                    const TNDArray1D &qns, 
                    const TNDArray1D &pn1)  = 0;

    virtual int invert(const double &taun1, 
           const TNDArray1D &qns, 
           const TNDArray1D &pn1) = 0;

    virtual void source(const double &taun1, const TNDArray1D &pn1, const TNDArray1D &srcim)=0;

} ;


//! \brief TBRSSSStep class performs step update of the gird.
//!
//! typical usage is
//! 
//! brsssstep.step(gr,dt);
//!
//! note that TBRSSSStep does not change grid time so after calling 'step', one
//! has to call TGrid3D::updateTime(const double &dt);
//! TBRSSSTep stores an array of acceleration corrections, which is carried from
//! one step to another, separately from primitives stored in TGrid3D.
class TBRSSSStep : public virtual TStep3D {

    private:

        TRNavier3DBj *fRNavier ;
        TGrid3D *fGrid ;
        TBRSSS *fModel ;

        //! Used to select the stepping algorithm
        std::string fStepType; 

        //! A log for bad states which we failed to invert.
        static std::ofstream fInversionLog; 
        // and two counters.
        int fBadStates ;
        int fRealBadStates ;

        double fCFL; 
        bool fIsIdealFluid ;

        std::vector<double> fTau_rk;
        std::vector<TNDArray4D> fP_rk;
        TNDArray4D fQstar; //!< workspace for stepStage 
        TNDArray4D fDaI; //!< space for acceleration correction 
        

        //! Pointer to solver method
        std::unique_ptr<TBRSSSSolver> fSolver ; 

//------------------ stepper algorithm ----------------------------------------
        //! 1st order IMEX Euler update
        void eulerUpdateGrid(const double &dt);
        
        //! 2nd order ARS222 update
        void ars222UpdateGrid(const double &dt);

        //! Performs single update stage in a multistage algorithm
        void stepStage(const int &nstage,
                const double &dt,
                const double &taun0,
                const double &taun1,
                const TNDArray4D &pn,
                const TNDArray4D &qs,
                const double &aex,
                const double &aim,
                const TNDArray4D &pn1);        

        //!
        void fillqs(const int nstage, 
                const double &dt,
                const TNDArray2D &aex, 
                const TNDArray2D &aim,
                const TNDArray4D &qs) ;

        //! Updates the boundary of a specified stage
        void fillBoundary(const TNDArray4D &pn1) ;

        //! Calculates correction to acceleration for the next stage
        void accelerationCorrection(const int &nstage);
        

//------------------ stepper stage functions ----------------------------------
//------------------ KT scheme conserved fluxes -------------------------------
        //! Conserved flux update in given direction
        void conservedFlux(const char &direction, 
                const double &tau, 
                const TNDArray2D &p,
                const int  &ilo, 
                const int  &ihi,
                const TNDArray2D &dH,
                const TNDArray2D &pl, 
                const TNDArray2D &pr);

        //! Reconstructs left and right cell boundaries with minmod limiter
        void reconstruct( const int  &ilo, 
                const int  &ihi,
                const TNDArray2D &p,
                const TNDArray2D &pl,
                const TNDArray2D &pr);

        //! Conserved flux in given direction (calls fluxX,fluxY, fluxZ from TBRSSS)
        void flux(const char &direction, const double &tau, const TNDArray1D &auxi, const TNDArray1D &f) ;

        //! Selects the propagation model
        enum EPropagation : int { kIdealHydro=0, kSpeedOfLight=1} ;

        //! Switch for selecting the type of propgation in propagationVelocity
        EPropagation fPropagationFlag;
        
        //! Find max and min propagation speed in given direction
        void propagationVelocity(const char &direction, const TNDArray1D &pl, const TNDArray1D &pr, double &ap, double &am);
        //! \brief Finds max and min eigenvalues of dJ^k/dJ^t where J^t are conserved
        //! charges of ideal fluid and J^k is conserved flux in given direcion
        void idealPropagationVelocity(const char &direction, const TNDArray1D &aux, double &ap, double &am);

//------------------ stepper stage functions ----------------------------------
//------------------------- sources -------------------------------------------
        //! Calculates the geometric source terms.
        void findGeometricSource(const double &tau, 
                const TNDArray4D &pn, 
                const int &i, const int &j, const int &k,
                const TNDArray1D &src);

        //! \brief Calculates the explicit part of the source terms for the 
        //! viscous equations.
        void findViscousSource(
                const int &nstage,
                const double &taun0, 
                const TNDArray4D &pn, 
                const double &taun1, 
                const TNDArray4D &pn1, 
                const int &i, const int &j, const int &k,
                const TNDArray1D &src);
        //! \brief Finds limited difference between of primitive variable m in specified
        //! direction. 
        double DIp(const int &direction, 
                const TNDArray4D &p, 
                const int &i, const int &j, const int &k,
                const int &m);

        //! Calculates time derivative at a point given value pm at nstage+1.
        double Dtp(const int &nstage, 
                const double &pm,
                const int &i, const int &j, const int &k,
                const int &m);


        //! Choices for selecting the  why the acceleration is computed
        enum EAcceleration : int {kIdealAcceleration=0, kNonIdealAcceleration=1};
        
        //! The choice of the acceleration model
        EAcceleration fAccelerationFlag ;

        void acceleration(const int &nstage, 
                 const TNDArray4D &pn,
                 const TNDArray4D &pn1, 
                 const int &i, const int &j, const int &k, 
                 const double (&dxI)[3],
                 const double (&u)[4],
                 const double (&DIuJ)[3][3],
                 const TNDArray4D &DaI,
                 double (&a)[4],
                 double &theta) ; 

        void accelerationNonIdeal(const int &nstage, 
                 const TNDArray4D &pn,
                 const TNDArray4D &pn1, 
                 const int &i, const int &j, const int &k, 
                 const double (&dxI)[3],
                 const double (&u)[4],
                 const double (&DIuJ)[3][3],
                 const TNDArray4D &DaI,
                 double (&a)[4],
                 double &theta) ; 

        void accelerationIdeal(const int &nstage, 
                 const TNDArray4D &pn,
                 const TNDArray4D &pn1, 
                 const int &i, const int &j, const int &k, 
                 const double (&dxI)[3],
                 const double (&u)[4],
                 const double (&DIuJ)[3][3],
                 const TNDArray4D &DaI,
                 double (&a)[4],
                 double &theta) ; 

//--------------------- handling bad states -----------------------------------
        void handleBadState(const int &nstage, 
                const double &taun0, 
                const TNDArray4D &pn, 
                const int &i, const int &j, const int &k, 
                const double &aimdt, 
                const TNDArray4D &qs, 
                const double &taun1, 
                const TNDArray4D &pn1) ;

        void find_averaged_state(const int &i, const int &j, const int &k, 
                const TNDArray4D &pn, 
                const TNDArray1D &pn1i) ;

        void log_bad_state(const int &nstage,
                const double &aimdt,
                const int &i, const int &j, const int &k,
                const TNDArray1D &aux);

    public:

        ~TBRSSSStep() {;}
        TBRSSSStep(TRNavier3DBj *rn, TGrid3D *gr);

        int step(TGrid3D &gr, const double &dt) override ;
        double getDt(TGrid3D &gr) override;

        TBRSSSSolver *getSolver() { return fSolver.get() ; }

};


#endif



