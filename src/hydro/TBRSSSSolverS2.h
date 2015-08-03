#ifndef rnavier_TBRSSSSolverS2_h
#define rnavier_TBRSSSSolverS2_h

#include "TNDArray.h"
#include "TBRSSS.h"
#include "TBRSSSStep.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multiroots.h"


//! The purpose of this class is to implement a TBRSSSSolver (see
//! documentation).
//!
//! The  current implementation is a straightforward  procedure
//! solving the for the full system of TBRSSS::Ncharge() equations
//! equations for the Ncharge() primitives.
class TBRSSSSolverS2 : public TBRSSSSolver
{
    private:

        TBRSSS *fModel ;
        gsl_vector *fX;    //!< vector to query the solver state
        TNDArray1D fAux;    //!< vector to store all primitives
        gsl_multiroot_fsolver *fFSolver; //!< gsl solver without Jacobian
        gsl_multiroot_fdfsolver *fFDFSolver; //!< gsl solver with explicit Jacobian calculation

        bool fIsIdealFluid ; 

    public:

        static constexpr const double epsabs() { return 1.e-8; }
        static constexpr const double epsrel() { return 1.e-6; }

        TBRSSSSolverS2(TBRSSS *model, const bool &isideal=false)  ;
        virtual ~TBRSSSSolverS2() ;
        int solveSource(const double &aimdt, 
                const double &taun1, 
                const TNDArray1D &qns, 
                const TNDArray1D &pn1) override ;

        int invert(const double &taun1, 
                const TNDArray1D &qns, 
                const TNDArray1D &pn1) override;

        void source(const double &taun1, 
                const TNDArray1D &pn1,
                const TNDArray1D &srcim) override;

        gsl_vector *getX() {return fX;}
        gsl_multiroot_fsolver *getFSolver () {return fFSolver;}
        gsl_multiroot_fdfsolver *getFDFSolver () {return fFDFSolver ;}

} ;

//! Parameters for tbrssssolvers2_roots_f() and its Jacobian
//!
//! f = q(p^{n+1}) - q^* + a_im * dt * S_im(p^{n+1})
//!
struct tbrssssolvers2_params{
    TBRSSS       *model ; 
    const TNDArray1D  &qs;    //!< q^*
    const double   &aimdt ;   //!< a_im*dt
    const double   &taun1 ;   //!< tau^{n+1}
    const TNDArray1D   &aux ;   //!< p at each iteration stage
} ;


//! Returns the function 
//!
//! f = q(p^{n+1}) - q^* + a_im * dt * S_im(p^{n+1})
int tbrssssolvers2_roots_f (const gsl_vector * x, void *params, gsl_vector *f) ;

//! df= df/dp^{n+1}
int tbrssssolvers2_roots_df (const gsl_vector * x, void *params, gsl_matrix *J) ;

//! f and df calculated at the same time
int tbrssssolvers2_roots_fdf (const gsl_vector * x, void *params, gsl_vector *f, gsl_matrix *J) ;

//! Prints the current values of the iteration variables
void tbrssssolvers2_print_state (size_t iter, const gsl_vector *x, const gsl_vector *f) ;

//! Checks if function has converged
int tbrssssolvers2_gsl_multiroot_residual(gsl_vector *f, const  TNDArray1D &qn, const double &epsabs, const double &epsrel)  ;

#endif

