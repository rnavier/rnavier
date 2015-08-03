/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
/*! \file TBRSSSSolverS2.cxx 
\brief Defines the implicit solver.

It uses gsl hybrid solvers to do the inversion. The flag
#define TBRSSSSolverS2_SOLVESOURCE_F
controls if Jacobian should be computed explicitly (if defined, it is not computed).
*/
#include "TBRSSSSolverS2.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multiroots.h"

//#define TBRSSSSolverS2_SOLVESOURCE_F
//allocate space for solver
TBRSSSSolverS2::TBRSSSSolverS2(TBRSSS *model, const bool &isideal) : fModel(model), fIsIdealFluid (isideal)  
{
    fX = gsl_vector_alloc(3) ;
    const gsl_multiroot_fsolver_type *Tf = gsl_multiroot_fsolver_hybrids;
    fFSolver = gsl_multiroot_fsolver_alloc (Tf, 3);
    const gsl_multiroot_fdfsolver_type *Tdf = gsl_multiroot_fdfsolver_hybridsj;
    fFDFSolver = gsl_multiroot_fdfsolver_alloc (Tdf, 3) ;
    fAux = ndarray_alloc(TBRSSS::Naux());
}
TBRSSSSolverS2::~TBRSSSSolverS2() 
{
    gsl_vector_free(fX) ;
    gsl_multiroot_fsolver_free(fFSolver) ;
    gsl_multiroot_fdfsolver_free(fFDFSolver) ;
}

void TBRSSSSolverS2::source(const double &taun1,
        const TNDArray1D &pn1,
        const TNDArray1D &srcim) 
{
    using namespace TBRSSSConsts; //NPi33, kPiIJ

    // Zero out all non-essential
    srcim.deep() = 0. ;

    if (fIsIdealFluid) return ;

    double s = fModel->getS(pn1) ;
    double vc[20] ;
    fModel->getIMVSC(pn1, vc) ; 

    // Solve for piIJ and fill up the fill auxiliary
    for (int i = 0 ; i < NPi33; i++) {
        srcim[kPiIJ + i] = -TBRSSS::g33(taun1)*s*( vc[0] * pn1[kPiIJ + i] );
    }
}


int tbrssssolvers2_roots_f (const gsl_vector * x, void *params, gsl_vector *f)  
{
    TBRSSS *model    = ((struct tbrssssolvers2_params *) params)->model;
    const TNDArray1D &qs  = ((struct tbrssssolvers2_params *) params)->qs;
    const double &tau  = ((struct tbrssssolvers2_params *) params)->taun1;
    const double &aimdt= ((struct tbrssssolvers2_params *) params)->aimdt;
    const TNDArray1D &aux    = ((struct tbrssssolvers2_params *) params)->aux;

    //! Takes input parameters and x and fills these variables
#include "TBRSSSSolverS2_p_util.inc" 
    //! Takes variables above and fills ff
    double ff[3] ;
    model->getTtI(tau, aux,ff);

    gsl_vector_set(f, 0, ff[0]-qs[1] );
    gsl_vector_set(f, 1, ff[1]-qs[2] );
    gsl_vector_set(f, 2, ff[2]-qs[3] );

    return GSL_SUCCESS ;
}


int tbrssssolvers2_roots_df (const gsl_vector * x, void *params, gsl_matrix *J)  
{
    TBRSSS *model    = ((struct tbrssssolvers2_params *) params)->model;
    const TNDArray1D &qs  = ((struct tbrssssolvers2_params *) params)->qs;
    const double &tau  = ((struct tbrssssolvers2_params *) params)->taun1;
    const double &aimdt= ((struct tbrssssolvers2_params *) params)->aimdt;
    const TNDArray1D &aux    = ((struct tbrssssolvers2_params *) params)->aux;
    int i, j;

    //! Takes input parameters and x and fills aux and cs
#include "TBRSSSSolverS2_p_util.inc"
    
    double JJ[3][3];
//! Takes variables above and fills JJ
#include "TBRSSSSolverS2_J_util.inc"
    for (i = 0 ; i < 3 ; i++) {
        for (j = 0 ; j < 3 ; j++) {
            gsl_matrix_set(J,i,j, JJ[i][j] ) ;
        }
    }
    return GSL_SUCCESS ;
}


int tbrssssolvers2_roots_fdf (const gsl_vector * x, void *params, gsl_vector *f, gsl_matrix *J)  
{
    TBRSSS *model    = ((struct tbrssssolvers2_params *) params)->model;
    const TNDArray1D &qs  = ((struct tbrssssolvers2_params *) params)->qs;
    const double &tau  = ((struct tbrssssolvers2_params *) params)->taun1;
    const double &aimdt= ((struct tbrssssolvers2_params *) params)->aimdt;
    const TNDArray1D &aux    = ((struct tbrssssolvers2_params *) params)->aux;
    int i, j;

    //! Takes input parameters and x and fills aux and cs
#include "TBRSSSSolverS2_p_util.inc"
    double ff[3] ;
    //! Takes variables above and fills ff
    model->getTtI(tau, aux,ff);
    gsl_vector_set(f, 0, ff[0]-qs[1] );
    gsl_vector_set(f, 1, ff[1]-qs[2] );
    gsl_vector_set(f, 2, ff[2]-qs[3] );

    double JJ[3][3];
    //! Takes variables above and fills JJ
#include "TBRSSSSolverS2_J_util.inc"
    for (i = 0 ; i < 3 ; i++) {
        for (j = 0 ; j < 3 ; j++) {
            gsl_matrix_set(J,i,j, JJ[i][j] ) ;
        }
    }
    return GSL_SUCCESS ;
}


//! solveSource is stepper class member, it has to access privite class variable
//! struct tbrss_params fSolverData to access allocated memory for the solver.
//! solveSource does the implicit time step
//! q^{n+1}-q^*+a_im*dt*S_im(p^{n+1}) = 0
//! For single stage q^*=q^n-a_ex*dt*(F(q^n)+S_ex(q^n))
int TBRSSSSolverS2::solveSource(
        const double &aimdt, //implicit time step a_im*dt
        const double &taun1,//final time at the end of implicit step tau^{n+1}
        const TNDArray1D &qns, //q^*
        const TNDArray1D &pn1) //initial guess and output of p^{n+1} 
{
    // copy in state
    fAux.deep() = pn1;
//pick solver
#ifdef TBRSSSSolverS2_SOLVESOURCE_F
    gsl_multiroot_fsolver *s = getFSolver() ;
#else
    gsl_multiroot_fdfsolver *sdf = getFDFSolver() ;
#endif
    gsl_vector *x = getX() ;
    struct tbrssssolvers2_params params = {fModel, qns, aimdt, taun1, fAux}  ;

#ifdef TBRSSSSolverS2_SOLVESOURCE_F
    gsl_multiroot_function f = {tbrssssolvers2_roots_f, 3, &params};
#else
    gsl_multiroot_function_fdf fdf= {tbrssssolvers2_roots_f, tbrssssolvers2_roots_df, tbrssssolvers2_roots_fdf, 3, &params } ;
#endif
    // set iterator initial conditions
    gsl_vector_set(x,0, fModel->getU1(fAux)) ;
    gsl_vector_set(x,1, fModel->getU2(fAux)) ;
    gsl_vector_set(x,2, fModel->getU3(fAux)) ;

    int status;
#ifdef TBRSSSSolverS2_SOLVESOURCE_F
    gsl_multiroot_fsolver_set (s, &f, x);
    status = tbrssssolvers2_gsl_multiroot_residual(s->f, 
             qns, GSL_ROOT6_DBL_MIN, GSL_DBL_EPSILON) ;
#else
    gsl_multiroot_fdfsolver_set (sdf, &fdf, x);
    status = tbrssssolvers2_gsl_multiroot_residual(sdf->f, 
             qns, GSL_ROOT6_DBL_MIN, GSL_DBL_EPSILON) ;
#endif
    //Work around GSL BUG#42220. Check that the state is not already the root
    //Calling fdfsolver_set evaluates the state. 
    if (status == GSL_SUCCESS) {
       pn1.deep() = fAux; 
       return 1;
    }

    size_t  iter = 0;

    do {
        iter++;
#ifdef TBRSSSSolverS2_SOLVESOURCE_F
        status = gsl_multiroot_fsolver_iterate (s);
#else
        status = gsl_multiroot_fdfsolver_iterate (sdf);
#endif
        if (status)  { 
            // Solver is stuck  return to sender with code 0
            // indicating erro
            return 0;
        }

#ifdef TBRSSSSolverS2_SOLVESOURCE_F
        status = tbrssssolvers2_gsl_multiroot_residual(s->f, 
                qns, epsabs(), epsrel()) ;
#else
        status = tbrssssolvers2_gsl_multiroot_residual(sdf->f, 
                qns, epsabs(), epsrel()) ;
#endif

    }
    while (status == GSL_CONTINUE && iter < 100);
    if (status !=GSL_SUCCESS) {
        return  0;
    }

    //copy back iterator state
    pn1.deep() = fAux; 
    return iter; //tbrssssolvers2_check_output_state(fModel, pnt, qi, auxi) ;

}  

//! Gets the primitives from the q's for the BRSSS Model
int TBRSSSSolverS2::invert(const double &tau, const TNDArray1D &q, const TNDArray1D &p) 
{
    return  solveSource(0.0, tau,  q, p);
}
int tbrssssolvers2_gsl_multiroot_residual(gsl_vector *f, const TNDArray1D &qs, const double &epsabs, const double &epsrel) 
{

    using TBRSSSConsts::kUI;
    double qi,ff ;
    int i ;

    for (i = 0 ; i < 3; i++) {
       ff = fabs(gsl_vector_get(f,i)) ;
       if (gsl_isnan(ff)) {
          return GSL_EBADFUNC ;
       }
    }
    for (i = 0 ; i < 3; i++) {
        ff = fabs(gsl_vector_get(f,i)) ;
        qi = fabs(qs[kUI + i]) ;
        if ( ff > (epsabs + qi*epsrel) ) {
            return GSL_CONTINUE;
        } 
    }
    return GSL_SUCCESS ;
}

void tbrssssolvers2_print_state (size_t iter, const gsl_vector *x, const gsl_vector *f)
{
    printf ("i = %3lu \n", iter) ;

    unsigned int j ;
    for (j= 0 ; j < x->size; j++) {
        printf ("x = % .20e \n", gsl_vector_get (x, j)) ;
    }
    printf("\n") ;
    for (j= 0 ; j < f->size; j++) {
        printf ("f = % .20e \n", gsl_vector_get (f, j)) ;
    }
    printf("\n") ;
}


