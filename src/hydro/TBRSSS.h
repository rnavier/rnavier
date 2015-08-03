#ifndef RN_TBRSSS_h
#define RN_TBRSSS_h

#include <fstream>
#include <array>
#include "TNDArray.h"
#include "TEOS.h"

class TRNavier3DBj ;

//! \brief This class describes the BRSSS lattice site
//! 
//! The layout of each lattice site consists of Naux() real numbers. (Naux is a
//! static integer=12) The purpose of this class is to manipulate these real
//! numbers to produce physical quantities, such as pi^{\mu\nu}, and 
//! pressure
//!
//! The Naux() real numbers stored at each lattice site are
//!
//! - <tt> aux[0..3] = (e, ux, uy, uz) </tt>
//! - <tt> aux[4..8] = (pixx, pixy, pixz, piyy, piyz) </tt>
//! - <tt> aux[9 ] = (n) </tt>
//! - <tt> aux[10..11] = (p, ut) </tt>
//!
//! The members of this class are used to get access to the grid, which should
//! never be accessed directly (if at all possible), e.g.
//!
//! double e = model->getE(aux) ;
//!
//! The first Ncharge() real numbers are the ones that are actually needed to
//! specify the hydro state, while the remaining Naux()-Ncharge() numbers are
//! computed for convenience. 
namespace TBRSSSConsts {
    const int kUI=1;
    const int kPiIJ=4;
    const int kNb=9;

    const int NPi33=5;
    const int NPi44=10;
    const int NOmega44=6;
}

class TBRSSS 
{

    private:

        TEOS  *fEOS;
        //! Toggles between cartesian and bjorken coordinates
        static bool fIsCartesian ;
        double evaluateP(const TNDArray1D &aux) ;

    public:

        static const size_t fNcharge = 10  ;
        static const size_t fNaux = 12 ;

        TBRSSS(TRNavier3DBj *rn) ; 
        TBRSSS(TEOS *eos, const bool &iscartesian) : fEOS(eos) {
            fIsCartesian = iscartesian;
        } 
        ~TBRSSS() {;}

        static constexpr size_t Ncharge() { return 10 ; }
        static constexpr size_t Naux() { return 12 ; }

        //! metric component sqrt(g33)
        static double g33(const double &tau) { return (fIsCartesian ? 1. : tau); }
        //! switches on/off Bjorken Christofell symbol Gamma^3_{31}
        static double gamma33() { return (fIsCartesian ? 0. : 1.); }

        static const TNDArray1D view(double *aux) 
        { return ndarray::external(aux, ndarray::makeVector(fNaux));}  

        static const TNDArray1D view(std::array<double,fNaux> aux) 
        { return ndarray::external(aux.data(), ndarray::makeVector(fNaux));}  

        void fillaux(const TNDArray1D &aux) ;
        double getE(const TNDArray1D &aux) const { return aux[0] ;}
        double getN(const TNDArray1D &aux) const { return aux[9] ;}
        double getP(const TNDArray1D &aux) const { return aux[10]; }
        double getCs(const TNDArray1D &aux) ;
        double getS(const TNDArray1D &aux) ;
        double getT(const TNDArray1D &aux) ;

        //EXVSC = Explicit Viscous Source Coefficients
        void getEXVSC(const TNDArray1D&aux, double *vc) ;
        //IMVSC = Implicit Viscous Source Coeffcieints
        void getIMVSC(const TNDArray1D &aux, double *vc);

        double getU0(const TNDArray1D &aux) const { return aux[11];}
        double getU1(const TNDArray1D &aux) const { return aux[1] ;}
        double getU2(const TNDArray1D &aux) const { return aux[2] ;}
        double getU3(const TNDArray1D &aux) const { return aux[3] ;}

        double getPi11(const TNDArray1D &aux) const { return aux[4] ; }
        double getPi12(const TNDArray1D &aux) const { return aux[5] ; }
        double getPi13(const TNDArray1D &aux) const { return aux[6] ; }
        double getPi22(const TNDArray1D &aux) const { return aux[7] ; } 
        double getPi23(const TNDArray1D &aux) const { return aux[8] ; } 

        double getPi33(const TNDArray1D &aux) const {
            const double pixx = getPi11(aux) ;
            const double pixy = getPi12(aux) ;
            const double pixz = getPi13(aux) ;
            const double piyy = getPi22(aux) ; 
            const double piyz = getPi23(aux) ; 
            const double ut = getU0(aux) ;
            const double ux = getU1(aux) ;
            const double uy = getU2(aux) ;
            const double uz = getU3(aux) ;
            return (-(ut*ut-ux*ux)*pixx - (ut*ut-uy*uy)*piyy+2*(ux*uy*pixy+ux*uz*pixz+uy*uz*piyz))/(ut*ut-uz*uz) ; 
        }

        double getPi00(const TNDArray1D &aux) const {
            const double pixx = getPi11(aux) ;
            const double piyy = getPi22(aux) ;
            const double pizz = getPi33(aux) ;
            return (pixx + piyy + pizz) ;
        }

        double getPi01(const TNDArray1D &aux) const {
            const double pixx = getPi11(aux) ;
            const double piyx = getPi12(aux) ; //symmetric
            const double pizx = getPi13(aux) ; //symmetric
            const double ut = getU0(aux) ;
            const double ux = getU1(aux) ;
            const double uy = getU2(aux) ;
            const double uz = getU3(aux) ;
            return (ux*pixx + uy*piyx + uz*pizx)/ut ; 
        }
        double getPi02(const TNDArray1D &aux) const {
            const double pixy = getPi12(aux) ;
            const double piyy = getPi22(aux) ;
            const double pizy = getPi23(aux) ;
            const double ut = getU0(aux) ;
            const double ux = getU1(aux) ;
            const double uy = getU2(aux) ;
            const double uz = getU3(aux) ;
            return (ux*pixy + uy*piyy + uz*pizy)/ut ; 
        }
        double getPi03(const TNDArray1D &aux) const {
            const double pixz = getPi13(aux) ;
            const double piyz = getPi23(aux) ; 
            const double pizz = getPi33(aux) ;
            const double ut = getU0(aux) ;
            const double ux = getU1(aux) ;
            const double uy = getU2(aux) ;
            const double uz = getU3(aux) ;
            return (ux*pixz + uy*piyz + uz*pizz)/ut ; 
        }

        //! Returns the full pimunu in an size 10 array, with elements
        //!
        //!  Pi00, Pi01, Pi02, Pi03, Pi11, Pi12, Pi13, Pi22, Pi23, Pi33
        void   getPi(const TNDArray1D &aux, double *pi) const;

        //! Returns umu as an array with elements u[0], u[1], u[2], u[3]
        void   getU(const TNDArray1D &aux, double *u) const;

        TEOS *getEOS() { return fEOS ; }

//        static constexpr int size_pi44() { return 10 ; }
//        static constexpr int size_om44() { return 6 ; }
//        static constexpr int size_pi33() { return 5; }
//        static constexpr int index_nb() { return 9; }
//        static constexpr int index_pi() { return 4; }

        void setaux(const double &e0, const double &n0, const double uI[3], const double piIJ[5], const TNDArray1D &aux) ;

        void setauxIdeal(const double &e0, const double &n0, const double uI[3], const TNDArray1D &aux); 

        //! Conserved flux in 'x' direction
        void fluxX(const double &tau, const TNDArray1D &auxi, const TNDArray1D &fx) ;
        //! Conserved flux in 'y' direction
        void fluxY(const double &tau, const TNDArray1D &auxi, const TNDArray1D &fy) ;
        //! Conserved flux in 'z' direction
        void fluxZ(const double &tau, const TNDArray1D &auxi, const TNDArray1D &fz) ;
        //! Calculates charges 'qi' from primitives 'auxi'
        void charges(double const &tau, const TNDArray1D &auxi, const TNDArray1D &qi) ;
        //! Calculates momentum density 'TtI' from primitives 'aux'
        void getTtI(const double &tau, const TNDArray1D &aux, double TtI[3]) ;

}  ;
#endif



