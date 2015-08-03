//This code is based on TGlauberMC_v1.1
//Original code can be found at
//http://www.hepforge.org/downloads/tglaubermc
//write up/citation reference: http://arxiv.org/abs/0805.4411
#ifndef RNAVIER_TPHOBOSMC_H
#define RNAVIER_TPHOBOSMC_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TF1.h>
#include <complex>
#endif


#include "TGlbEvent.h"

//---------------------------------------------------------------------------------
class TPhbsNucleus : public TNamed
{
   private:
      Int_t      fN;          //Number of nucleons
      Double_t   fR;          //Parameters of function
      Double_t   fA;          //Parameters of function
      Double_t   fW;          //Parameters of function
      Double_t   fMinDist;    //Minimum separation distance
      Int_t      fF;          //Type of radial distribution
      Int_t      fTrials;     //Store trials needed to complete nucleus
      TF1*       fFunction;   //Probability density function rho(r)
      TObjArray* fNucleons;   //Array of nucleons
      void       Lookup(const std::string &name);

   public:
      TPhbsNucleus(const std::string &iname, Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
      virtual ~TPhbsNucleus();

      using      TObject::Draw;
      void       Draw(Double_t xs, Int_t col);
      Int_t      GetN()             const {return fN;}
      Double_t   GetR()             const {return fR;}
      Double_t   GetA()             const {return fA;}
      Double_t   GetW()             const {return fW;}
      TObjArray *GetNucleons()      const {return fNucleons;}
      Int_t      GetTrials()        const {return fTrials;}
      void       SetN(Int_t in)           {fN=in;}
      void       SetR(Double_t ir);
      void       SetA(Double_t ia);
      void       SetW(Double_t iw);
      void       SetMinDist(Double_t min) {fMinDist=min;}
      void       ThrowNucleons(Double_t xshift=0.);

      ClassDef(TPhbsNucleus,1)
};

//---------------------------------------------------------------------------------
//This is a small class describing rapidity distribution of colliding nucleon density
class TRapidityModel;

class TPhobosMC : public TNamed 
{
   private:
       TRapidityModel *fRapidity ;
   public:

      TPhbsNucleus *fANucleus;  //Nucleus A
      TPhbsNucleus *fBNucleus;  //Nucleus B
      Double_t     fXSect;    //Nucleon-nucleon cross section

      bool       fUseRandomB ; //Flag to delect random B (in range fBMin .. fBMax)  from 2pi b db
      Double_t     fB_MC;     //Impact parameter (b)
      Double_t     fBMin;    //Minimum impact parameter to be generated
      Double_t     fBMax;    //Maximum impact parameter to be generated
      Int_t        fNpartMin ; //The minimum number of particpants

      std::string  fEntropyModel; //Controls how the entropy  is  distributed.

      // Smear the distribution of entropy by a guassian
      double fSmearingSigma ; //The smearing width = sqrt(<x^2> + <y^2>)
      double fEntropyPerNW ;  //The entropy per participant
      double fAlphaBC ;  //The entropy per binary collision

      Int_t        fEvents;   //Number of events with at least one collision
      Int_t        fTotalEvents;  //All events within  impact parameter range
      Int_t        fNpart;   //Number of wounded  nucleons in current event
      Int_t        fNcoll;   //Number of binary collisions in current event

      // A distance beyond which a partipant can be ignored
      double ParticipantRange() { return 6.*fSmearingSigma; }
      void WriteInputs(ostream &out) ;
      void SetWeight(TGlbNucleon *nucl) ;

      Bool_t       CalcResults(Double_t bgen);
      Bool_t       CalcEvent(Double_t bgen);

   public:

      TPhobosMC() { ; }
      TPhobosMC(istream &in);

      virtual     ~TPhobosMC() ;

      void         Draw(Option_t* option);
      Double_t     GetB()               const {return fB_MC;}
      Double_t     GetBMin()            const {return fBMin;}
      Double_t     GetBMax()            const {return fBMax;}
      Int_t        GetNcoll()           const {return fNcoll;}
      Int_t        GetNpart()           const {return fNpart;}
      Int_t        GetAN() const {return fANucleus->GetN();}
      Int_t        GetBN() const {return fBNucleus->GetN();}
      TObjArray   *GetNucleonsA()       const {return fANucleus->GetNucleons(); } 
      TObjArray   *GetNucleonsB()       const {return fBNucleus->GetNucleons(); }
      Double_t     GetTotXSect()        const;
      Double_t     GetTotXSectErr()     const;
      Bool_t       NextEvent(Double_t bgen=-1);
      void         SetBmin(Double_t bmin)      {fBMin = bmin;}
      void         SetBmax(Double_t bmax)      {fBMax = bmax;}
      void         SetMinDistance(Double_t d)  {fANucleus->SetMinDist(d); fBNucleus->SetMinDist(d);}

      // Start of interface for proposed base clase
      // void fillGlbEvent(TGlbEvent &evnt) ;
      void entropy(const double &time, const double &x, const double &y,const double &z, double &s, double &n) ;
      void nextEvent(const double &bgen=-1) ;
      double getB() {return GetB(); }
      int getNpart()  {return GetNpart() ; }

      ClassDef(TPhobosMC,1)
};

class TRapidityModel : public TNamed
{
    private:
        double fEtap;
        double fEtaBeam;
        double fSigmaEta;
        double f(double const &z);
        double H(double const &z);
        void WriteInputs(ostream &out); 
    public:
     TRapidityModel() { ; }
     TRapidityModel(istream &in);
     double ForwardWeight(const double &time, const double &xo, const double &yo, const double &zo) {return H(zo)*f(zo);};
     double BackwardWeight(const double &time, const double &xo, const double &yo, const double &zo){ return H(zo)*f(-zo);};
};



#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_
#endif
#endif // _run_glauber_
#endif // RNAVIER_TPHOBOSMC_H
