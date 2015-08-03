/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
//This code is based on TGlauberMC_v1.1
//Original code can be found at
//http://www.hepforge.org/downloads/tglaubermc
//write up/citation reference: http://arxiv.org/abs/0805.4411
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TNamed.h>
#include <TF1.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TEllipse.h>
#include <TH2D.h>
#include <complex>
#endif


#include "TPhobosMC.h"
#include "TIOBuffer.h"

//---------------------------------------------------------------------------------
ClassImp(TPhbsNucleus)
//---------------------------------------------------------------------------------

TPhbsNucleus::TPhbsNucleus(const string &iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
   fN(iN),fR(iR),fA(ia),fW(iw),fMinDist(-1),
   fF(0),fTrials(0),fFunction(ifunc),
   fNucleons(0)
{
   if (fN==0) {
      cout << "# Setting up nucleus " << iname << endl;
      Lookup(iname);
   }
}



TPhbsNucleus::~TPhbsNucleus()
{
   if (fNucleons) {
      delete fNucleons;
   }
   delete fFunction;
}

void TPhbsNucleus::Draw(Double_t xs, Int_t col)
{
   Double_t r = 0.5*sqrt(xs/TMath::Pi()/10.);
   TEllipse e;
   e.SetLineColor(col);
   e.SetFillColor(0);
   e.SetLineWidth(1);

   for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
      TGlbNucleon* gn = (TGlbNucleon*) fNucleons->At(i);
      e.SetLineStyle(1);
      if (gn->IsSpectator()) e.SetLineStyle(3);
      e.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
   }
}

void TPhbsNucleus::Lookup(const string &name)
{
   SetName(name.c_str());

   if      (TString(name) == "p")    {fN = 1;   fR = 0.6;   fA = 0;      fW =  0;      fF = 0;}
   else if (TString(name) == "d")    {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 1;}
   else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 3;}
   else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 4;}
   else if (TString(name) == "O")    {fN = 16;  fR = 2.608; fA = 0.513;  fW = -0.051;  fF = 1;}
   else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;  fA = 0.580;  fW = -0.233;  fF = 1;}
   else if (TString(name) == "S")    {fN = 32;  fR = 2.54;  fA = 2.191;  fW =  0.16;   fF = 2;}
   else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766; fA = 0.586;  fW = -0.161;  fF = 1;}
   else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309; fA = 0.517;  fW = -0.1308; fF = 1;}
   else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;   fA = 0.596;  fW =  0;      fF = 1;}
   else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
   else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
   else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
   else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;}
   else {
      cout << "Could not find nucleus " << name << endl;
      return;
   }

   switch (fF)
   {
      case 0: // Proton
         fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
         fFunction->SetNpx(200) ;
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction = new TF1(name.c_str(),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
         fFunction->SetNpx(200) ;
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 2: // 3pG
         fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
         fFunction->SetNpx(200) ;
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 3: // Hulthen
         fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
         fFunction->SetNpx(200) ;
         fFunction->SetParameters(1/4.38,1/.85);
         break;
      case 4: // Hulthen HIJING
         fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
         fFunction->SetNpx(200) ;
         fFunction->SetParameters(2/4.38,2/.85);
         break;
      default:
         cerr << "Could not find function type " << fF << endl;
         return;
   }
}

void TPhbsNucleus::SetR(Double_t ir)
{
   fR = ir;
   switch (fF)
   {
      case 0: // Proton
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction->SetParameter(0,fR);
         break;
      case 2: // 3pG
         fFunction->SetParameter(0,fR);
         break;
   }
}

void TPhbsNucleus::SetA(Double_t ia)
{
   fA = ia;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(1,fA);
         break;
      case 2: // 3pG
         fFunction->SetParameter(1,fA);
         break;
   }
}

void TPhbsNucleus::SetW(Double_t iw)
{
   fW = iw;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(2,fW);
         break;
      case 2: // 3pG
         fFunction->SetParameter(2,fW);
         break;
   }
}

void TPhbsNucleus::ThrowNucleons(Double_t xshift)
{
   if (fNucleons==0) {
      fNucleons=new TObjArray(fN);
      fNucleons->SetOwner();
      for(Int_t i=0;i<fN;i++) {
         TGlbNucleon *nucleon=new TGlbNucleon(); 
         fNucleons->Add(nucleon); 
      }
   } 
   
   fTrials = 0;

   Double_t sumx=0;       
   Double_t sumy=0;       
   Double_t sumz=0;       

   Bool_t hulthen = (TString(GetName())=="dh");
   if (fN==2 && hulthen) { //special treatmeant for Hulten

      Double_t r = fFunction->GetRandom()/2;
      Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
      Double_t ctheta = 2*gRandom->Rndm() - 1 ;
      Double_t stheta = sqrt(1-ctheta*ctheta);
     
      TGlbNucleon *nucleon1=(TGlbNucleon*)(fNucleons->At(0));
      TGlbNucleon *nucleon2=(TGlbNucleon*)(fNucleons->At(1));
      nucleon1->Reset();
      nucleon1->SetXYZ(r * stheta * cos(phi) + xshift,
		       r * stheta * sin(phi),
		       r * ctheta);
      nucleon2->Reset();
      nucleon2->SetXYZ(-nucleon1->GetX() + 2*xshift,
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());
      fTrials = 1;
      return;
   }

   for (Int_t i = 0; i<fN; i++) {
      TGlbNucleon *nucleon=(TGlbNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while(1) {
         fTrials++;
         Double_t r = fFunction->GetRandom();
         Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
         Double_t ctheta = 2*gRandom->Rndm() - 1 ;
         Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
         Double_t x = r * stheta * cos(phi) + xshift;
         Double_t y = r * stheta * sin(phi);      
         Double_t z = r * ctheta;      
         nucleon->SetXYZ(x,y,z);
         if(fMinDist<0) break;
         Bool_t test=1;
         for (Int_t j = 0; j<i; j++) {
            TGlbNucleon *other=(TGlbNucleon*)fNucleons->At(j);
            Double_t xo=other->GetX();
            Double_t yo=other->GetY();
            Double_t zo=other->GetZ();
            Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
                                       (y-yo)*(y-yo)+
                                       (z-zo)*(z-zo));
	       
            if(dist<fMinDist) {
               test=0;
               break;
            }
         }
         if (test) break; //found nucleuon outside of mindist
      }
           
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
   }
      
//Dont shift the center of mass 
//   if(1) { // set the centre-of-mass to be at zero (+xshift)
//      sumx = sumx/fN;  
//      sumy = sumy/fN;  
//      sumz = sumz/fN;  
//      for (Int_t i = 0; i<fN; i++) {
//         TGlbNucleon *nucleon=(TGlbNucleon*)(fNucleons->At(i));
//         nucleon->SetXYZ(nucleon->GetX()-sumx-xshift,
//                         nucleon->GetY()-sumy,
//                         nucleon->GetZ()-sumz);
//      }
//   }

}

//---------------------------------------------------------------------------------
ClassImp(TPhobosMC)
//---------------------------------------------------------------------------------

TPhobosMC::~TPhobosMC() 
{
   delete fANucleus ;
   delete fBNucleus ;
   delete fRapidity ;
}

TPhobosMC::TPhobosMC(istream &in)
{
   string NA, NB ;
   TIOBuffer b ;
   b.read_section(in, "TPhobosMC") ;
   b.getS("NucleusA", NA) ;
   b.getS("NucleusB", NB) ;
   b.getD("fXSect", fXSect) ;
   b.getB("fUseRandomB", fUseRandomB) ;
   b.getD("fB_MC", fB_MC) ;
   b.getD("fBMin", fBMin) ;
   b.getD("fBMax", fBMax) ;
   b.getI("fNpartMin", fNpartMin) ; //Required number of participants

   b.getS("fEntropyModel", fEntropyModel) ;
   b.getD("fSmearingSigma", fSmearingSigma) ;
   b.getD("fEntropyPerNW", fEntropyPerNW) ;
   b.getD("fAlphaBC", fAlphaBC) ;

   fANucleus = new TPhbsNucleus(NA.c_str()) ;
   fBNucleus = new TPhbsNucleus(NB.c_str()) ;

   fEvents = 0;   
   fTotalEvents = 0;  
   fNpart = 0;   
   fNcoll = 0;   
   
   TString name(Form("Glauber_%s_%s",fANucleus->GetName(),fBNucleus->GetName()));
   TString title(Form("Glauber %s+%s",fANucleus->GetName(),fBNucleus->GetName()));
   SetName(name);
   SetTitle(title);

   WriteInputs(clog) ;
//! read in the rapidity weight
    fRapidity = new TRapidityModel(in)  ;
}
void TPhobosMC::WriteInputs(ostream &out) 
{
   out << "# " << GetName() << endl;
   out << "# " << GetTitle() << endl;

   out << "[TPhobosMC]" << endl;
   TIOBuffer b ;
   b.writeLineS(out, fANucleus->GetName(), "NA; The name of nucleus A") ;
   b.writeLineS(out, fBNucleus->GetName(), "NB; The name of nucleus A") ;
   b.writeLineD(out, fXSect, "fXSect; The nucleon nucleon cross section") ;
   b.writeLineB(out, fUseRandomB, "fUseRandomB; selects wheter to use a random impact parameter of a fixed impact parameter set by fB_MC") ;
   b.writeLineD(out, fB_MC, "fB_MC; The value of the last impact parameter") ;
   b.writeLineD(out, fBMin, "fBMin; The minimum value of B") ;
   b.writeLineD(out, fBMax, "fBMax; The maximum value of B") ;
   b.writeLineI(out, fNpartMin, "fNpartMin; The minimum number of particpants") ;

   b.writeLineS(out, fEntropyModel, "fEntropyModel; How two model local entropy production, allowed values TwoComponent") ;
   b.writeLineD(out, fSmearingSigma, "fSmearingSigma; Entropy is smeared over scale sigma") ;
   b.writeLineD(out, fEntropyPerNW, "fEntropyPerNW; Mean entropy per participant") ;
   b.writeLineD(out, fAlphaBC, "fAlphaBC; The faction controlling the contiabution of bindary coll") ;
   out << "\n" << endl;
}

//void TPhobosMC::fillGlbEvent(TGlbEvent &glb) 
//{
//   glb.clear() ;
//   glb.setSmearingSigma(fSmearingSigma) ;
//   glb.setB(fB_MC) ;
//   for(Int_t iNucl=0;iNucl<GetAN();iNucl++) {
//      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsA()->At(iNucl);
//      if(!nucl->IsWounded()) continue ;
//      glb.addNucleon(*nucl) ;
//   }
//   for(Int_t iNucl=0;iNucl<GetBN();iNucl++) {
//      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsB()->At(iNucl);
//      if(!nucl->IsWounded()) continue ;
//      glb.addNucleon(*nucl) ;
//   }
//}

Bool_t TPhobosMC::CalcEvent(Double_t bgen)
{
   fB_MC = bgen;
   // prepare event
   fANucleus->ThrowNucleons(-bgen/2.);
   for (Int_t i = 0; i<GetAN(); i++) {
      TGlbNucleon *nucleonA=(TGlbNucleon*)(GetNucleonsA()->At(i));
      nucleonA->SetInNucleusA();
   }
   fBNucleus->ThrowNucleons(bgen/2.);
   for (Int_t i = 0; i<GetBN(); i++) {
      TGlbNucleon *nucleonB=(TGlbNucleon*)(GetNucleonsB()->At(i));
      nucleonB->SetInNucleusB();
   }

   // "ball" diameter = distance at which two balls interact
   Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

   // for each of the A nucleons in nucleus B
   for (Int_t i = 0; i<GetBN(); i++) {
      TGlbNucleon *nucleonB=(TGlbNucleon*)(GetNucleonsB()->At(i));
      for (Int_t j = 0 ; j < GetAN() ;j++) {
         TGlbNucleon *nucleonA=(TGlbNucleon*)(GetNucleonsA()->At(j));
         Double_t dx = nucleonB->GetX()-nucleonA->GetX();
         Double_t dy = nucleonB->GetY()-nucleonA->GetY();
         Double_t dij = dx*dx+dy*dy;
         if (dij < d2) {
            nucleonB->Collide();
            nucleonA->Collide();
         }
      }
   }
   for(Int_t iNucl=0;iNucl<GetAN();iNucl++) {
      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsA()->At(iNucl);
      if(!nucl->IsWounded()) continue ;
      SetWeight(nucl) ;
   }
   for(Int_t iNucl=0;iNucl<GetBN();iNucl++) {
      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsB()->At(iNucl);
      if(!nucl->IsWounded()) continue ;
      SetWeight(nucl) ;
   }

   return CalcResults(bgen);
}

void TPhobosMC::SetWeight(TGlbNucleon *nucl) 
{
   double w = fEntropyPerNW*((1.-fAlphaBC)/2.+fAlphaBC/2.*nucl->GetNColl());
   nucl->SetWeight(w) ;
}

Bool_t TPhobosMC::CalcResults(Double_t bgen)
{
   // calc results for the given event
   fNpart=0;
   fNcoll=0;
   for (Int_t i = 0; i<GetAN(); i++) {
      TGlbNucleon *nucleonA=(TGlbNucleon*)(GetNucleonsA()->At(i));
      if(nucleonA->IsWounded()) {
         fNpart++;
      }
   }

   for (Int_t i = 0; i<GetBN(); i++) {
      TGlbNucleon *nucleonB=(TGlbNucleon*)(GetNucleonsB()->At(i));
      if(nucleonB->IsWounded()) {
         fNpart++;
         fNcoll += nucleonB->GetNColl();
      }
   }

   fTotalEvents++;
   if (fNpart>0) fEvents++;
   if (fNpart==0) return kFALSE;

   return kTRUE;
}


void TPhobosMC::entropy(const double &time, const double &xo, const double &yo, const double &zo, double &s, double &n) 
{
   n=  0. ;
   s = 0. ;
   // Loop over A
   for(Int_t iNucl=0;iNucl<GetAN();iNucl++) {
      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsA()->At(iNucl);
      if(!nucl->IsWounded()) continue ;

      Double_t x=nucl->GetX()  ;
      Double_t y=nucl->GetY()  ;
      double dx = x-xo ;
      double dy = y-yo ;
      double deltar = sqrt( dx*dx + dy*dy ) ;
      if (deltar > ParticipantRange()) {
         continue ;
      }
      double w = nucl->GetWeight() ; 
      double sig = fSmearingSigma/sqrt(2.) ;
      double fp = fRapidity->ForwardWeight(time,xo,yo,zo);
      s += w * 1./(2*M_PI*sig*sig) * exp( -0.5*deltar*deltar/(sig*sig) )*fp ;
   }
   // Loop over B
   for(Int_t iNucl=0;iNucl<GetBN();iNucl++) {
      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsB()->At(iNucl);
      if(!nucl->IsWounded()) continue ;

      Double_t x=nucl->GetX()  ;
      Double_t y=nucl->GetY()  ;
      double dx = x-xo ;
      double dy = y-yo ;
      double deltar = sqrt( dx*dx + dy*dy ) ;
      if (deltar > ParticipantRange()) {
         continue ;
      }
      double w = nucl->GetWeight() ; 
      double sig = fSmearingSigma/sqrt(2.) ;
      double fm = fRapidity->BackwardWeight(time,xo,yo,zo);
      s += w * 1./(2*M_PI*sig*sig) * exp( -0.5*deltar*deltar/(sig*sig) )*fm ;
   }
   // We have computed the entropy per rapidity we wish to 
   // compute the entropy density.
   s = s / time ;
   n = n / time ;
}

void TPhobosMC::Draw(Option_t* /*option*/)
{
   fANucleus->Draw(fXSect, 2);
   fBNucleus->Draw(fXSect, 4);

   TEllipse e;
   e.SetFillColor(0);
   e.SetLineColor(1);
   e.SetLineStyle(2);
   e.SetLineWidth(1);
   e.DrawEllipse(GetB()/2,0,fBNucleus->GetR(),fBNucleus->GetR(),0,360,0);
   e.DrawEllipse(-GetB()/2,0,fANucleus->GetR(),fANucleus->GetR(),0,360,0);
}


Double_t TPhobosMC::GetTotXSect() const
{
   return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100.;
}

Double_t TPhobosMC::GetTotXSectErr() const
{
   return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
      TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

Bool_t TPhobosMC::NextEvent(Double_t bgen)
{
   // If b is passed a parameter calculate the event for that B
   // If b is not passed as a parameter then if fUseRandomB is True
   // select B according to 2pi b db. s
   double b = GetB() ;
   if (bgen > 0) {
      b = bgen ;
   } else {
      if(fUseRandomB) { 
         b= TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);
      } else {
         b = GetB() ;
      }
   } 
   return CalcEvent(b);
}
void TPhobosMC::nextEvent(const double &bgen) 
{
   static const int  itrymax = 10000 ;
   int itry = 0  ;
   do { 
      if (NextEvent(bgen)) {
         if (getNpart() > fNpartMin)  return ;
      }
      itry++ ;
   } while (itry <  itrymax)  ; 
   clog << "** TPhobosMC::nextEvent ** Failed to produce a successful event ater " << itrymax << "tries. " << endl;
   exit(EXIT_FAILURE) ;
}

//////////////////////////////////////////

TRapidityModel::TRapidityModel(istream &in)
{

   SetName("TRapidityModel");
   TIOBuffer b ;
   b.read_section(in, "TRapidityModel") ;
   b.getD("fEtap", fEtap) ;
   b.getD("fEtaBeam", fEtaBeam) ;
   b.getD("fSigmaEta", fSigmaEta) ;
   WriteInputs(clog) ;
}
void TRapidityModel::WriteInputs(ostream &out) 
{
   out << "[TRapidityModel]" << endl;
   TIOBuffer b ;
   b.writeLineD(out, fEtap, "fEtap;") ;
   b.writeLineD(out, fEtaBeam, "fEtap;") ;
   b.writeLineD(out, fSigmaEta, "fSigmaEta;") ;
   out << "\n" << endl;
}
double TRapidityModel::H(double const &z)
{ 
    return ((fabs(z) < fEtap) ? 1.0 : exp(-(fabs(z) -fEtap)*(fabs(z)-fEtap)/2.0/fSigmaEta/fSigmaEta));
}

double TRapidityModel::f(double const &z)
{
    if (z <= -fEtaBeam) return 0.0;
    else if (fabs(z) < fEtaBeam) return (z+fEtaBeam)/fEtaBeam;
    else return 2.0;
}

//using namespace std ;

////! Returns the expectation values: 
////!
////! If shift == true then we  measure these moments
////! with respect to the center of mass. Otherwise these
////! are measured with respect to the unshifted center of mass.
////! 
////! rn[0] = Nparts
////! rn[1] = <r>
////! rn[2] = <r**2>
////! rn[3] = <r**3>
////! .
////!
////! zn[0] = <z>
////! zn[1] = <r**2 z>
////! zn[2] = <z^2>
//void TPhobosMC::CalcMoments(double *rn, complex<double> *zn, bool shift) 
//{
//   Double_t ave_r1 = 0. ;
//   Double_t ave_r2 = 0. ;
//   Double_t ave_r3 = 0. ;
//   Double_t ave_r4 = 0. ;
//   Double_t ave_r5 = 0. ;
//   Double_t ave_r6 = 0. ;

//   complex<double> ave_z = 0 ;
//   complex<double> ave_z1 = 0. ;   // < r^2 z >
//   complex<double> ave_z2 = 0. ;
//   complex<double> ave_z3 = 0. ;
//   complex<double> ave_z4 = 0. ;
//   complex<double> ave_z5 = 0. ;
//   complex<double> ave_z6 = 0. ;

//   for(Int_t iNucl=0;iNucl<GetAN();iNucl++) {
//      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsA()->At(iNucl);
//      if(!nucl->IsWounded()) continue ;

//      Double_t x=nucl->GetX()  ;
//      Double_t y=nucl->GetY()  ;
////      if (shift) {
////         x -= fMeanXParts ;
////         y -= fMeanYParts ;
////      }
//      Double_t r=sqrt(x*x + y*y) ;

//      ave_r1  += r ;
//      ave_r2 += r*r ;
//      ave_r3 += r*r*r ;
//      ave_r4 += r*r*r*r ;
//      ave_r5 += r*r*r*r*r ;
//      ave_r6 += r*r*r*r*r*r ;

//      complex<double> z(x,y) ;
//      ave_z  += z ;
//      ave_z1 += r*r*z ;
//      ave_z2 += z*z ;
//      ave_z3 += z*z*z ;
//      ave_z4 += z*z*z*z ;
//      ave_z5 += z*z*z*z*z ;
//      ave_z6 += z*z*z*z*z*z ;
//   }
//   for(Int_t iNucl=0;iNucl< GetBN();iNucl++) {
//      TGlbNucleon *nucl=(TGlbNucleon *) GetNucleonsB()->At(iNucl);
//      if(!nucl->IsWounded()) continue ;

//      Double_t x=nucl->GetX() ;
//      Double_t y=nucl->GetY() ;
////      if (shift) {
////         x -= fMeanXParts ;
////         y -= fMeanYParts ;
////      }
//      Double_t r=sqrt(x*x + y*y) ;

//      ave_r1 += r ;
//      ave_r2 += r*r ;
//      ave_r3 += r*r*r ;
//      ave_r4 += r*r*r*r ;
//      ave_r5 += r*r*r*r*r ;
//      ave_r6 += r*r*r*r*r*r ;

//      complex<double> z(x,y) ;
//      ave_z += z ;
//      ave_z1 += r*r*z ;
//      ave_z2 += z*z ;
//      ave_z3 += z*z*z ;
//      ave_z4 += z*z*z*z ;
//      ave_z5 += z*z*z*z*z ;
//      ave_z6 += z*z*z*z*z*z ;
//   }

//   ave_r1 /= (double) GetNpart() ;
//   ave_r2 /= (double) GetNpart() ;
//   ave_r3 /= (double)  GetNpart() ;
//   ave_r4 /= (double) GetNpart() ;
//   ave_r5 /= (double)  GetNpart() ;
//   ave_r6 /= (double)  GetNpart() ;


//   ave_z  /= (double) GetNpart() ;
//   ave_z1 /= (double) GetNpart() ;
//   ave_z2 /= (double) GetNpart() ;
//   ave_z3 /= (double) GetNpart() ;
//   ave_z4 /= (double) GetNpart() ;
//   ave_z5 /= (double) GetNpart() ; 
//   ave_z6 /= (double) GetNpart() ; 


//   rn[0] = GetNpart() ;
//   rn[1] = ave_r1 ;
//   rn[2] = ave_r2 ;
//   rn[3] = ave_r3 ; 
//   rn[4] = ave_r4 ;
//   rn[5] = ave_r5 ;
//   rn[6] = ave_r6 ;

//   zn[0] = ave_z ;
//   zn[1] = ave_z1 ;
//   zn[2] = ave_z2 ;
//   zn[3] = ave_z3 ; 
//   zn[4] = ave_z4 ;
//   zn[5] = ave_z5 ;
//   zn[6] = ave_z6 ;

//}
