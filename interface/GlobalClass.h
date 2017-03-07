#ifndef __GLOBALCLASS__
#define __GLOBALCLASS__

#include "ntuple.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "ParameterConfig_SM.h"

class GlobalClass{
 public:
  GlobalClass();
  ~GlobalClass();

  int is1DCategories(TString category);
  int is2DCategories(TString category);
  double get2DVar(TString sub);
  int returnBin(vector<double> bins, double value);

  int getNjets();
  float getMjj();
  float getJdeta();
  float getMT();
  float getMT2();
  float getMTTOT();
  float CalcJdeta();
  float CalcHPt();

  int Baseline(TString sign, TString cat);
  int passMTCut();
  int passIso(TString type);
  int Vetos();
  int CategorySelection(TString cat, TString mtcut = "");

  double getWSFUncertainty( TString cat );
  double getQCDSFUncertainty( TString cat );
  double getRenormScale( TString cat );
  double getZmumuWeight( TString cat );
  double applyZmumuUncertainty( TString cat );
  TString return2DString( TString cat );

  int jet2_mvis();
  int VBF_low(TString mtcut = "");
  int VBF_high(TString mtcut = "");
  int Jet1_low(TString mtcut = "");
  int Jet1_high(TString mtcut = "");
  int Jet0_low(TString mtcut = "");
  int Jet0_high(TString mtcut = "");
  int Jet0(TString mtcut = "");
  int Boosted(TString mtcut = "");
  int VBF(TString mtcut = "");
  int PUJetIdSelection(TString wp);
  float PUIdCutParamsTight(float eta);
  float PUIdCutParamsMedium(float eta);
  float PUIdCutParamsLoose(float eta);

  double QCD_OSSS(TString cat);
  int OS_W(TString cat);
  int SS_W(TString cat);
  int relaxed_W(TString cat, TString mt);
  int SS_Low(TString cat);
  int SS_Low_relaxed(TString cat);

  TH1D* GetHistbyName(TString name, TString strVar="");
  TH1D* JITHistoCreator(TString name, TString strVar);
  void returnBinning(double*, vector<double> input);
  int returnBins(vector<double> input);
  TH1D* getBinnedHisto(TString name,vector<double> input);
  void resetZeroBins(TString hist, TString var);
  void resetZeroBins(TH1D* hist);

  vector<TH1D*> histos;
  vector<TString> histo_names = {};

  vector<TString> cats;
  vector<TString> vars;
  vector< vector<TString> > files; 

  ntuple *NtupleView;
  Int_t isJEC=0; //0->no correction; 1->jecUp; -1->jecDown;
};



#endif
