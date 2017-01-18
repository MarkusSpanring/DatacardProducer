#ifndef __CREATEHISTOS__
#define __CREATEHISTOS__

#include "ntuple.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "ParameterConfig_SM.h"
#include "interface/FFCalculator.h"

class CreateHistos{
 public:




//  vector<TString> cats ={"splitlowv",
// "splithighv",
// "lowv",
// "highv",
// "inclusive",
// "PUId_tight",
// "PUId_tight_fail",
// "PUId_med",
// "PUId_med_fail",
// "PUId_loose",
// "PUId_loose_fail",
// "2jet_mvis5080_mt40",
// "2jet_mvis5080_mt40_PUId_tight",
// "2jet_mvis5080_mt40_PUId_tight_fail",
// "2jet_mvis5080_mt40_PUId_med",
// "2jet_mvis5080_mt40_PUId_med_fail",
// "2jet_mvis5080_mt40_PUId_loose",
// "2jet_mvis5080_mt40_PUId_loose_fail"};

  vector<TString> cats;
  vector<TString> vars;
  vector< vector<TString> > files; 

  CreateHistos();
  ~CreateHistos();


  void initFakeFactors();
  int is1DCategories(TString category);
  int is2DCategories(TString category);
  void DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void VVSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void signalSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend="");

  float getAntiLep_tauscaling();
  float CalcJdeta();
  float CalcHPt();
  void getFFInputs(vector<double>&inputs);
  void getFF1Inputs(vector<double>&inputs);
  void getFF2Inputs(vector<double>&inputs);

  int Baseline(TString sign, TString cat);
  int Vetos();


  double QCD_OSSS(TString cat);
  void CreateQCD_osw(TString strVar, TString cat, TString extend="");
  void CreateW(TString strVar, TString cat, TString extend="");
  void CreateQCD(TString strVar, TString cat, TString extend="");
  void Estimate_W_QCD(TString strVar, TString cat, TString extend="");
  void EstimateFF(TString strVar, TString cat, TString extend="");


  int OS_W(TString cat);
  int SS_W(TString cat);
  int relaxed_W(TString cat, TString mt);
  int SS_Low(TString cat);
  int SS_Low_relaxed(TString cat);

  int CategorySelection(TString cat, TString mtcut = "");
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

  void loadFile(TString filename);
  void run(TString isTest);
  void clearHistos();
  void writeHistos(TString channel, vector<TString> cats, vector<TString> vars);
  TH1D* GetHistbyName(TString name, TString strVar="");
  TH1D* JITHistoCreator(TString name, TString strVar);
  void returnBinning(double*, vector<double> input);
  int returnBins(vector<double> input);
  TH1D* getBinnedHisto(TString name,vector<double> input);

  double get2DVar(TString sub);
  int returnBin(vector<double> bins, double value);
  double getMT();
  double getMT2();
  double getMTTOT();
  int passMTCut();
  int passIso(TString type);


  TFile *outfile;
  vector<TH1D*> histos;
  vector<TString> histo_names = {};

  vector<Double_t> FFinputs;
  map< TString, TFile*> FFfile;
  map< TString, FakeFactor*> FFObj;
  map< TString, vector<string> > FFsyst;

  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
 private:
  ntuple *NtupleView;

};



#endif
