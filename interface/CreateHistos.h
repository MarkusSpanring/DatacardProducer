#ifndef __CREATEHISTOS__
#define __CREATEHISTOS__

#include "ntuple.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "ParameterConfig_SM.h"

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


  void DYSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void VVSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void TSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void WSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void dataSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void signalSelections(float var, float weight, TString cat, TString strVar, TString fname);

  float getAntiLep_tauscaling();
  float CalcJdeta();
  float CalcHPt();

  int Baseline(TString sign, TString cat);
  int Vetos();


  double QCD_OSSS(TString cat);
  void CreateQCD_osw(TString strVar, TString cat);
  void CreateW(TString strVar, TString cat);
  void CreateQCD(TString strVar, TString cat);
  void Estimate_W_QCD(TString strVar, TString cat);


  int OS_W(TString cat);
  int SS_W(TString cat);
  int relaxed_W(TString cat, TString mt);
  int SS_Low(TString cat);
  int SS_Low_relaxed(TString cat);

  int CategorySelection(TString cat);
  int jet2_mvis();
  int VBF_low();
  int VBF_high();
  int Jet1_low();
  int Jet1_high();
  int Jet0_low();
  int Jet0_high();

  int PUJetIdSelection(TString wp);
  float PUIdCutParamsTight(float eta);
  float PUIdCutParamsMedium(float eta);
  float PUIdCutParamsLoose(float eta);

  void loadFile(TString filename);
  void run();
  void clearHistos();
  void writeHistos(TString channel, vector<TString> cats, vector<TString> vars);
  TH1D* GetHistbyName(TString name, TString strVar="");
  TH1D* JITHistoCreator(TString name, TString strVar);
  void returnBinning(double*, vector<double> input);
  int returnBins(vector<double> input);
  TH1D* getBinnedHisto(TString name,vector<double> input);

  double getMT();
  int passMTCut();
  int passIso(TString type);


  TFile *outfile;
  vector<TH1D*> histos;
  vector<TString> histo_names = {};

  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
 private:
  ntuple *NtupleView;

};



#endif
