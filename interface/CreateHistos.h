#ifndef __CREATEHISTOS__
#define __CREATEHISTOS__

#include "ntuple.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

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

 vector<TString> cats ={
"inclusive"
};

 // vector<TString> vars = {"m_vis","mt_1","jpt_1","jpt_2","jeta_1","jeta_2","jdeta","mjj","jeta1eta2"};
 vector<TString> vars = {
"jeta_1"
};
                        

  vector< vector<TString> > files = {{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_mt_v3.root","Z"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_mt_v3_TauPtUp.root","ZtauUp"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_mt_v3_TauPtDown.root","ZtauDown"},
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_mt_v3.root","TT"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_mt_v3_TauPtUp.root","TTtauUp"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_mt_v3_TauPtDown.root","TTtauDown"},
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_WXJets_merged_MCSpring16_mt_v3.root","W"},
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_VV_MCSpring16_mt_v3.root","VV"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_VV_MCSpring16_mt_v3_TauPtUp.root","VVtauUp"},
                                    //{"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_VV_MCSpring16_mt_v3_TauPtDown.root","VVtauDown"},
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_SingleMuonRun2016BCD_mt_v3.root","data"},                                    
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSpring16_reHLT_160921_mt_v3.root","qqH"},
                                    {"/data/higgs/data_2016/ntuples_v3/mt/ntuples_woSVFIT_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSpring16_reHLT_160921_mt_v3.root","ggH"}};


  CreateHistos();
  ~CreateHistos();


  void DYSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void VVSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void TSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void WSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void dataSelections(float var, float weight, TString cat, TString strVar, TString fname);
  void signalSelections(float var, float weight, TString cat, TString strVar, TString fname);

  float GetWRatio();
  float GetQCDRatio();
 
  void ExtractQCD(float SS_OS_ratio, TString cat, TString strVar);

  float CalcJdeta();
  float CalcHPt();

  int Common(TString sign, TString cat);
  int Vetos();
  
  int CR_W();
  int SR_QCD();
  int CR_QCD();

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


  TFile *outfile;
  vector<TH1D*> histos;
  vector<TString> histo_names = {};

  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  
 private:
  ntuple *NtupleView;

};



#endif
