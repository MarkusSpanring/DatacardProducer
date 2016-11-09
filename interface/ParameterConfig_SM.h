#include "Settings.h"

const struct Parameter{

  struct dataset{
    TString data="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_SingleMuonRun2016BCD_"+channel+"_"+version+".root";
    TString Z="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+".root";
    TString ZtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+"_TauPtUp.root";
    TString ZtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+"_TauPtDown.root";
    TString TT="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_"+channel+"_"+version+".root";
    TString TTtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_"+channel+"_"+version+"_TauPtUp.root";
    TString TTtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_160919_"+channel+"_"+version+"_TauPtDown.root";
    TString W="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_WXJets_merged_MCSpring16_"+channel+"_"+version+".root";
    TString VV="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+".root";
    TString VVtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+"_TauPtUp.root";
    TString VVtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+"_TauPtDown.root";
    TString qqH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSpring16_reHLT_160921_"+channel+"_"+version+".root";
    TString ggH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSpring16_reHLT_160921_"+channel+"_"+version+".root";    
    
  } dataset;
  struct variable{
    
    struct m_vis{
      int nbins = 30;
      double nmin = 50;
      double nmax = 80;
      int doVarBins = 1;
      vector<double> varBins = {50,55,60,65,80};
    } m_vis;

  } variable;
  struct analysisCut{
    double mTHigh = 80;
    double mTLow  =50;
  } analysisCut;

} Parameter;

