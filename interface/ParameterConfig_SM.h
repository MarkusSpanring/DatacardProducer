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
    struct jpt{
      int nbins = 70;
      double nmin = 0;
      double nmax = 350;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jpt;
    struct jeta{
      int nbins = 40;
      double nmin = -5;
      double nmax = 5;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jeta;
    struct jdeta{
      int nbins = 40;
      double nmin = 0;
      double nmax = 10;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jdeta;
    struct mt_1{
      int nbins = 40;
      double nmin = 0;
      double nmax = 40;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } mt_1;
    struct CR_mt_1{
      int nbins = 40;
      double nmin = 80;
      double nmax = 1000;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } CR_mt_1;
    struct iso_1{
      int nbins = 25;
      double nmin = 0;
      double nmax = 0.5;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } iso_1;
    struct mjj{
      int nbins = 70;
      double nmin = 0;
      double nmax = 1400;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } mjj;
    struct jeta1eta2{
      int nbins = 40;
      double nmin = -10;
      double nmax = 10;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jeta1eta2;
    
    
    

  } variable;
  struct analysisCut{
    double mTHigh = 80;
    double mTLow  =50;
  } analysisCut;

} Parameter;

