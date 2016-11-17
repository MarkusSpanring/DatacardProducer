#include "Settings.h"

const struct Parameter{

  struct dataset{
    TString data="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_SingleMuonRun2016BCD_"+channel+"_"+version+".root";
    TString Z="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+".root";
    TString ZtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+"_TauPtUp.root";
    TString ZtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSpring16_"+channel+"_"+version+"_TauPtDown.root";
    TString TT="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_"+channel+"_"+version+".root";
    TString TTtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_"+channel+"_"+version+"_TauPtUp.root";
    TString TTtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_TT_powheg_MCSpring16_"+channel+"_"+version+"_TauPtDown.root";
    TString W="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_WXJets_merged_MCSpring16_"+channel+"_"+version+".root";
    TString VV="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+".root";
    TString VVtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+"_TauPtUp.root";
    TString VVtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VV_MCSpring16_"+channel+"_"+version+"_TauPtDown.root";
    TString qqH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+".root";
    TString qqHtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+"_TauPtUp.root";
    TString qqHtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+"_TauPtDown.root";
    TString ggH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+".root";
    TString ggHtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+"_TauPtUp.root";
    TString ggHtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+svfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSpring16_reHLT_"+channel+"_"+version+"_TauPtDown.root";
    
  } dataset;
  struct variable{
    
    struct m_vis{
      int nbins = 50;
      double nmin = 0;
      double nmax = 250;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } m_vis;
    struct jpt{
      int nbins = 70;
      double nmin = 0;
      double nmax = 350;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jpt;
    struct pt{
      int nbins = 30;
      double nmin = 20;
      double nmax = 120;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } pt;
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
      int nbins = 24;
      double nmin = 0;
      double nmax = 1200;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } mt_1;
    struct Low_mt_1{
      int nbins = 25;
      double nmin = 0;
      double nmax = 50;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } Low_mt_1;
    struct High_mt_1{
      int nbins = 40;
      double nmin = 80;
      double nmax = 1000;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } High_mt_1;
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
    double mTLow  = 50;
    double muIso_base = 0.15;
    double muIso_relaxed = 0.3;
    double elIso_base = 0.1;
    double elIso_relaxed = 0.3;
  } analysisCut;

  struct category{
    vector<TString> D1categories = {"inclusive","0Jet_low","0Jet_high","1Jet_low","1Jet_high","VBF_low","VBF_high"};
    //FIXME: add D2 categories
    //FIXME: add PU jet categories
  } category;

  struct FFsystematics{
    struct mt{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_stat_up","ff_qcd_stat_down","ff_w_syst_up","ff_w_syst_down","ff_w_stat_up","ff_w_stat_down","ff_tt_syst_up","ff_tt_syst_down","ff_tt_stat_up","ff_tt_stat_down"};
    } mt;
    struct et{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_stat_up","ff_qcd_stat_down","ff_w_syst_up","ff_w_syst_down","ff_w_stat_up","ff_w_stat_down","ff_tt_syst_up","ff_tt_syst_down","ff_tt_stat_up","ff_tt_stat_down"};
    } et;
    struct tt{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_stat_up","ff_qcd_stat_down"};
      //FIXME: add systematics for tautau if available
    } tt;
  } FFsystematics;
  
  struct QCD_OSSS
  {
    struct mt
    {
      double ZeroJet = 1.02;
      double Boosted = 1.22;
      double VBF = 1.13;
    }mt;
    struct et
    {
      double ZeroJet = 0.74;
      double Boosted = 1.0;
      double VBF = 1.15;
      
    }et;
  } QCD_OSSS ;
  
} Parameter;

 
