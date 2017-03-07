#include "Settings.h"

const struct Parameter{

  struct dataset{
    TString data_mt="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_SingleMuon_"+channel+"_"+version+reduced+".root";
    TString MCsum_mt="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_MCsum_merged_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString data_et="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_SingleElectron_"+channel+"_"+version+reduced+".root";
    TString data_tt="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_Tau_"+channel+"_"+version+reduced+".root";
    TString Z="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString ZtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString ZtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_DYXJetsToLL_lowMass_merged_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";
    TString EWKZ="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_EWK_merged_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString EWKZtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_EWK_merged_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString EWKZtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_EWK_merged_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";
    TString TT="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_TT_merged_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString TTtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_TT_merged_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString TTtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_TT_merged_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";
    TString W="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_WXJets_merged_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString VV="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VV_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString VVtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VV_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString VVtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VV_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";
    TString qqH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString qqHtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString qqHtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_VBFHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";
    TString ggH="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+reduced+".root";
    TString ggHtauUp="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+"_TauPtUp"+reduced+".root";
    TString ggHtauDown="/data/higgs/data_2016/ntuples_"+version+"/"+channel+"/ntuples_"+doSvfit+"_merged/BASIS_ntuple_GluGluHToTauTau_M125_powheg_MCSummer16_"+channel+"_"+version+"_TauPtDown"+reduced+".root";

  } dataset;
  struct variable{
    
    struct m_vis{
      int nbins = 20;
      double nmin = 0;
      double nmax = 200;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } m_vis;
    struct m_sv{
      int nbins = 25;
      double nmin = 0;
      double nmax = 250;
      int doVarBins = 0;
      vector<double> varBins = {0,20,25,30};
    } m_sv;
    struct pt_sv{
      int nbins = 20;
      double nmin = 0;
      double nmax = 200;
      int doVarBins = 0;
      vector<double> varBins = {0,20,25,30};
    } pt_sv;
    struct jpt{
      int nbins = 70;
      double nmin = 0;
      double nmax = 350;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jpt;
    struct pt{
      int nbins = 25;
      double nmin = 30;
      double nmax = 100;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } pt;
    struct eta{
      int nbins = 25;
      double nmin = -2.5;
      double nmax = 2.5;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } eta;
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
      int nbins = 25;
      double nmin = 0;
      double nmax = 250;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } mt_1;
    struct met{
      int nbins = 25;
      double nmin = 0;
      double nmax = 250;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } met;
    struct mttot{
      int nbins = 25;
      double nmin = 0;
      double nmax = 250;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } mttot;
    struct Hpt{
      int nbins = 20;
      double nmin = 0;
      double nmax = 200;
      int doVarBins = 0;
      vector<double> varBins = {0,200};
    } Hpt;
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
      int doVarBins = 1;
      vector<double> varBins = {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,625,650,675,700,725,750,775,800,850,900,950,1000,1100,1200,1400};
    } mjj;
    struct jeta1eta2{
      int nbins = 40;
      double nmin = -10;
      double nmax = 10;
      int doVarBins = 0;
      vector<double> varBins = {50,55,60,65,80};
    } jeta1eta2;

  } variable;

  struct variable2D{
    struct D2_0Jet{
      TString name = s_0jet;
      vector<double> binsX = {30,35,40,45,50,55,300}; //pt_2
      vector<double> binsY = {0,60,65,70,75,80,85,90,95,100,105,110,400}; //m_vis
    } D2_0Jet;
    struct D2_boosted{
      TString name = s_boosted;
      vector<double> binsX = {0,100,150,200,250,300,5000}; //Higgs_pt
      vector<double> binsY_svfit = {0,80,90,100,110,120,130,140,150,160,300}; //svfit_mass
      vector<double> binsY_mvis = {0,60,70,80,90,100,110,120,130,140,300}; //m_vis instead of svfit_mass
      //binsY_mvis must have same number of bins
    } D2_boosted;
    struct D2_vbf{
      TString name = s_vbf;
      vector<double> binsX = {300,700,1100,1500,10000}; //m_jj
      vector<double> binsY_svfit = {0,95,115,135,155,400}; //svfit_mass
      vector<double> binsY_mvis = {0,75,95,115,135,400}; //m_vis instead of svfit_mass
      //binsY_mvis must have same number of bins
    } D2_vbf;
    struct antiiso{
      vector<double> binsX = {40,80,120,160,200};
    } antiiso;
  } variable2D;
  struct variable2D_tt{
    struct D2_0Jet{
      TString name = s_0jet;
      vector<double> binsX = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300}; //svfit_mass of m_vis
    } D2_0Jet;
    struct D2_boosted{
      TString name = s_boosted;
      vector<double> binsX = {0,100,170,300,1000}; //Higgs_pt
      vector<double> binsY_svfit = {0,40,60,70,80,90,100,110,120,130,150,200,250}; //svfit_mass
      vector<double> binsY_mvis = {0,30,50,60,70,80,90,100,110,120,140,200,250}; //m_vis instead of svfit_mass
      //binsY_mvis must have same number of bins
    } D2_boosted;
    struct D2_vbf{
      TString name = s_vbf;
      vector<double> binsX = {300,500,800,10000}; //m_jj
      vector<double> binsY_svfit = {0,40,60,70,80,90,100,110,120,130,150,200,250}; //svfit_mass
      vector<double> binsY_mvis = {0,30,50,60,70,80,90,100,110,120,140,200,250}; //m_vis instead of svfit_mass
      //binsY_mvis must have same number of bins
    } D2_vbf;
  } variable2D_tt;  

  struct analysisCut{
    double mTHigh = 80;
    double mTLow  = 50;
    double muIso_base = 0.15;
    double muIso_relaxed = 0.3;
    double muIso_antiIsoLow = 0.15;
    double muIso_antiIsoHigh = 0.3;
    double elIso_base = 0.1;
    double elIso_relaxed = 0.3;
    double elIso_antiIsoLow = 0.1;
    double elIso_antiIsoHigh = 0.3;
  } analysisCut;

  struct category{
    vector<TString> D1categories = {s_inclusive,s_0jet_low,s_0jet_high,s_1jet_low,s_1jet_high,s_vbf_low,s_vbf_high};
    vector<TString> D2categories = {s_0jet,s_boosted,s_vbf,s_wjets_0jet_cr,s_wjets_boosted_cr,s_wjets_vbf_cr,s_antiiso_0jet_cr,s_antiiso_boosted_cr,s_antiiso_vbf_cr};
    //FIXME: add PU jet categories
  } category;

  struct FFsystematics{
    struct XSvariation{
      double DYXSvariation = 0.08;
      double TTXSvariation = 0.12;
    }XSvariation;
    struct mt{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_dm0_njet0_stat_up","ff_qcd_dm0_njet0_stat_down","ff_qcd_dm0_njet1_stat_up","ff_qcd_dm0_njet1_stat_down","ff_qcd_dm1_njet0_stat_up","ff_qcd_dm1_njet0_stat_down","ff_qcd_dm1_njet1_stat_up","ff_qcd_dm1_njet1_stat_down","ff_w_syst_up","ff_w_syst_down","ff_w_dm0_njet0_stat_up","ff_w_dm0_njet0_stat_down","ff_w_dm0_njet1_stat_up","ff_w_dm0_njet1_stat_down","ff_w_dm1_njet0_stat_up","ff_w_dm1_njet0_stat_down","ff_w_dm1_njet1_stat_up","ff_w_dm1_njet1_stat_down","ff_tt_syst_up","ff_tt_syst_down","ff_tt_dm0_njet0_stat_up","ff_tt_dm0_njet0_stat_down","ff_tt_dm0_njet1_stat_up","ff_tt_dm0_njet1_stat_down","ff_tt_dm1_njet0_stat_up","ff_tt_dm1_njet0_stat_down","ff_tt_dm1_njet1_stat_up","ff_tt_dm1_njet1_stat_down"};
    } mt;
    struct et{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_dm0_njet0_stat_up","ff_qcd_dm0_njet0_stat_down","ff_qcd_dm0_njet1_stat_up","ff_qcd_dm0_njet1_stat_down","ff_qcd_dm1_njet0_stat_up","ff_qcd_dm1_njet0_stat_down","ff_qcd_dm1_njet1_stat_up","ff_qcd_dm1_njet1_stat_down","ff_w_syst_up","ff_w_syst_down","ff_w_dm0_njet0_stat_up","ff_w_dm0_njet0_stat_down","ff_w_dm0_njet1_stat_up","ff_w_dm0_njet1_stat_down","ff_w_dm1_njet0_stat_up","ff_w_dm1_njet0_stat_down","ff_w_dm1_njet1_stat_up","ff_w_dm1_njet1_stat_down","ff_tt_syst_up","ff_tt_syst_down","ff_tt_dm0_njet0_stat_up","ff_tt_dm0_njet0_stat_down","ff_tt_dm0_njet1_stat_up","ff_tt_dm0_njet1_stat_down","ff_tt_dm1_njet0_stat_up","ff_tt_dm1_njet0_stat_down","ff_tt_dm1_njet1_stat_up","ff_tt_dm1_njet1_stat_down"};
    } et;
    struct tt{
      vector<string> syst = {"ff_qcd_syst_up","ff_qcd_syst_down","ff_qcd_dm0_njet0_stat_up","ff_qcd_dm0_njet0_stat_down","ff_qcd_dm0_njet1_stat_up","ff_qcd_dm0_njet1_stat_down","ff_qcd_dm1_njet0_stat_up","ff_qcd_dm1_njet0_stat_down","ff_qcd_dm1_njet1_stat_up","ff_qcd_dm1_njet1_stat_down","ff_w_syst_up","ff_w_syst_down","ff_tt_syst_up","ff_tt_syst_down"};
    } tt;
  } FFsystematics;
  
  struct QCD_OSSS
  {
    struct mt
    {
      double ZeroJet = 1.00;
      double Boosted = 1.15;
      double VBF = 1.20;
    }mt;
    struct et
    {
      double ZeroJet = 1.00;
      double Boosted = 1.15;
      double VBF = 1.20;
      
    }et;
  } QCD_OSSS ;
  
} Parameter;

 
