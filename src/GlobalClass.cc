#include "interface/GlobalClass.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TObject.h>
#include <algorithm>

using namespace std;

GlobalClass::GlobalClass()
{  
}
GlobalClass::~GlobalClass()
{
}

int GlobalClass::getNjets(){
  if(isJEC==1) return NtupleView->njets;
  else if(isJEC==-1) return NtupleView->njets;
  else return NtupleView->njets;
}

float GlobalClass::getMjj(){
  if(isJEC==1) return NtupleView->mjj;
  else if(isJEC==-1) return NtupleView->mjj;
  else return NtupleView->mjj;
}

float GlobalClass::getJdeta(){
  if(isJEC==1) return abs(NtupleView->jdeta);
  else if(isJEC==-1) return abs(NtupleView->jdeta);
  else return abs(NtupleView->jdeta);
}

float GlobalClass::getMT(){
  if(useMVAMET) return NtupleView->mt_1;
  return NtupleView->pfmt_1;
}

float GlobalClass::getMT2(){
  if(useMVAMET) return NtupleView->mt_2;
  return NtupleView->pfmt_2;
}

float GlobalClass::getMTTOT(){
  return TMath::Sqrt( TMath::Power(this->getMT(),2) + TMath::Power(this->getMT2(),2) + 2*NtupleView->pt_1*NtupleView->pt_2*( 1-TMath::Cos( TVector2::Phi_mpi_pi( NtupleView->phi_1-NtupleView->phi_2 ) ) ) );
}

float GlobalClass::CalcJdeta(){
    if(NtupleView->jeta_1 != -999 && NtupleView->jeta_2 != -999 ){
        return fabs( NtupleView->jeta_1 - NtupleView->jeta_2 );
    }
    else return -999;

}
float GlobalClass::CalcHPt(){

  TLorentzVector tau;
  TLorentzVector mu;
  TLorentzVector met;

  mu.SetPtEtaPhiM(NtupleView->pt_1, NtupleView->eta_1, NtupleView->phi_1, NtupleView->m_1);
  tau.SetPtEtaPhiM(NtupleView->pt_2, NtupleView->eta_2, NtupleView->phi_2, NtupleView->m_2);
  if(useMVAMET)met.SetPtEtaPhiM(NtupleView->mvamet,0.,NtupleView->mvametphi,0.);
  else met.SetPtEtaPhiM(NtupleView->met,0.,NtupleView->metphi,0.);

  return (tau + mu + met).Pt();
}

int GlobalClass::Baseline(TString sign, TString cat){
    
    if( NtupleView->q_1 * NtupleView->q_2 < 0
       && this->passMTCut()
       && this->passIso("base")
       && this->Vetos()
       && this->CategorySelection(cat)
        ){
      if(channel == "tt"
         && sign == "OS"
         && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_1
         && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
         ) return 1;
      else if( channel != "tt"
               && sign=="OS"
               && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
               ) return 1;
      if( sign == "FF"
          && !NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
          && NtupleView->byVLooseIsolationMVArun2v1DBoldDMwLT_2
          ) return 1;
      if( sign == "FF1"
          && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
          && !NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_1
          && NtupleView->byVLooseIsolationMVArun2v1DBoldDMwLT_1
          ) return 1;
      if( sign == "FF2"
          && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_1
          && !NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
          && NtupleView->byVLooseIsolationMVArun2v1DBoldDMwLT_2
          ) return 1;
    } 
    
   
  return 0;
}

int GlobalClass::passMTCut(){
  if( !applyMTCut ) return 1;
  if( this->getMT() < Parameter.analysisCut.mTLow ) return 1;
  return 0;
}

int GlobalClass::passIso(TString type){
  if(channel == "tt") return 1;
  if(type == "base"){
    if(channel == "et" && NtupleView->iso_1 < Parameter.analysisCut.elIso_base) return 1;
    if(channel == "mt" && NtupleView->iso_1 < Parameter.analysisCut.muIso_base) return 1;
    //if(channel == "tt" && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_1 && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2 ) return 1;
  }
  else if(type == "relaxed"){
    if(channel == "et" && NtupleView->iso_1 < Parameter.analysisCut.elIso_relaxed) return 1;
    if(channel == "mt" && NtupleView->iso_1 < Parameter.analysisCut.muIso_relaxed) return 1;
  }
  return 0;
}

//FIXME: temporary fix, has to be undone as soon as tauLepVeto is correct for etau channel
int GlobalClass::Vetos(){
  if(NtupleView->passesThirdLepVeto){
    if(channel == "et" && NtupleView->againstMuonLoose3_2 && NtupleView->againstElectronTightMVA6_2 && NtupleView->passesDiElectronVeto) return 1;
    if(channel == "mt" && NtupleView->passesTauLepVetos && NtupleView->passesDiMuonVeto) return 1;
    if(channel == "tt" && NtupleView->passesTauLepVetos) return 1; 
  }
  return 0;
}

int GlobalClass::CategorySelection(TString cat, TString mtcut){

//////////////////  Vertex dependendcy ////////////////////////////////////
  if(cat == "splitlowv"
     || cat == "splithighv"
     || cat == "lowv"
     || cat == "highv"){

    if( this->getNjets() > 1
        && this->getMT() < Parameter.analysisCut.mTLow
        && NtupleView->m_vis > 50
        && NtupleView->m_vis < 80
      ){

      if(cat == "splitlowv" && NtupleView->npv < 16) return 1;
      if(cat == "splithighv" && NtupleView->npv >= 16) return 1;
      if(cat == "lowv" && NtupleView->npv < 13) return 1;
      if(cat == "highv" && NtupleView->npv > 18) return 1;
    }  
    return 0;
  }
/////////////////////////  Inclusive ////////////////////////////////////////
  if(cat == "inclusive") return 1;
  if(cat == "PUId_tight")      return  this->PUJetIdSelection("tight");
  if(cat == "PUId_tight_fail") return !this->PUJetIdSelection("tight");
  if(cat == "PUId_med")      return  this->PUJetIdSelection("medium");
  if(cat == "PUId_med_fail") return !this->PUJetIdSelection("medium");
  if(cat == "PUId_loose")      return  this->PUJetIdSelection("loose");
  if(cat == "PUId_loose_fail") return !this->PUJetIdSelection("loose");

/////////////////////////  2jet mvis5080 mt40 ///////////////////////////////
  if(cat == "2jet_mvis5080_mt40") return this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_tight")      return  this->PUJetIdSelection("tight") & this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_tight_fail") return !this->PUJetIdSelection("tight") & this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_med")      return  this->PUJetIdSelection("medium") & this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_med_fail") return !this->PUJetIdSelection("medium") & this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_loose")      return  this->PUJetIdSelection("loose") & this->jet2_mvis();
  if(cat == "2jet_mvis5080_mt40_PUId_loose_fail") return !this->PUJetIdSelection("loose") & this->jet2_mvis();

/////////////////////////  VBF low category ///////////////////////////////
  if(cat == "VBF_low") return this->VBF_low(mtcut);
  if(cat == "PUId_lo_VBF_low") return ( this->PUJetIdSelection("loose") & this->VBF_low(mtcut) );
  if(cat == "PUId_me_VBF_low") return ( this->PUJetIdSelection("medium") & this->VBF_low(mtcut) );
  if(cat == "PUId_ti_VBF_low") return ( this->PUJetIdSelection("tight") & this->VBF_low(mtcut) );

/////////////////////////  VBF high category ///////////////////////////////
  if(cat == "VBF_high") return this->VBF_high(mtcut);
  if(cat == "PUId_VBF_high") return ( this->PUJetIdSelection("tight") & this->VBF_high(mtcut) );

/////////////////////////  1Jet low category ///////////////////////////////
  if(cat == "1Jet_low") return this->Jet1_low(mtcut);
  if(cat == "PUId_1Jet_low") return ( this->PUJetIdSelection("tight") & this->Jet1_low(mtcut) );

/////////////////////////  1Jet high category ///////////////////////////////
  if(cat == "1Jet_high") return this->Jet1_high(mtcut);
  if(cat == "PUId_1Jet_high") return ( this->PUJetIdSelection("tight") & this->Jet1_high(mtcut) );

/////////////////////////  2Jet low category ///////////////////////////////
  if(cat == "0Jet_low") return this->Jet0_low(mtcut);
  if(cat == "PUId_0Jet_low") return ( this->PUJetIdSelection("tight") & this->Jet0_low(mtcut) );

/////////////////////////  Jet high category ///////////////////////////////
  if(cat == "0Jet_high") return this->Jet0_high(mtcut);  
  if(cat == "PUId_0Jet_high") return ( this->PUJetIdSelection("tight") & this->Jet0_high(mtcut) );

  ///////////////////////  0jet category     ///////////////////////////////
  if(cat == "0jet") return this->Jet0(mtcut);
  
  ///////////////////////  boosted category     ///////////////////////////////
  if(cat == "boosted") return this->Boosted(mtcut);

  ///////////////////////  vbf category     ///////////////////////////////
  if(cat == "vbf") return this->VBF(mtcut);
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::jet2_mvis(){
  
  if( this->getNjets() > 1
      && NtupleView->jpt_1 > 30
      && NtupleView->jpt_2 > 30
      && this->getMT() < Parameter.analysisCut.mTLow
      && NtupleView->m_vis > 50
      && NtupleView->m_vis < 80)  return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::VBF_low(TString mtcut){


  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
      && this->getNjets() == 2
      && NtupleView->pt_2 > 20
      && this->getMjj() > 500
      && (this->getMjj() < 800
          || this->CalcHPt() < 100
          )
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::VBF_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
      && this->getNjets() == 2
      && NtupleView->pt_2 > 20
      && this->getMjj() > 800
      && this->CalcHPt() > 100
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Jet1_low(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
      && (this->getNjets() == 1 
          || (this->getNjets() == 2
              && this->getMjj() < 500
              ) 
          )
      && ( (NtupleView->pt_2 > 30
            && NtupleView->pt_2 <40
            )
           || (NtupleView->pt_2 > 40
               && this->CalcHPt() < 140
               )
           )
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Jet1_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
      && (this->getNjets() == 1 
          || (this->getNjets() == 2
              && this->getMjj() < 500
              ) 
          )
      && NtupleView->pt_2 > 40
      && this->CalcHPt() > 140             
      
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Jet0_low(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
       )
      && this->getNjets() == 0
      && NtupleView->pt_2 > 20
      && NtupleView->pt_2 < 50
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Jet0_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
      && this->getNjets() == 0
      && NtupleView->pt_2 > 50
      )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Jet0(TString mtcut){

  if( channel != "tt" ){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && this->getNjets() == 0
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if( this->getNjets() == 0
        && NtupleView->pt_1 > 50
        )return 1;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::Boosted(TString mtcut){

  if( channel != "tt" ){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && (this->getNjets() == 1
            || (this->getNjets()==2 && this->getMjj()<300)
            || this->getNjets()>2
            )
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if(NtupleView->pt_1 > 50 
       && ( this->getNjets() == 1
            || (this->getNjets()>=2 && !( abs(this->getJdeta()) > 2.5 && this->CalcHPt() > 100 ) )
            )
       )return 1;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::VBF(TString mtcut){

  if( channel != "tt"){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && this->getNjets() == 2
        && this->getMjj() > 300
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if(NtupleView->pt_1 > 50
       && this->getNjets() >= 2
       && abs(this->getJdeta()) > 2.5
       && this->CalcHPt() > 100
       )return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GlobalClass::PUJetIdSelection(TString wp){
  if( wp == "tight" && NtupleView->jmva_1 > this->PUIdCutParamsTight(NtupleView->jeta_1) )    return 1;
  if( wp == "medium" && NtupleView->jmva_1 > this->PUIdCutParamsMedium(NtupleView->jeta_1) )  return 1;
  if( wp == "loose" && NtupleView->jmva_1 > this->PUIdCutParamsLoose(NtupleView->jeta_1) )    return 1;
  return 0;
}

float GlobalClass::PUIdCutParamsTight(float eta){
  if( fabs(eta) < 2.5) return 0.62;
  else if( fabs(eta) < 2.75 ) return -0.21;
  else if( fabs(eta) < 3.0 )  return -0.07;
  else if( fabs(eta) < 5.0 ) return  -0.03;

  return -1.;
}

float GlobalClass::PUIdCutParamsMedium(float eta){
  if( fabs(eta) < 2.5) return -0.06;
  else if( fabs(eta) < 2.75 ) return -0.42;
  else if( fabs(eta) < 3.0 )  return -0.3;
  else if( fabs(eta) < 5.0 ) return  -0.23;

  return -1.;
}

float GlobalClass::PUIdCutParamsLoose(float eta){
  if( fabs(eta) < 2.5) return -0.92;
  else if( fabs(eta) < 2.75 ) return -0.56;
  else if( fabs(eta) < 3.0 )  return -0.44;
  else if( fabs(eta) < 5.0 ) return  -0.39;

  return -1.;
}

TH1D* GlobalClass::GetHistbyName(TString name, TString strVar){

   for(int i = 0;i<histo_names.size();i++){
      if(histo_names.at(i) == name) return histos.at(i);
   }

   return this->JITHistoCreator(name, strVar);

}

void GlobalClass::returnBinning(double* returnHisto, vector<double> input){
  for(int i=0; i<input.size(); i++) returnHisto[i]=input.at(i);
}

int GlobalClass::returnBins(vector<double> input){
  return input.size();
}

TH1D* GlobalClass::getBinnedHisto(TString name,vector<double> input){
  double binning[this->returnBins(input)];
  this->returnBinning(binning,input);
  TH1D* tmp=new TH1D(name,"",this->returnBins(input)-1,&binning[0]);
  return tmp;
}

TH1D* GlobalClass::JITHistoCreator(TString name, TString strVar){

  int nbins = 1;
  double nmin = 0;
  double nmax = 1;

  int usingVarBins = 0;

  if( channel != "tt" && name.Contains("2D") ){
    if( name.Contains(Parameter.variable2D.D2_0Jet.name) ){
      nbins = ( Parameter.variable2D.D2_0Jet.binsX.size()-1 ) * ( Parameter.variable2D.D2_0Jet.binsY.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    else if( name.Contains(Parameter.variable2D.D2_boosted.name) ){
      nbins = (doSvfit=="SVFIT") ? ( Parameter.variable2D.D2_boosted.binsX.size()-1 ) * ( Parameter.variable2D.D2_boosted.binsY_svfit.size()-1 ) : ( Parameter.variable2D.D2_boosted.binsX.size()-1 ) * ( Parameter.variable2D.D2_boosted.binsY_mvis.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    else if( name.Contains(Parameter.variable2D.D2_vbf.name) ){
      nbins = (doSvfit=="SVFIT") ? ( Parameter.variable2D.D2_vbf.binsX.size()-1 ) * ( Parameter.variable2D.D2_vbf.binsY_svfit.size()-1 ) : ( Parameter.variable2D.D2_vbf.binsX.size()-1 ) * ( Parameter.variable2D.D2_vbf.binsY_mvis.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    histos.push_back( new TH1D(name,"", nbins, nmin, nmax  ) );
    histos.back()->Sumw2();
    histo_names.push_back(name);

    return histos.back();
  }
  if( channel == "tt" && name.Contains("2D") ){
    if( name.Contains(Parameter.variable2D_tt.D2_0Jet.name) ){
      nbins = ( Parameter.variable2D_tt.D2_0Jet.binsX.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    else if( name.Contains(Parameter.variable2D_tt.D2_boosted.name) ){
      nbins = (doSvfit=="SVFIT") ? ( Parameter.variable2D_tt.D2_boosted.binsX.size()-1 ) * ( Parameter.variable2D_tt.D2_boosted.binsY_svfit.size()-1 ) : ( Parameter.variable2D_tt.D2_boosted.binsX.size()-1 ) * ( Parameter.variable2D_tt.D2_boosted.binsY_mvis.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    else if( name.Contains(Parameter.variable2D_tt.D2_vbf.name) ){
      nbins = (doSvfit=="SVFIT") ? ( Parameter.variable2D_tt.D2_vbf.binsX.size()-1 ) * ( Parameter.variable2D_tt.D2_vbf.binsY_svfit.size()-1 ) : ( Parameter.variable2D_tt.D2_vbf.binsX.size()-1 ) * ( Parameter.variable2D_tt.D2_vbf.binsY_mvis.size()-1 );
      nmin = 0;
      nmax = nbins;
    }
    histos.push_back( new TH1D(name,"", nbins, nmin, nmax  ) );
    histos.back()->Sumw2();
    histo_names.push_back(name);

    return histos.back();
  }

  if(strVar == "norm"){
    nbins = 4;
    nmin  = 0;
    nmax  = 4;
  }
  
  if(strVar == "m_vis"){
    if(Parameter.variable.m_vis.doVarBins) {
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.m_vis.varBins) );
    }
    else{
      nbins = Parameter.variable.m_vis.nbins;
      nmin  = Parameter.variable.m_vis.nmin;;
      nmax  = Parameter.variable.m_vis.nmax;
    }
  }

  if(strVar == "m_sv"){
    if(Parameter.variable.m_sv.doVarBins) {
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.m_sv.varBins) );
    }
    else{
      nbins = Parameter.variable.m_sv.nbins;
      nmin  = Parameter.variable.m_sv.nmin;;
      nmax  = Parameter.variable.m_sv.nmax;
    }
  }

  if(strVar == "pt_sv"){
    if(Parameter.variable.pt_sv.doVarBins) {
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.pt_sv.varBins) );
    }
    else{
      nbins = Parameter.variable.pt_sv.nbins;
      nmin  = Parameter.variable.pt_sv.nmin;;
      nmax  = Parameter.variable.pt_sv.nmax;
    }
  }

  else if(strVar == "jpt_1"
     || strVar == "jpt_2"
     || strVar == "jpt_1_2"
     || strVar == "jpt_1_2p5"
     || strVar == "jpt_1_3"
     || strVar == "jpt_2_2"
     || strVar == "jpt_2_2p5"
     || strVar == "jpt_2_3"){
    if(Parameter.variable.jpt.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.jpt.varBins) );
    }
    else{
      nbins = Parameter.variable.jpt.nbins;
      nmin  = Parameter.variable.jpt.nmin;
      nmax  = Parameter.variable.jpt.nmax;
    }
  }
  else if(strVar == "pt_1" || strVar == "pt_2"){
    if(Parameter.variable.pt.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.pt.varBins) );
    }
    else{
      nbins = Parameter.variable.pt.nbins;
      nmin  = Parameter.variable.pt.nmin;
      nmax  = Parameter.variable.pt.nmax;
    }
  }
  else if(strVar == "eta_1" || strVar == "eta_2"){
    if(Parameter.variable.eta.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.eta.varBins) );
    }
    else{
      nbins = Parameter.variable.eta.nbins;
      nmin  = Parameter.variable.eta.nmin;
      nmax  = Parameter.variable.eta.nmax;
    }
  }
  else if(strVar == "jeta_1" || strVar == "jeta_2"){
    if(Parameter.variable.jeta.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.jeta.varBins) );
    }
    else{
      nbins = Parameter.variable.jeta.nbins;
      nmin  = Parameter.variable.jeta.nmin;
      nmax  = Parameter.variable.jeta.nmax;
    }
  }
  else if(strVar == "jdeta"){
    if(Parameter.variable.jdeta.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.jdeta.varBins) ); 
    }
    else{
      nbins = Parameter.variable.jdeta.nbins;
      nmin  = Parameter.variable.jdeta.nmin;
      nmax  = Parameter.variable.jdeta.nmax;
    }
  }
  else if(strVar == "mt_1"){
    if(Parameter.variable.mt_1.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.mt_1.varBins) );
    }
    else{
      nbins = Parameter.variable.mt_1.nbins;
      nmin  = Parameter.variable.mt_1.nmin;
      nmax  = Parameter.variable.mt_1.nmax;
    }
  }
  else if(strVar == "met"){
    if(Parameter.variable.met.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.met.varBins) );
    }
    else{
      nbins = Parameter.variable.met.nbins;
      nmin  = Parameter.variable.met.nmin;
      nmax  = Parameter.variable.met.nmax;
    }
  }
  else if(strVar == "mttot"){
    if(Parameter.variable.mttot.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.mttot.varBins) );
    }
    else{
      nbins = Parameter.variable.mttot.nbins;
      nmin  = Parameter.variable.mttot.nmin;
      nmax  = Parameter.variable.mttot.nmax;
    }
  }
  else if(strVar == "Hpt"){
    if(Parameter.variable.Hpt.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.Hpt.varBins) );
    }
    else{
      nbins = Parameter.variable.Hpt.nbins;
      nmin  = Parameter.variable.Hpt.nmin;
      nmax  = Parameter.variable.Hpt.nmax;
    }
  }
  else if(strVar == "High_mt_1"){
    if(Parameter.variable.High_mt_1.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.High_mt_1.varBins) );
    }
    else{
      nbins = Parameter.variable.High_mt_1.nbins;
      nmin  = Parameter.variable.High_mt_1.nmin;
      nmax  = Parameter.variable.High_mt_1.nmax;
    }
  }
  else if(strVar == "Low_mt_1"){
    if(Parameter.variable.Low_mt_1.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.Low_mt_1.varBins) );
    }
    else{
      nbins = Parameter.variable.Low_mt_1.nbins;
      nmin  = Parameter.variable.Low_mt_1.nmin;
      nmax  = Parameter.variable.Low_mt_1.nmax;
    }
  }
  else if(strVar == "iso_1"){
    if(Parameter.variable.iso_1.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.iso_1.varBins) );
    }
    else{
      nbins = Parameter.variable.iso_1.nbins;
      nmin  = Parameter.variable.iso_1.nmin;
      nmax  = Parameter.variable.iso_1.nmax;
    }
  }
  else if(strVar == "mjj"){
    if(Parameter.variable.mjj.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.mjj.varBins) );
    }
    else{
      nbins = Parameter.variable.mjj.nbins;
      nmin  = Parameter.variable.mjj.nmin;
      nmax  = Parameter.variable.mjj.nmax;
    }
  }
  else if(strVar == "jeta1eta2"){
    if(Parameter.variable.jeta1eta2.doVarBins){
      usingVarBins = 1;
      histos.push_back( this->getBinnedHisto(name,Parameter.variable.jeta1eta2.varBins) );
    }
    else{
      nbins = Parameter.variable.jeta1eta2.nbins;
      nmin  = Parameter.variable.jeta1eta2.nmin;
      nmax  = Parameter.variable.jeta1eta2.nmax;
    }
  }
  
  //else throw std::invalid_argument( "Cannot create histo: " + name + ". Binning not found"  );

  if(!usingVarBins) histos.push_back( new TH1D(name,"", nbins, nmin, nmax  ) );
  histos.back()->Sumw2();
  histo_names.push_back(name);

  return histos.back();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GlobalClass::QCD_OSSS(TString cat){
  if(channel == "mt"){
    if(cat == "0jet"
       || cat == "0Jet_low"
       || cat == "0Jet_high")  return Parameter.QCD_OSSS.mt.ZeroJet;

    if(cat == "boosted"
       || cat == "1Jet_low"
       || cat == "1Jet_high")  return Parameter.QCD_OSSS.mt.Boosted;

    if(cat == "vbf"
       || cat == "VBF_low"
       || cat == "VBF_high")   return Parameter.QCD_OSSS.mt.VBF;
  }
  if(channel == "et"){
    if(cat == "0jet"
       || cat == "0Jet_low"
       || cat == "0Jet_high")  return Parameter.QCD_OSSS.et.ZeroJet;

    if(cat == "boosted"
       || cat == "1Jet_low"
       || cat == "1Jet_high")  return Parameter.QCD_OSSS.et.Boosted;

    if(cat == "vbf"
       || cat == "VBF_low"
       || cat == "VBF_high")   return Parameter.QCD_OSSS.et.VBF;
  }
  return 1.0;
}

int GlobalClass::OS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int GlobalClass::SS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int GlobalClass::relaxed_W(TString cat, TString mt){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && this->passIso("relaxed")
     && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
     && this->Vetos()
     && this->CategorySelection(cat, "wo")){

    if(mt == "low" && this->getMT() < Parameter.analysisCut.mTLow) return 1;
    if(mt == "high" && this->getMT() > Parameter.analysisCut.mTHigh) return 1;
  }
  return 0;
}

int GlobalClass::SS_Low(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("base")
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int GlobalClass::SS_Low_relaxed(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("relaxed")
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int GlobalClass::is1DCategories(TString category){
  for( auto cat : Parameter.category.D1categories ){
    if(category == cat) return 1;
  }
  return 0;
  
}

int GlobalClass::is2DCategories(TString category){
  for( auto cat : Parameter.category.D2categories ){
    if(category == cat) return 1;
  }
  return 0;
  
}

int GlobalClass::returnBin(vector<double> bins, double value){

  if( value<bins.at(0) ) return -1;
  else if( value > bins.at( bins.size() -1 ) ) return -1;
  for(int i=0; i<bins.size()-1; i++){
    if( value>bins.at(i) && value<bins.at(i+1) ) return i;
  }
  return -1;
}

double GlobalClass::get2DVar(TString sub){

  if( channel != "tt" ){
    if( sub.Contains(Parameter.variable2D.D2_0Jet.name) ){
      int binX=this->returnBin(Parameter.variable2D.D2_0Jet.binsX,NtupleView->pt_2);
      int binY=this->returnBin(Parameter.variable2D.D2_0Jet.binsY,NtupleView->m_vis);
      //if(binX>=0 && binY>=0) return binX+binY*(Parameter.variable2D.D2_0Jet.binsX.size()-1);
      if(binX>=0 && binY>=0) return binY+binX*(Parameter.variable2D.D2_0Jet.binsY.size()-1);
    }
    else if( sub.Contains(Parameter.variable2D.D2_boosted.name) ){
      int binX= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D.D2_boosted.binsX,NtupleView->pt_sv) : this->returnBin(Parameter.variable2D.D2_boosted.binsX,this->CalcHPt());
      int binY= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D.D2_boosted.binsY_svfit,NtupleView->m_sv) : this->returnBin(Parameter.variable2D.D2_boosted.binsY_mvis,NtupleView->m_vis); 
      if(binX>=0 && binY>=0) {
        if(doSvfit=="SVFIT") return binY+binX*(Parameter.variable2D.D2_boosted.binsY_svfit.size()-1);
        else return binY+binX*(Parameter.variable2D.D2_boosted.binsY_mvis.size()-1);
        //return binX+binY*(Parameter.variable2D.D2_boosted.binsX.size()-1);
      }
    }
    else if( sub.Contains(Parameter.variable2D.D2_vbf.name) ){
      int binX=this->returnBin(Parameter.variable2D.D2_vbf.binsX,this->getMjj());
      int binY= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D.D2_vbf.binsY_svfit,NtupleView->m_sv) : this->returnBin(Parameter.variable2D.D2_vbf.binsY_mvis,NtupleView->m_vis); 
      if(binX>=0 && binY>=0) {
        if(doSvfit=="SVFIT") return binY+binX*(Parameter.variable2D.D2_vbf.binsY_svfit.size()-1);
        else return binY+binX*(Parameter.variable2D.D2_vbf.binsY_mvis.size()-1);
        //return binX+binY*(Parameter.variable2D.D2_vbf.binsX.size()-1);
      }
    }
  }

  else if( channel == "tt" ){
    if( sub.Contains(Parameter.variable2D_tt.D2_0Jet.name) ){
      int binX= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D_tt.D2_0Jet.binsX,NtupleView->m_sv) : this->returnBin(Parameter.variable2D_tt.D2_0Jet.binsX,NtupleView->m_vis);
      if( binX>=0 ) return binX;
    }
    else if( sub.Contains(Parameter.variable2D_tt.D2_boosted.name) ){
      int binX= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D_tt.D2_boosted.binsX,NtupleView->pt_sv) : this->returnBin(Parameter.variable2D_tt.D2_boosted.binsX,this->CalcHPt());
      int binY= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D_tt.D2_boosted.binsY_svfit,NtupleView->m_sv) : this->returnBin(Parameter.variable2D_tt.D2_boosted.binsY_mvis,NtupleView->m_vis); 
      if(binX>=0 && binY>=0) {
        if(doSvfit=="SVFIT") return binY+binX*(Parameter.variable2D_tt.D2_boosted.binsY_svfit.size()-1);
        else return binY+binX*(Parameter.variable2D_tt.D2_boosted.binsY_mvis.size()-1);
        //return binX+binY*(Parameter.variable2D.D2_boosted.binsX.size()-1);
      }
    }
    else if( sub.Contains(Parameter.variable2D_tt.D2_vbf.name) ){
      int binX=this->returnBin(Parameter.variable2D_tt.D2_vbf.binsX,this->getMjj());
      int binY= (doSvfit=="SVFIT") ? this->returnBin(Parameter.variable2D_tt.D2_vbf.binsY_svfit,NtupleView->m_sv) : this->returnBin(Parameter.variable2D_tt.D2_vbf.binsY_mvis,NtupleView->m_vis); 
      if(binX>=0 && binY>=0) {
        if(doSvfit=="SVFIT") return binY+binX*(Parameter.variable2D_tt.D2_vbf.binsY_svfit.size()-1);
        else return binY+binX*(Parameter.variable2D_tt.D2_vbf.binsY_mvis.size()-1);
        //return binX+binY*(Parameter.variable2D.D2_vbf.binsX.size()-1);
      }
    }
  }
  
  else return -1;
  
}
