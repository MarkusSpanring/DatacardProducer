#include "interface/FFCalculator.h"

using namespace std;

FFCalculator::FFCalculator(){  
}

FFCalculator::~FFCalculator(){
}

void FFCalculator::initFakeFactors(){
  for(auto cat : cats){
    if( !this->is1DCategories(cat) && !this->is2DCategories(cat) ) continue;
    TString catSuffix = cat;
    if(cat == s_wjets_0jet_cr) catSuffix = s_0jet;
    if(cat == s_wjets_boosted_cr) catSuffix = s_boosted;
    if(cat == s_wjets_vbf_cr) catSuffix = s_vbf;
    if(cat == s_antiiso_0jet_cr) catSuffix = s_0jet;
    if(cat == s_antiiso_boosted_cr) catSuffix = s_boosted;
    if(cat == s_antiiso_vbf_cr) catSuffix = s_vbf;
    FFfile[cat] = TFile::Open("HTTutilities/Jet2TauFakes/data/"+channel+"/"+catSuffix+"/"+FFversion);
    FFObj[cat] = (FakeFactor*)FFfile[cat]->Get("ff_comb");
  }
  FFsyst["mt"] = Parameter.FFsystematics.mt.syst;
  FFsyst["et"] = Parameter.FFsystematics.et.syst;
  FFsyst["tt"] = Parameter.FFsystematics.tt.syst;
}

void FFCalculator::applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  float usedVar=var;
  if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isJEC==0){
    if( channel != "tt" ){
      if( this->Baseline("FF",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFFInputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, (1.+Parameter.FFsystematics.XSvariation.DYXSvariation)*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, (1.-Parameter.FFsystematics.XSvariation.DYXSvariation)*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, (1.+Parameter.FFsystematics.XSvariation.TTXSvariation)*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, (1.+Parameter.FFsystematics.XSvariation.TTXSvariation)*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
            TString tmp=syst;
            this->getCorrectUncertaintyString( tmp );
            if(doMC && tmp.Contains("qcd")) continue;
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs, syst) );
          }
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
      if( this->Baseline("FF1",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_1 < 6 ){
          FFinputs.clear();
          this->getFF1Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, 1.08*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, 0.92*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, 1.12*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, 0.88*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
          TString tmp=syst;
          this->getCorrectUncertaintyString( tmp );
          if(doMC && tmp.Contains("qcd")) continue;
          this->GetHistbyName( fname+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs, syst) );
          }
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, 1.08*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, 0.92*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_up+sub,strVar)->Fill(usedVar, 1.12*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+s_down+sub,strVar)->Fill(usedVar, 0.88*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
            TString tmp=syst;
            this->getCorrectUncertaintyString( tmp );
            if(doMC && tmp.Contains("qcd")) continue;
            this->GetHistbyName( fname+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs, syst) );
          }
        }
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(isJEC==1){
    if( channel != "tt" ){
      if( this->Baseline("FF",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFFInputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
      if( this->Baseline("FF1",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_1 < 6 ){
          FFinputs.clear();
          this->getFF1Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if(isJEC==-1){
    if( channel != "tt" ){
      if( this->Baseline("FF",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFFInputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
      if( this->Baseline("FF1",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_1 < 6 ){
          FFinputs.clear();
          this->getFF1Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_"+s_jetFakes+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
  }  

}

void FFCalculator::getFFInputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( this->getMT() );
  inputs.push_back( NtupleView->iso_1 );
}
void FFCalculator::getFF1Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_1 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}
void FFCalculator::getFF2Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}


void FFCalculator::doUpDownReplace( TString &replaceString ){

  replaceString.ReplaceAll("_down",s_down);
  replaceString.ReplaceAll("_up",s_up);
  
}

void FFCalculator::doStatUncertaintyReplace( TString &replaceString ){

  if( replaceString.Contains("ff_qcd")
      || replaceString.Contains("ff_w") ) replaceString.ReplaceAll("stat",channel+"_stat");
  
}

void FFCalculator::doSystUncertaintyReplace( TString &replaceString ){

  if( replaceString.Contains("ff_qcd") ) replaceString.ReplaceAll("syst",channel+"_syst");
  if( channel == "tt"
      && ( replaceString.Contains("ff_w")
           || replaceString.Contains("ff_tt") ) ) replaceString.ReplaceAll("syst",channel+"_syst");
  
}
 
void FFCalculator::getCorrectUncertaintyString( TString &replaceString ){

  this->doUpDownReplace( replaceString );
  this->doStatUncertaintyReplace( replaceString );
  this->doSystUncertaintyReplace( replaceString);
  
}
