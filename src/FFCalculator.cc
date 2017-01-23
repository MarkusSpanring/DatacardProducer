#include "interface/FFCalculator.h"

using namespace std;

FFCalculator::FFCalculator(){  
}

FFCalculator::~FFCalculator(){
}

void FFCalculator::initFakeFactors(){
  for(auto cat : cats){
    if( !this->is1DCategories(cat) && !this->is2DCategories(cat) ) continue;
    FFfile[cat] = TFile::Open("HTTutilities/Jet2TauFakes/data/"+channel+"/"+cat+"/"+FFversion);
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
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.08*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.92*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.12*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.88*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
            TString tmp=syst; tmp.ReplaceAll("_down","Down"); tmp.ReplaceAll("_up","Up");
            this->GetHistbyName( fname+"_jetFakes_"+tmp+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs, syst) );
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
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.08*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.92*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.12*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.88*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
          TString tmp=syst; tmp.ReplaceAll("_down","Down"); tmp.ReplaceAll("_up","Up");
          this->GetHistbyName( fname+"_jetFakes_"+tmp+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs, syst) );
          }
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
          if(fname=="Z"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.08*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.92*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          if(fname=="TT"){
            this->GetHistbyName( fname+"_jetFakes_Up"+sub,strVar)->Fill(usedVar, 1.12*0.5*weight*FFObj[cat]->value(FFinputs) );
            this->GetHistbyName( fname+"_jetFakes_Down"+sub,strVar)->Fill(usedVar, 0.88*0.5*weight*FFObj[cat]->value(FFinputs) );
          }
          
          for( auto syst : FFsyst[channel] ){
            TString tmp=syst; tmp.ReplaceAll("_down","Down"); tmp.ReplaceAll("_up","Up");
            this->GetHistbyName( fname+"_jetFakes_"+tmp+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs, syst) );
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
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
      if( this->Baseline("FF1",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_1 < 6 ){
          FFinputs.clear();
          this->getFF1Inputs(FFinputs);
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
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
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, weight*FFObj[cat]->value(FFinputs) );
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else{
      if( this->Baseline("FF1",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_1 < 6 ){
          FFinputs.clear();
          this->getFF1Inputs(FFinputs);
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
        }
      }
      if( this->Baseline("FF2",cat) &&  ( is1DCategories(cat) || is2DCategories(cat) ) ){
        if( isData || NtupleView->gen_match_2 < 6 ){
          FFinputs.clear();
          this->getFF2Inputs(FFinputs);
          this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(usedVar, 0.5*weight*FFObj[cat]->value(FFinputs) );
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
