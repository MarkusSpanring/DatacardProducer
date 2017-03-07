#include "interface/SelectionAnalyzer.h"

using namespace std;

SelectionAnalyzer::SelectionAnalyzer(){  
}

SelectionAnalyzer::~SelectionAnalyzer(){
}

void SelectionAnalyzer::initDYSelections(TString cat, TString strVar, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  this->GetHistbyName(s_Z+sub,strVar);
  this->GetHistbyName(s_ZL+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSdyShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSdyShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSZLShape+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSZLShape+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZLjecUp+sub,strVar);
  this->GetHistbyName(s_ZLjecDown+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSzmumu+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZL+s_CMSzmumu+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZLL+sub,strVar);
  this->GetHistbyName(s_ZJ+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSdyShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSdyShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSjetToTauFake+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSjetToTauFake+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSzmumu+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZJ+s_CMSzmumu+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZJjecUp+sub,strVar);
  this->GetHistbyName(s_ZJjecDown+sub,strVar);    
  this->GetHistbyName(s_ZTT+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMSdyShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMSdyShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMSzmumu+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ZTT+s_CMSzmumu+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ZTTjecUp+sub,strVar);
  this->GetHistbyName(s_ZTTjecDown+sub,strVar);
  
}

void SelectionAnalyzer::initEWKZSelections(TString cat, TString strVar, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  this->GetHistbyName(s_EWKZ+sub,strVar);
  this->GetHistbyName(s_EWKZ+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_EWKZ+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_EWKZjecUp+sub,strVar);
  this->GetHistbyName(s_EWKZjecDown+sub,strVar);
  
}

void SelectionAnalyzer::initTSelections(TString cat, TString strVar, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  this->GetHistbyName(s_TTJ+sub,strVar);
  this->GetHistbyName(s_TTJ+s_CMSttShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TTJ+s_CMSttShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TTJ+s_CMSjetToTauFake+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TTJ+s_CMSjetToTauFake+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TTJjecUp+sub,strVar);
  this->GetHistbyName(s_TTJjecDown+sub,strVar);    
  this->GetHistbyName(s_TTT+sub,strVar);
  this->GetHistbyName(s_TTT+s_CMSttShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TTT+s_CMSttShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TTT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TTT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TTTjecUp+sub,strVar);
  this->GetHistbyName(s_TTTjecDown+sub,strVar);
  this->GetHistbyName(s_TT+sub,strVar);
  this->GetHistbyName(s_TT+s_CMSttShape+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TT+s_CMSttShape+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_TT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_TTjecUp+sub,strVar);
  this->GetHistbyName(s_TTjecDown+sub,strVar);

}

void SelectionAnalyzer::initVVSelections(TString cat, TString strVar, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  this->GetHistbyName(s_VVJ+sub,strVar);
  this->GetHistbyName(s_VVJjecUp+sub,strVar);
  this->GetHistbyName(s_VVJjecDown+sub,strVar);    
  this->GetHistbyName(s_VVT+sub,strVar);
  this->GetHistbyName(s_VVTjecUp+sub,strVar);
  this->GetHistbyName(s_VVTjecDown+sub,strVar);
  this->GetHistbyName(s_VV+sub,strVar);
  this->GetHistbyName(s_VVjecUp+sub,strVar);
  this->GetHistbyName(s_VVjecDown+sub,strVar);

}

void SelectionAnalyzer::initSignalSelections(TString cat, TString strVar, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  this->GetHistbyName(s_qqH+"125"+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+sub,strVar);
  this->GetHistbyName(s_qqH+"125"+s_jecUp+sub,strVar);
  this->GetHistbyName(s_qqH+"125"+s_jecDown+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_jecUp+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_jecDown+sub,strVar);
  this->GetHistbyName(s_qqH+"125"+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_qqH+"125"+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_CMSscalegg+s_13TeVUp+sub,strVar);
  this->GetHistbyName(s_ggH+"125"+s_CMSscalegg+s_13TeVDown+sub,strVar);

}

void SelectionAnalyzer::DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == s_Z){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_Z+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->LSelection() ){
          this->GetHistbyName(s_ZLL+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_ZL+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_ZL+s_CMSdyShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName(s_ZL+s_CMSdyShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );
          this->GetHistbyName(s_ZL+s_CMSzmumu+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * this->applyZmumuUncertainty( cat ) );
          this->GetHistbyName(s_ZL+s_CMSzmumu+s_13TeVDown+sub,strVar)->Fill(usedVar, weight / this->applyZmumuUncertainty( cat ) );
        }
        else if( this->TSelection() ){
          this->GetHistbyName(s_ZTT+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_ZTT+s_CMSdyShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName(s_ZTT+s_CMSdyShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );
          this->GetHistbyName(s_ZTT+s_CMSzmumu+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * this->applyZmumuUncertainty( cat ) );
          this->GetHistbyName(s_ZTT+s_CMSzmumu+s_13TeVDown+sub,strVar)->Fill(usedVar, weight / this->applyZmumuUncertainty( cat ) );
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_ZLL+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_ZJ+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_ZJ+s_CMSdyShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName(s_ZJ+s_CMSdyShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );
          this->GetHistbyName(s_ZJ+s_CMSjetToTauFake+s_13TeVUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          this->GetHistbyName(s_ZJ+s_CMSjetToTauFake+s_13TeVDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
          this->GetHistbyName(s_ZJ+s_CMSzmumu+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * this->applyZmumuUncertainty( cat ) );
          this->GetHistbyName(s_ZJ+s_CMSzmumu+s_13TeVDown+sub,strVar)->Fill(usedVar, weight / this->applyZmumuUncertainty( cat ) );
        }
        ////////////////////////////////////////////////////////////////
        if( NtupleView->NUP == 0)          this->GetHistbyName("Z_0Jets"+sub,strVar)->Fill(usedVar, weight);
        else if( NtupleView->NUP == 1)     this->GetHistbyName("Z_1Jets"+sub,strVar)->Fill(usedVar, weight);  
        else if(  NtupleView->NUP > 1)     this->GetHistbyName("Z_ge2Jets"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_ZJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_ZJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_Z+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_Z+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_Z+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_Z+sub,strVar)->Fill(usedVar, weight);
     
        if( this->JSelection() ){
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
        }
        else{
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
        }

      }
      
    }
    else if(fname == s_ZtauUp){
      if( this->Baseline("OS",cat)
          && this->TSelection() ){
        this->GetHistbyName(s_ZTT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
      }
      if( this->Baseline("OS",cat)
          && this->LSelection() ){
        this->GetHistbyName(s_ZL+s_CMSZLShape+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
      }
    }
    else if(fname == s_ZtauDown){
      if( this->Baseline("OS",cat)
          && this->TSelection() ){
        this->GetHistbyName(s_ZTT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
      }
      if( this->Baseline("OS",cat)
          && this->LSelection() ){
        this->GetHistbyName(s_ZL+s_CMSZLShape+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
      }
    }
    else if(fname == s_ZjecUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_ZjecUp+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->LSelection() ){
          this->GetHistbyName(s_ZLjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->TSelection() ){
          this->GetHistbyName(s_ZTTjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_ZJjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_ZJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_ZJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_Z+s_jecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_Z+s_jecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_Z+s_jecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_Z+s_jecUp+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == s_ZjecDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_ZjecDown+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->LSelection() ){
          this->GetHistbyName(s_ZLjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->TSelection() ){
          this->GetHistbyName(s_ZTTjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_ZJjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_ZJjecDown+"_+"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_ZJjecDown+"_+"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_Z+s_jecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_Z+s_jecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_Z+s_jecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_Z+s_jecDown+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == s_EWKZ){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_EWKZ+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_EWKZ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_EWKZ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_EWKZ+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_EWKZ+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_EWKZ+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_EWKZ+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == s_EWKZtauUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_EWKZ+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
      }
    }
    else if(fname == s_EWKZtauDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_EWKZ+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
      }
    }
    else if(fname == s_EWKZjecUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_EWKZjecUp+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_EWKZjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_EWKZjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_EWKZjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_EWKZjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_EWKZjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecUp+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == s_EWKZjecDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_EWKZjecDown+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_EWKZjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_EWKZjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_"+s_EWKZjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_"+s_EWKZjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_"+s_EWKZjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecDown+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
    TString sub = extend+ "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == s_TT){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_TT+sub,strVar)->Fill(usedVar, weight);
        this->GetHistbyName(s_TT+s_CMSttShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight * NtupleView->topWeight);
        this->GetHistbyName(s_TT+s_CMSttShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
        if( this->LSelection() || this->TSelection() ){ // for TT TSelection and L Selection are combined
          this->GetHistbyName(s_TTT+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_TTT+s_CMSttShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName(s_TTT+s_CMSttShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_TTJ+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName(s_TTJ+s_CMSttShape+s_13TeVUp+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName(s_TTJ+s_CMSttShape+s_13TeVDown+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
          this->GetHistbyName(s_TTJ+s_CMSjetToTauFake+s_13TeVUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          this->GetHistbyName(s_TTJ+s_CMSjetToTauFake+s_13TeVDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_TTJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_TTJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TT+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_TT+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_TT+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_TT+sub,strVar)->Fill(usedVar, weight);

        if( this->JSelection() ){
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
        }
        else{
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
        }
        
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_TTtauUp){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName(s_TT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( this->LSelection()
               || this->TSelection() ) )             this->GetHistbyName(s_TTT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_TTtauDown){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName(s_TT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( this->LSelection()
               || this->TSelection() ) )             this->GetHistbyName(s_TTT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_TTjecUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_TTjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->LSelection() || this->TSelection() ){
          this->GetHistbyName(s_TTTjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_TTJjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_TTJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_TTJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TTjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_TTjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_TTjecUp+sub,strVar)->Fill(usedVar, weight);

        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_TTjecUp+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_TTjecDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_TTjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->LSelection() || this->TSelection() ){
          this->GetHistbyName(s_TTTjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        else if ( this->JSelection() ){
          this->GetHistbyName(s_TTJjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_TTJjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_TTJjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_TTjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_TTjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_TTjecDown+sub,strVar)->Fill(usedVar, weight);

        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_TTjecDown+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;
    
    if(fname == s_W){
      if( this->Baseline("OS",cat) ){ 
        this->GetHistbyName("SR_MC"+s_W+sub,strVar)->Fill(usedVar, weight);

        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_W+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_W+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_W+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_W+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_"+s_W+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_"+s_W+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_W+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_W+sub, strVar)->Fill(usedVar, weight);

        if( this->JSelection() ){
          if( this->Baseline("OS",cat) )      this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->Baseline("OS",cat) )      this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
          if( this->OS_W(cat) )               this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->OS_W(cat) )               this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
          if( this->relaxed_W(cat, "low") )   this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->relaxed_W(cat, "low") )   this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
          if( this->relaxed_W(cat, "high") )  this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeUp(NtupleView->pt_2) );
          if( this->relaxed_W(cat, "high") )  this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight*this->getJetToTauFakeDown(NtupleView->pt_2) );
        }
        else{
          if( this->Baseline("OS",cat) )      this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->Baseline("OS",cat) )      this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
          if( this->OS_W(cat) )               this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->OS_W(cat) )               this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
          if( this->relaxed_W(cat, "low") )   this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->relaxed_W(cat, "low") )   this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
          if( this->relaxed_W(cat, "high") )  this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeUp+sub,strVar)->Fill(usedVar, weight);
          if( this->relaxed_W(cat, "high") )  this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeDown+sub,strVar)->Fill(usedVar, weight);
        }

      }
      
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_WjecUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("SR_MC"+s_WjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_WjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_WjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
      
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_WjecUp+sub, strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_WjecDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("SR_MC"+s_WjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_WjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_WjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
      
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_"+s_WjecDown+sub, strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::VVSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;
    
    if(fname == s_VV){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_VV+sub,strVar)->Fill(usedVar, weight);
        if( this->LSelection() || this->TSelection() ){ //// for VV TSelection and L Selection are combined
          this->GetHistbyName(s_VVT+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_VVJ+sub,strVar)->Fill(usedVar, weight);
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_VVJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_VVJ+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_"+s_VV+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_"+s_VV+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_"+s_VV+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_"+s_VV+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_VVtauUp){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName(s_VV+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          &&( this->LSelection()
              || this->TSelection() ) )              this->GetHistbyName(s_VVT+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_VVtauDown){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName(s_VV+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          &&( this->LSelection()
              || this->TSelection() ) )              this->GetHistbyName(s_VVT+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_VVjecUp){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_VVjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->LSelection() || this->TSelection() ){
          this->GetHistbyName(s_VVTjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_VVJjecUp+sub,strVar)->Fill(usedVar, weight);
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_VVJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_VVJjecUp+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_"+s_VVjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_"+s_VVjecUp+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_"+s_VVjecUp+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_"+s_VVjecUp+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == s_VVjecDown){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName(s_VVjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->LSelection() || this->TSelection() ){
          this->GetHistbyName(s_VVTjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        else if( this->JSelection() ){
          this->GetHistbyName(s_VVJjecDown+sub,strVar)->Fill(usedVar, weight);
        }
        if( this->FFRest() ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName(s_VVJjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName(s_VVJjecDown+"_"+s_rest+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_"+s_VVjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_"+s_VVjecDown+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_"+s_VVjecDown+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_"+s_VVjecDown+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::signalSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;
    
    if(fname == s_ggH
       || fname == s_qqH){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname+"125"+sub,strVar)->Fill(usedVar, weight);
      if(fname == s_ggH){
        if( this->Baseline("OS",cat) )                 this->GetHistbyName(fname+"125"+s_CMSscalegg+s_13TeVUp+sub,strVar)->Fill(usedVar, weight*this->getRenormScale(cat) );
        if( this->Baseline("OS",cat) )                 this->GetHistbyName(fname+"125"+s_CMSscalegg+s_13TeVDown+sub,strVar)->Fill(usedVar, weight*(2-this->getRenormScale(cat) ) );
      }
    }
    else if(fname == s_ggHjecUp
       || fname == s_qqHjecUp){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname.ReplaceAll(s_jecUp,"125"+s_jecUp)+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_ggHjecDown
       || fname == s_qqHjecDown){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname.ReplaceAll(s_jecDown,"125"+s_jecDown)+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_ggHtauUp){
      if( this->Baseline("OS",cat) )                   this->GetHistbyName(s_ggH+"125"+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_ggHtauDown){
      if( this->Baseline("OS",cat) )                   this->GetHistbyName(s_ggH+"125"+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_qqHtauUp){
      if( this->Baseline("OS",cat) )                   this->GetHistbyName(s_qqH+"125"+s_CMStauScale+channel+"_"+s_13TeVUp+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == s_qqHtauDown){
      if( this->Baseline("OS",cat) )                   this->GetHistbyName(s_qqH+"125"+s_CMStauScale+channel+"_"+s_13TeVDown+sub,strVar)->Fill(usedVar, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
 
    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if( this->Baseline("OS",cat) )          this->GetHistbyName("data_obs"+sub,strVar)->Fill(usedVar, weight);

    if( channel != "tt" ){
      if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_"+s_data+sub,strVar)->Fill(usedVar, weight);
      if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_"+s_data+sub,strVar)->Fill(usedVar, weight);
      if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_"+s_data+sub,strVar)->Fill(usedVar, weight);
    }

    if(calcFF) this->applyFF(var,weight,cat,strVar,fname,1,extend);
    
}

int SelectionAnalyzer::LSelection(){

  if( ( channel != "tt"
        && NtupleView->gen_match_2 < 5 ) 
      || ( channel == "tt" &&
           !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
           && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
      ) return 1;
  else return 0;
  
}

int SelectionAnalyzer::TSelection(){

  if( ( channel != "tt"
        && NtupleView->gen_match_2 == 5 )
      || ( channel == "tt"
           && NtupleView->gen_match_1 == 5
           && NtupleView->gen_match_2 == 5 )
      ) return 1;
  else return 0;
  
}

int SelectionAnalyzer::JSelection(){

  if( ( channel != "tt"
        && NtupleView->gen_match_2 == 6 )
      || ( channel == "tt"
           && ( NtupleView->gen_match_1 == 6
                || NtupleView->gen_match_2 == 6 ) )
      ) return 1;
  else return 0;
  
}

int SelectionAnalyzer::FFRest(){

  if( calcFF
      && channel == "tt"
      && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ) return 1;
  else return 0;

}

double SelectionAnalyzer::getJetToTauFakeUp( Float_t inputPt ){

  if(inputPt<200) return 1-(0.2*inputPt/100);
  else return 0.6;
  
}

double SelectionAnalyzer::getJetToTauFakeDown( Float_t inputPt ){

  if(inputPt<200) return 1+(0.2*inputPt/100);
  else return 1.4;
  
}


