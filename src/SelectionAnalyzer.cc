#include "interface/SelectionAnalyzer.h"

using namespace std;

SelectionAnalyzer::SelectionAnalyzer(){  
}

SelectionAnalyzer::~SelectionAnalyzer(){
}

void SelectionAnalyzer::DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == "Z"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("Z"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
           ){
          this->GetHistbyName("ZLL"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("ZL"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("ZL_CMS_htt_dyShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName("ZL_CMS_htt_dyShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );          
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("ZTT"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("ZTT_CMS_htt_dyShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName("ZTT_CMS_htt_dyShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("ZLL"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("ZJ"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("ZJ_CMS_htt_dyShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight * NtupleView->ZWeight );
          this->GetHistbyName("ZJ_CMS_htt_dyShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->ZWeight );
        }
        ////////////////////////////////////////////////////////////////
        if( NtupleView->NUP == 0)          this->GetHistbyName("Z_0Jets"+sub,strVar)->Fill(usedVar, weight);
        else if( NtupleView->NUP == 1)     this->GetHistbyName("Z_1Jets"+sub,strVar)->Fill(usedVar, weight);  
        else if(  NtupleView->NUP > 1)     this->GetHistbyName("Z_ge2Jets"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("ZJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("ZJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_Z"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_Z"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_Z"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == "ZtauUp"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )       this->GetHistbyName("ZTT_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ZtauDown"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )       this->GetHistbyName("ZTT_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ZjecUp"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("ZjecUp"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
           ){
          this->GetHistbyName("ZLjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("ZTTjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("ZJjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("ZJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("ZJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_ZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_ZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_ZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_ZjecUp"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == "ZjecDown"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("ZjecDown"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
           ){
          this->GetHistbyName("ZLjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("ZTTjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("ZJjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("ZJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("ZJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_ZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_ZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_ZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_ZjecDown"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == "EWKZ"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("EWKZ"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("EWKZJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("EWKZJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_EWKZ"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_EWKZ"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_EWKZ"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_EWKZ"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == "EWKZtauUp"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )       this->GetHistbyName("EWKZ_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "EWKZtauDown"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )       this->GetHistbyName("EWKZ_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "EWKZjecUp"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("EWKZjecUp"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("EWKZJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("EWKZJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_EWKZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_EWKZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_EWKZjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_EWKZjecUp"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
    else if(fname == "EWKZjecDown"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("EWKZjecDown"+sub,strVar)->Fill(usedVar, weight);
        ////////////////////////////////////////////////////////////////
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("EWKZJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("EWKZJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);

      if( channel != "tt" ){
        if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_EWKZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_EWKZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_EWKZjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_EWKZjecDown"+sub,strVar)->Fill(usedVar, weight);
      }
      
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
    TString sub = extend+ "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if(fname == "TT"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("TT"+sub,strVar)->Fill(usedVar, weight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight * NtupleView->topWeight * NtupleView->topWeight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("TTT"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("TTT_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName("TTT_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("TTT"+sub,strVar)->Fill(usedVar, weight); 
          this->GetHistbyName("TTT_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName("TTT_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("TTJ"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("TTJ_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName("TTJ_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("TTJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("TTJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_TT"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_TT"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_TT"+sub,strVar)->Fill(usedVar, weight);

        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "TTtauUp"){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName("TT_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                          this->GetHistbyName("TTT_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "TTtauDown"){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName("TT_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                          this->GetHistbyName("TTT_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "TTjecUp"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("TTjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("TTLjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("TTTjecUp"+sub,strVar)->Fill(usedVar, weight); 
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("TTJjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("TTJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("TTJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_TTjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_TTjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_TTjecUp"+sub,strVar)->Fill(usedVar, weight);

        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_TTjecUp"+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "TTjecDown"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("TTjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("TTLjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("TTTjecDown"+sub,strVar)->Fill(usedVar, weight); 
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("TTJjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("TTJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("TTJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_TTjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_TTjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_TTjecDown"+sub,strVar)->Fill(usedVar, weight);

        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_TTjecDown"+sub,strVar)->Fill(usedVar, weight);
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
    
    if(fname == "W"){
      if( this->Baseline("OS",cat) ){ 
        //this->GetHistbyName("W_MC"+sub,strVar)->Fill(usedVar, weight);
        
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("W_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("W_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_W" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_W" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_W" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_W" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_W" + sub, strVar)->Fill(usedVar, weight);
      
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_W" + sub, strVar)->Fill(usedVar, weight);
      }
      
      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "WjecUp"){
      if( this->Baseline("OS",cat) ){         
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("WjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("WjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_WjecUp" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_WjecUp" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_WjecUp" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_WjecUp" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_WjecUp" + sub, strVar)->Fill(usedVar, weight);
      
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_WjecUp" + sub, strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "WjecDown"){
      if( this->Baseline("OS",cat) ){         
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("WjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("WjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_WjecDown" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_WjecDown" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_WjecDown" + sub, strVar)->Fill(usedVar, weight);
        if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_WjecDown" + sub, strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_WjecDown" + sub, strVar)->Fill(usedVar, weight);
      
        if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_WjecDown" + sub, strVar)->Fill(usedVar, weight);
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
    
    if(fname == "VV"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("VV"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("VVT"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("VVT"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("VVJ"+sub,strVar)->Fill(usedVar, weight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("VVJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("VVJ_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_VV"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_VV"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_VV"+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "VVtauUp"){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName("VV_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                          this->GetHistbyName("VVT_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "VVtauDown"){
      if( this->Baseline("OS",cat) )                 this->GetHistbyName("VV_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                          this->GetHistbyName("VVT_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "VVjecUp"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("VVjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("VVLjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("VVTjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("VVJjecUp"+sub,strVar)->Fill(usedVar, weight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("VVJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("VVJjecUp_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_VVjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_VVjecUp"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_VVjecUp"+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_VVjecUp"+sub,strVar)->Fill(usedVar, weight);
      }

      if(calcFF) this->applyFF(usedVar,weight,cat,strVar,fname,0,extend);
      
    }
    else if(fname == "VVjecDown"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("VVjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( ( channel != "tt" && NtupleView->gen_match_2 < 5 ) 
            || ( channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) && NtupleView->gen_match_1 < 6 && NtupleView->gen_match_2 < 6  )
            ){
          this->GetHistbyName("VVLjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
                 || ( channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
                 ){
          this->GetHistbyName("VVTjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        else if( ( channel != "tt" && NtupleView->gen_match_2 == 6 )
                 || ( channel == "tt" && ( NtupleView->gen_match_1 == 6 || NtupleView->gen_match_2 == 6 ) )
                 ){
          this->GetHistbyName("VVJjecDown"+sub,strVar)->Fill(usedVar, weight);
        }
        if( calcFF  && channel == "tt" && !( NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 ) ){
          if( NtupleView->gen_match_1 < 6 ) this->GetHistbyName("VVJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
          if( NtupleView->gen_match_2 < 6 ) this->GetHistbyName("VVJjecDown_rest"+sub,strVar)->Fill(usedVar, weight*0.5);
        }
      }

      if( channel != "tt" ){
        if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_VVjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_VVjecDown"+sub,strVar)->Fill(usedVar, weight);
        if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_VVjecDown"+sub,strVar)->Fill(usedVar, weight);
        
        if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_VVjecDown"+sub,strVar)->Fill(usedVar, weight);
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
    
    if(fname == "ggH"
       || fname == "qqH"){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname+"125"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ggHjecUp"
       || fname == "qqHjecUp"){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname.ReplaceAll("jecUp","125jecUp")+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ggHjecDown"
       || fname == "qqHjecDown"){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname.ReplaceAll("jecDown","125jecDown")+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ggHtauUp"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                            this->GetHistbyName("ggH125_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "ggHtauDown"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                            this->GetHistbyName("ggH125_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "qqHtauUp"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                            this->GetHistbyName("qqH125_CMS_scale_t_"+channel+"_13TeVUp"+sub,strVar)->Fill(usedVar, weight);
    }
    else if(fname == "qqHtauDown"){
      if( this->Baseline("OS",cat)
          && ( ( channel != "tt" && NtupleView->gen_match_2 == 5 )
               || (channel == "tt" && NtupleView->gen_match_1 == 5 && NtupleView->gen_match_2 == 5 )
               )
          )                                            this->GetHistbyName("qqH125_CMS_scale_t_"+channel+"_13TeVDown"+sub,strVar)->Fill(usedVar, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SelectionAnalyzer::dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
 
    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

    if( this->Baseline("OS",cat) )        this->GetHistbyName("data_obs"+sub,strVar)->Fill(usedVar, weight);

    if( channel != "tt" ){
      if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_data"+sub,strVar)->Fill(usedVar, weight);
      if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_data"+sub,strVar)->Fill(usedVar, weight);
      if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_data"+sub,strVar)->Fill(usedVar, weight);
    }

    if(calcFF) this->applyFF(var,weight,cat,strVar,fname,1,extend);
    
}
