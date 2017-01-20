#include "interface/CreateHistos.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TObject.h>
#include <algorithm>

using namespace std;

CreateHistos::CreateHistos(){

  histos.clear();
  files.clear();
  histo_names.clear();
  
  files.push_back({Parameter.dataset.Z,"Z"});
  files.push_back({Parameter.dataset.EWKZ,"EWKZ"});
  files.push_back({Parameter.dataset.W,"W"});
  files.push_back({Parameter.dataset.TT,"TT"});
  files.push_back({Parameter.dataset.VV,"VV"});
  if(channel=="mt")files.push_back({Parameter.dataset.data_mt,"data"});
  if(channel=="et")files.push_back({Parameter.dataset.data_et,"data"});
  if(channel=="tt")files.push_back({Parameter.dataset.data_tt,"data"});
  files.push_back({Parameter.dataset.ggH,"ggH"});
  files.push_back({Parameter.dataset.qqH,"qqH"});
  if(ptShift){
    files.push_back({Parameter.dataset.ZtauUp,"ZtauUp"});
    files.push_back({Parameter.dataset.ZtauDown,"ZtauDown"});
    files.push_back({Parameter.dataset.EWKZtauUp,"EWKZtauUp"});
    files.push_back({Parameter.dataset.EWKZtauDown,"EWKZtauDown"});
    files.push_back({Parameter.dataset.TTtauUp,"TTtauUp"});
    files.push_back({Parameter.dataset.TTtauDown,"TTtauDown"});
    files.push_back({Parameter.dataset.VVtauUp,"VVtauUp"});
    files.push_back({Parameter.dataset.VVtauDown,"VVtauDown"});
    files.push_back({Parameter.dataset.ggHtauUp,"ggHtauUp"});
    files.push_back({Parameter.dataset.ggHtauDown,"ggHtauDown"});
    files.push_back({Parameter.dataset.qqHtauUp,"qqHtauUp"});
    files.push_back({Parameter.dataset.qqHtauDown,"qqHtauDown"});
  }
  if(jecShift){
    files.push_back({Parameter.dataset.Z,"ZjecUp"});
    files.push_back({Parameter.dataset.Z,"ZjecDown"});
    files.push_back({Parameter.dataset.EWKZ,"EWKZjecUp"});
    files.push_back({Parameter.dataset.EWKZ,"EWKZjecDown"});
    files.push_back({Parameter.dataset.W,"WjecUp"});
    files.push_back({Parameter.dataset.W,"WjecDown"});
    files.push_back({Parameter.dataset.TT,"TTjecUp"});
    files.push_back({Parameter.dataset.TT,"TTjecDown"});
    files.push_back({Parameter.dataset.VV,"VVjecUp"});
    files.push_back({Parameter.dataset.VV,"VVjecDown"});
    files.push_back({Parameter.dataset.ggH,"ggHjecUp"});
    files.push_back({Parameter.dataset.ggH,"ggHjecDown"});
    files.push_back({Parameter.dataset.qqH,"qqHjecUp"});
    files.push_back({Parameter.dataset.qqH,"qqHjecDown"});
  }

  for(int i=0; i<variables.size(); i++)      vars.push_back(variables.at(i));
  for(int i=0; i<categories.size(); i++)     cats.push_back(categories.at(i));
  
  
}

CreateHistos::~CreateHistos(){

  histos.clear();

  cout << "Deleting instance of CreateHistos" << endl;

  FFinputs.clear();

  map<TString, FakeFactor*>::iterator itr;
  for (itr=FFObj.begin(); itr!=FFObj.end(); ++itr)
    {
      delete itr->second;
      FFObj.erase(itr);
    }
  map<TString, TFile*>::iterator fitr;
  for (fitr=FFfile.begin(); fitr!=FFfile.end(); ++fitr)
    {
      fitr->second->Close();
      FFfile.erase(fitr);
    }
}

void CreateHistos::loadFile(TString filename){

  TChain *tchain = new TChain("TauCheck");
  tchain->Add(filename);
  
  NtupleView = new ntuple(tchain);
  cout<<"File: "<<filename<<" loaded"<<endl;
 
}

int CreateHistos::is1DCategories(TString category){
  for( auto cat : Parameter.category.D1categories ){
    if(category == cat) return 1;
  }
  return 0;
  
}

int CreateHistos::is2DCategories(TString category){
  for( auto cat : Parameter.category.D2categories ){
    if(category == cat) return 1;
  }
  return 0;
  
}

void CreateHistos::initFakeFactors(){
  for(auto cat : cats){
    if( !this->is1DCategories(cat) && !this->is2DCategories(cat) ) continue;
    FFfile[cat] = TFile::Open("HTTutilities/Jet2TauFakes/data/"+channel+"/"+cat+"/"+FFversion);
    FFObj[cat] = (FakeFactor*)FFfile[cat]->Get("ff_comb");
  }
  FFsyst["mt"] = Parameter.FFsystematics.mt.syst;
  FFsyst["et"] = Parameter.FFsystematics.et.syst;
  FFsyst["tt"] = Parameter.FFsystematics.tt.syst;
}

void CreateHistos::run(TString isTest){


  //clearHistos();
  float weight = 1;
  float var = -999;
  
  initFakeFactors();

  cout << "Channel: " << channel << endl;
  for(auto strVar : vars) cout << "Variable " << strVar << endl;
  for(auto cat : cats) cout << "Category " << cat << endl;
  cout << endl;
  cout << "----Settings:-----" << endl;
  cout << "2D: " << do2DFit << endl;
  cout << "applyMTcut: " << applyMTCut << endl;
  cout << "useMVAMET: " << useMVAMET << endl;
  cout << "calcFF: " << calcFF << endl;
  cout << "FF version: " << FFversion << endl;
  cout << doSvfit << endl;
  cout << "Reduced string: " << reduced << endl;
  cout << endl;
  
  
  for(int i =0;i < files.size();i++){

    if(files[i][1].Contains("jecUp")) this->isJEC=1;
    else if(files[i][1].Contains("jecDown")) this->isJEC=-1;
    else this->isJEC=0;

    this->loadFile(files[i][0]);
    Int_t nentries=0;
    if(isTest=="test"){
      nentries = min( Int_t(NtupleView->fChain->GetEntries()), Int_t( 10000 ) );
    }else{
      nentries = Int_t(NtupleView->fChain->GetEntries());
    }
    cout<<"The input chain contains: "<<nentries<<" events."<<endl;
    float perc;
    for (Int_t jentry=0; jentry<nentries;jentry++){       

      if(jentry % 200000 == 0){
        if(nentries > 0){
          perc =  ( (float)jentry / (float)nentries ) * 100;
        }
        cout.precision(3);
        cout<< "                                                             \r"<< flush;
        cout<< jentry << "/" << nentries <<"\t\t" << perc << "%\r"<< flush;
      }

      NtupleView->GetEntry(jentry);    

      weight = NtupleView->stitchedWeight*NtupleView->puweight*NtupleView->effweight*NtupleView->genweight*NtupleView->antilep_tauscaling*usedLuminosity;
      if(files[i][1] == "Z"
         || files[i][1] == "ZtauUp"
         || files[i][1] == "ZtauDown"
         || files[i][1] == "ZjecUp"
         || files[i][1] == "ZjecDown"
         || files[i][1] == "EWKZ"
         || files[i][1] == "EWKZjecUp"
         || files[i][1] == "EWKZjecDown"
         || files[i][1] == "EWKZtauUp"
         || files[i][1] == "EWKZtauDown") weight *= NtupleView->ZWeight;
      if(files[i][1] == "TT"
         || files[i][1] == "TTtauUp"
         || files[i][1] == "TTtauDown"
         || files[i][1] == "TTjecUp"
         || files[i][1] == "TTjecDown") weight *= NtupleView->topWeight;

      //if(NtupleView->idisoweight_2 != 1) weight = weight * (0.9/0.83);
      if(channel == "et") weight = weight * this->getAntiLep_tauscaling();

      for(auto cat : cats){

        for(auto strVar : vars){

          var = -999;

          if(strVar == "m_vis")                              var = NtupleView->m_vis;
          else if(strVar == "m_sv")                          var = NtupleView->m_sv;
          else if(strVar == "pt_sv")                         var = NtupleView->pt_sv;
          else if(strVar == "mt_1")                          var = this->getMT();
          else if(strVar == "jpt_1")                         var = NtupleView->jpt_1;
          else if(strVar == "jpt_2")                         var = NtupleView->jpt_2;

          else if(strVar == "pt_1")                          var = NtupleView->pt_1;
          else if(strVar == "pt_2")                          var = NtupleView->pt_2;
          else if(strVar == "eta_1")                         var = NtupleView->eta_1;
          else if(strVar == "eta_2")                         var = NtupleView->eta_2;
          else if(strVar == "met")                           var = NtupleView->met;
          else if(strVar == "mttot")                         var = this->getMTTOT();
          else if(strVar == "Hpt")                           var = this->CalcHPt();
          else if(strVar == "pfmt_1")                        var = NtupleView->pfmt_1;


          else if(strVar == "jpt_1_2"
                  && fabs(NtupleView->jeta_1) < 2.5)         var = NtupleView->jpt_1;
          else if(strVar == "jpt_1_2p5"
                  && fabs(NtupleView->jeta_1) >= 2.5
                  && fabs(NtupleView->jeta_1) < 3.)          var = NtupleView->jpt_1;
          else if(strVar == "jpt_1_3"
                  && fabs(NtupleView->jeta_1) >= 3.)         var = NtupleView->jpt_1;



          else if(strVar == "jpt_2_2"
                  && fabs(NtupleView->jeta_1) < 2.5)         var = NtupleView->jpt_2;
          else if(strVar == "jpt_2_2p5"
                  && fabs(NtupleView->jeta_1) >= 2.5
                  && fabs(NtupleView->jeta_1) < 3.)          var = NtupleView->jpt_2;
          else if(strVar == "jpt_2_3"
                  && fabs(NtupleView->jeta_1) >= 3.)         var = NtupleView->jpt_2;



          else if(strVar == "jeta_1"
             && NtupleView->jpt_1 > 30)     var = NtupleView->jeta_1;

          else if(strVar == "jeta_2"
             && NtupleView->jpt_2 > 30)     var = NtupleView->jeta_2;

          /*else if(strVar == "mjj"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = NtupleView->mjj;*/

          else if(strVar == "mjj")          var = this->getMjj();

          else if(strVar == "jeta1eta2"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = NtupleView->jeta1eta2;

          else if(strVar == "jdeta"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = this->getJdeta();

          else continue;


          if(files[i][1] == "Z" 
             || files[i][1] == "ZtauUp"
             || files[i][1] == "ZtauDown"
             || files[i][1] == "ZjecUp"
             || files[i][1] == "ZjecDown")               this->DYSelections(var, weight, cat, strVar, files[i][1]);
          
          else if(files[i][1] == "EWKZ" 
                  || files[i][1] == "EWKZtauUp"
                  || files[i][1] == "EWKZtauDown"
                  || files[i][1] == "EWKZjecUp"
                  || files[i][1] == "EWKZjecDown")       this->EWKZSelections(var, weight, cat, strVar, files[i][1]);
          
          else if(files[i][1] == "TT"
                  || files[i][1] == "TTtauUp"
                  || files[i][1] == "TTtauDown"
                  || files[i][1] == "TTjecUp"
                  || files[i][1] == "TTjecDown")         this->TSelections(var, weight, cat, strVar, files[i][1]);
          
          else if(files[i][1] == "VV"
                  || files[i][1] == "VVtauUp"
                  || files[i][1] == "VVtauDown"
                  || files[i][1] == "VVjecUp"
                  || files[i][1] == "VVjecDown")         this->VVSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "W"
                  || files[i][1] == "WjecUp"
                  || files[i][1] == "WjecDown")          this->WSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "data" )                this->dataSelections(var, 1., cat, strVar, files[i][1]);

          else if(files[i][1] == "qqH"
                  || files[i][1] == "ggH"
                  || files[i][1] == "qqHtauUp"
                  || files[i][1] == "qqHtauDown"
                  || files[i][1] == "qqHjecUp"
                  || files[i][1] == "qqHjecDown"
                  || files[i][1] == "ggHtauUp"
                  || files[i][1] == "ggHtauDown"
                  || files[i][1] == "ggHjecUp"
                  || files[i][1] == "ggHjecDown")        this->signalSelections(var, weight, cat, strVar, files[i][1]);
          

          if(do2DFit){

            if ( !this->is2DCategories(cat) ) continue;
            
            if(files[i][1] == "Z" 
               || files[i][1] == "ZtauUp"
               || files[i][1] == "ZtauDown"
               || files[i][1] == "ZjecUp"
               || files[i][1] == "ZjecDown")             this->DYSelections(var, weight, cat, strVar, files[i][1], "2D");

            else if(files[i][1] == "EWKZ" 
                    || files[i][1] == "EWKZtauUp"
                    || files[i][1] == "EWKZtauDown"
                    || files[i][1] == "EWKZjecUp"
                    || files[i][1] == "EWKZjecDown")     this->EWKZSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if(files[i][1] == "TT"
                    || files[i][1] == "TTtauUp"
                    || files[i][1] == "TTtauDown"
                    || files[i][1] == "TTjecUp"
                    || files[i][1] == "TTjecDown")       this->TSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if(files[i][1] == "VV"
                    || files[i][1] == "VVtauUp"
                    || files[i][1] == "VVtauDown"
                    || files[i][1] == "VVjecUp"
                    || files[i][1] == "VVjecDown")       this->VVSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "W"
                    || files[i][1] == "WjecUp"
                    || files[i][1] == "WjecDown")         this->WSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "data" )               this->dataSelections(var, 1., cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "qqH"
                    || files[i][1] == "ggH"
                    || files[i][1] == "qqHtauUp"
                    || files[i][1] == "qqHtauDown"
                    || files[i][1] == "qqHjecUp"
                    || files[i][1] == "qqHjecDown"
                    || files[i][1] == "ggHtauUp"
                    || files[i][1] == "ggHtauDown"
                    || files[i][1] == "ggHjecUp"
                    || files[i][1] == "ggHjecDown")       this->signalSelections(var, weight, cat, strVar, files[i][1], "2D");
            
          }
        }
        
      }
    }
  }
  if( channel != "tt" ){
    for(auto cat : cats){
      for(auto strVar : vars){
        this->Estimate_W_QCD(strVar, cat);
        if( do2DFit && this->is2DCategories(cat) ) this->Estimate_W_QCD(strVar, cat, "2D"); 
      }
    }
  }
  if(calcFF){
    for(auto cat : cats){
      for(auto strVar : vars){
        this->EstimateFF(strVar, cat);
        if( do2DFit && this->is2DCategories(cat) ) this->EstimateFF(strVar, cat, "2D");
      }
    } 
  }
  
  cout << "Done running over events." << endl;
  writeHistos( channel, cats, vars );
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float CreateHistos::getAntiLep_tauscaling(){


    if(channel == "mt"){
      if(NtupleView->againstMuonTight3_2 > 0.5
         && ( NtupleView->gen_match_2 == 2 
             || NtupleView->gen_match_2 == 4)
         ){

        if(fabs(NtupleView->eta_2) < 1.2) return 1.28;       // +-0.06
        else if (fabs(NtupleView->eta_2) < 1.7) return 2.6;   // +-2.6
        else if (fabs(NtupleView->eta_2) < 2.3) return 2.1;  // +-0.9
      }
    }
    if(channel == "et"){
      if(NtupleView->againstElectronTightMVA6_2 > 0.5
         && ( NtupleView->gen_match_2 == 1 
             || NtupleView->gen_match_2 == 3)
         ){

        if(fabs(NtupleView->eta_2) < 1.46) return 1.42;       // +-0.06
        else if (fabs(NtupleView->eta_2) > 1.558) return 1.994;   // +-2.6
      }
    }
    return 1.0;

}

void CreateHistos::getFFInputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( this->getMT() );
  inputs.push_back( NtupleView->iso_1 );
}
void CreateHistos::getFF1Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_1 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}
void CreateHistos::getFF2Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( this->getNjets() );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}


void CreateHistos::applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend){

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

int CreateHistos::returnBin(vector<double> bins, double value){

  if( value<bins.at(0) ) return -1;
  else if( value > bins.at( bins.size() -1 ) ) return -1;
  for(int i=0; i<bins.size()-1; i++){
    if( value>bins.at(i) && value<bins.at(i+1) ) return i;
  }
  return -1;
}

double CreateHistos::get2DVar(TString sub){

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

void CreateHistos::DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

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
void CreateHistos::EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

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
void CreateHistos::TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
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
void CreateHistos::WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

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
void CreateHistos::VVSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

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
void CreateHistos::signalSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

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
void CreateHistos::dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){
 
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD_osw(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";

  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_data"+sub,strVar)   );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_VV"+sub,strVar), -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_Z"+sub,strVar),  -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_TT"+sub,strVar), -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_W"+sub,strVar),  -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_EWKZ"+sub,strVar),  -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Scale( this->QCD_OSSS(cat) );

  if(resetZero){
    int entries = this->GetHistbyName("QCD_OSW"+sub,strVar)->GetNbinsX();
    for(int i = 1; i <= entries; i++){
      if( this->GetHistbyName("QCD_OSW"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCD_OSW"+sub,strVar)->SetBinContent(i,0.);
    }
  }

  if(jecShift){
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_data"+sub,strVar)   );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_VVjecUp"+sub,strVar), -1 );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_ZjecUp"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_TTjecUp"+sub,strVar), -1 );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_WjecUp"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_EWKZjecUp"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->Scale( this->QCD_OSSS(cat) );
   
   if(resetZero){
     int entries = this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->GetNbinsX();
     for(int i = 1; i <= entries; i++){
       if( this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCDjecUp_OSW"+sub,strVar)->SetBinContent(i,0.);
     }
   }
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_data"+sub,strVar)   );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_VVjecDown"+sub,strVar), -1 );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_Z"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_TT"+sub,strVar), -1 );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_W"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_EWKZ"+sub,strVar),  -1 );
   this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->Scale( this->QCD_OSSS(cat) );
   
   if(resetZero){
     int entries = this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->GetNbinsX();
     for(int i = 1; i <= entries; i++){
       if( this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCDjecDown_OSW"+sub,strVar)->SetBinContent(i,0.);
     }
   }
   
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateW(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";

  //OLD Wjets estimation ala AN-16-355
  /*this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_data"+sub,strVar)   );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_VV"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_Z"+sub,strVar),  -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_TT"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("QCD_OSW"+sub,strVar), -1 );

  this->GetHistbyName("HL"+sub, strVar)->Add( this->GetHistbyName("relaxed_W_low_W"+sub,strVar) );
  this->GetHistbyName("HL"+sub, strVar)->Divide( this->GetHistbyName("relaxed_W_high_W"+sub,strVar) );

  this->GetHistbyName("W"+sub,strVar)->Add( this->GetHistbyName("W_OSW"+sub,strVar) );
  this->GetHistbyName("W"+sub,strVar)->Multiply( this->GetHistbyName("HL"+sub, strVar) );*/

  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_data"+sub,strVar)   );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_VV"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_Z"+sub,strVar),  -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_TT"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_EWKZ"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("QCD_OSW"+sub,strVar), -1 );

  double w_normFactor = this->GetHistbyName("W_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_W"+sub,strVar)->Integral();

  this->GetHistbyName("W"+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_W"+sub,strVar) );
  this->GetHistbyName("W"+sub,strVar)->Scale( w_normFactor );
  

  if(resetZero){
    int entries = this->GetHistbyName("W"+sub,strVar)->GetNbinsX();
    for(int i = 1; i <= entries; i++){
      if( this->GetHistbyName("W"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("W"+sub,strVar)->SetBinContent(i,0.);
    }
  }

  if(jecShift){
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_data"+sub,strVar)   );
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_VVjecUp"+sub,strVar), -1 );
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_ZjecUp"+sub,strVar),  -1 );
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_TTjecUp"+sub,strVar), -1 );
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_EWKZjecUp"+sub,strVar), -1 );
    this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Add( this->GetHistbyName("QCDjecUp_OSW"+sub,strVar), -1 );
    
    w_normFactor = this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_WjecUp"+sub,strVar)->Integral();
    
    this->GetHistbyName("WjecUp"+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_WjecUp"+sub,strVar) );
    this->GetHistbyName("WjecUp"+sub,strVar)->Scale( w_normFactor );
    
    
    if(resetZero){
      int entries = this->GetHistbyName("WjecUp"+sub,strVar)->GetNbinsX();
      for(int i = 1; i <= entries; i++){
        if( this->GetHistbyName("WjecUp"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("WjecUp"+sub,strVar)->SetBinContent(i,0.);
      }
    }
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_data"+sub,strVar)   );
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_VVjecDown"+sub,strVar), -1 );
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_ZjecDown"+sub,strVar),  -1 );
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_TTjecDown"+sub,strVar), -1 );
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_EWKZjecDown"+sub,strVar), -1 );
    this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Add( this->GetHistbyName("QCDjecDown_OSW"+sub,strVar), -1 );
    
    w_normFactor = this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_WjecDown"+sub,strVar)->Integral();
    
    this->GetHistbyName("WjecDown"+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_WjecDown"+sub,strVar) );
    this->GetHistbyName("WjecDown"+sub,strVar)->Scale( w_normFactor );
    
    
    if(resetZero){
      int entries = this->GetHistbyName("WjecDown"+sub,strVar)->GetNbinsX();
      for(int i = 1; i <= entries; i++){
        if( this->GetHistbyName("WjecDown"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("WjecDown"+sub,strVar)->SetBinContent(i,0.);
      }
    } 
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";
  
  double W_SF = ( this->GetHistbyName("W_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("W_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_W"+sub,strVar)->Integral() : 0;
  double W_rel = ( this->GetHistbyName("SS_Low_W"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_W"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Integral() : 0;
  double Z_rel = ( this->GetHistbyName("SS_Low_Z"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_Z"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Integral() : 0;
  double VV_rel = ( this->GetHistbyName("SS_Low_VV"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_VV"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Integral() : 0;
  double TT_rel = ( this->GetHistbyName("SS_Low_TT"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_TT"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Integral() : 0;
  double EWKZ_rel = ( this->GetHistbyName("SS_Low_EWKZ"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_EWKZ"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_EWKZ"+sub,strVar)->Integral() : 0;

  

  this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Scale(  W_SF   );
  this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Scale(  W_rel  );
  this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Scale( VV_rel );
  this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Scale(  Z_rel  );
  this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Scale( TT_rel );
  this->GetHistbyName("SS_Low_relaxed_EWKZ"+sub,strVar)->Scale( EWKZ_rel );

  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_data"+sub,strVar)  );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar), -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar),  -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar), -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar),  -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_EWKZ"+sub,strVar),  -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Scale( this->QCD_OSSS(cat) );

  if(resetZero){
    int entries = this->GetHistbyName("QCD"+sub,strVar)->GetNbinsX();
    for(int i = 1; i < entries; i++){
      if( this->GetHistbyName("QCD"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCD"+sub,strVar)->SetBinContent(i,0.);
    }
  }

  if(jecShift){
    W_SF = ( this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("WjecUp_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_WjecUp"+sub,strVar)->Integral() : 0;
    W_rel = ( this->GetHistbyName("SS_Low_WjecUp"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_WjecUp"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_WjecUp"+sub,strVar)->Integral() : 0;
    Z_rel = ( this->GetHistbyName("SS_Low_ZjecUp"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_ZjecUp"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_ZjecUp"+sub,strVar)->Integral() : 0;
    VV_rel = ( this->GetHistbyName("SS_Low_VVjecUp"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_VVjecUp"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_VVjecUp"+sub,strVar)->Integral() : 0;
    TT_rel = ( this->GetHistbyName("SS_Low_TTjecUp"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_TTjecUp"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_TTjecUp"+sub,strVar)->Integral() : 0;
    EWKZ_rel = ( this->GetHistbyName("SS_Low_EWKZjecUp"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_EWKZjecUp"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_EWKZjecUp"+sub,strVar)->Integral() : 0;
    
    this->GetHistbyName("SS_Low_relaxed_WjecUp"+sub,strVar)->Scale(  W_SF   );
    this->GetHistbyName("SS_Low_relaxed_WjecUp"+sub,strVar)->Scale(  W_rel  );
    this->GetHistbyName("SS_Low_relaxed_VVjecUp"+sub,strVar)->Scale( VV_rel );
    this->GetHistbyName("SS_Low_relaxed_ZjecUp"+sub,strVar)->Scale(  Z_rel  );
    this->GetHistbyName("SS_Low_relaxed_TTjecUp"+sub,strVar)->Scale( TT_rel );
    this->GetHistbyName("SS_Low_relaxed_EWKZjecUp"+sub,strVar)->Scale( EWKZ_rel );

    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_data"+sub,strVar)  );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_VVjecUp"+sub,strVar), -1 );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_ZjecUp"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_TTjecUp"+sub,strVar), -1 );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_WjecUp"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_EWKZjecUp"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecUp"+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  
    if(resetZero){
      int entries = this->GetHistbyName("QCDjecUp"+sub,strVar)->GetNbinsX();
      for(int i = 1; i < entries; i++){
        if( this->GetHistbyName("QCDjecUp"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCDjecUp"+sub,strVar)->SetBinContent(i,0.);
      }
    }
    W_SF = ( this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("WjecDown_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_WjecDown"+sub,strVar)->Integral() : 0;
    W_rel = ( this->GetHistbyName("SS_Low_WjecDown"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_WjecDown"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_WjecDown"+sub,strVar)->Integral() : 0;
    Z_rel = ( this->GetHistbyName("SS_Low_ZjecDown"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_ZjecDown"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_ZjecDown"+sub,strVar)->Integral() : 0;
    VV_rel = ( this->GetHistbyName("SS_Low_VVjecDown"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_VVjecDown"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_VVjecDown"+sub,strVar)->Integral() : 0;
    TT_rel = ( this->GetHistbyName("SS_Low_TTjecDown"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_TTjecDown"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_TTjecDown"+sub,strVar)->Integral() : 0;
    EWKZ_rel = ( this->GetHistbyName("SS_Low_EWKZjecDown"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_EWKZjecDown"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_EWKZjecDown"+sub,strVar)->Integral() : 0;

    this->GetHistbyName("SS_Low_relaxed_WjecDown"+sub,strVar)->Scale(  W_SF   );
    this->GetHistbyName("SS_Low_relaxed_WjecDown"+sub,strVar)->Scale(  W_rel  );
    this->GetHistbyName("SS_Low_relaxed_VVjecDown"+sub,strVar)->Scale( VV_rel );
    this->GetHistbyName("SS_Low_relaxed_ZjecDown"+sub,strVar)->Scale(  Z_rel  );
    this->GetHistbyName("SS_Low_relaxed_TTjecDown"+sub,strVar)->Scale( TT_rel );
    this->GetHistbyName("SS_Low_relaxed_EWKZjecDown"+sub,strVar)->Scale( EWKZ_rel );

    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_data"+sub,strVar)  );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_VVjecDown"+sub,strVar), -1 );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_ZjecDown"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_TTjecDown"+sub,strVar), -1 );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_WjecDown"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_EWKZjecDown"+sub,strVar),  -1 );
    this->GetHistbyName("QCDjecDown"+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  
    if(resetZero){
      int entries = this->GetHistbyName("QCDjecDown"+sub,strVar)->GetNbinsX();
      for(int i = 1; i < entries; i++){
        if( this->GetHistbyName("QCDjecDown"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCDjecDown"+sub,strVar)->SetBinContent(i,0.);
      }
    }
  }

  

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::Estimate_W_QCD(TString strVar, TString cat, TString extend){
  this->CreateQCD_osw(strVar, cat, extend);
  this->CreateW(strVar, cat, extend);
  this->CreateQCD(strVar, cat, extend);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::EstimateFF(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";

  double normUp_jetFakes=0;
  double normUp_jetFakes_syst=0;
  double normUp_data=0;
  double normUp_data_syst=0;
  double normUp_W=0;
  double normUp_W_syst=0;
  double normUp_Z=0;
  double normUp_Z_syst=0;
  double normUp_EWKZ=0;
  double normUp_EWKZ_syst=0;  
  double normUp_TT=0;
  double normUp_TT_syst=0;
  double normUp_VV=0;
  double normUp_VV_syst=0;
  double normDown_jetFakes=0;
  double normDown_jetFakes_syst=0;
  double normDown_data=0;
  double normDown_data_syst=0;
  double normDown_W=0;
  double normDown_W_syst=0;
  double normDown_Z=0;
  double normDown_Z_syst=0;
  double normDown_EWKZ=0;
  double normDown_EWKZ_syst=0;
  double normDown_TT=0;
  double normDown_TT_syst=0;
  double normDown_VV=0;
  double normDown_VV_syst=0;
  
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("EWKZ_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  if(jecShift){
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("WjecUp_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("TTjecUp_jetFakes"+sub,strVar),  -1 );
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("ZjecUp_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("EWKZjecUp_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecUp"+sub,strVar)->Add( this->GetHistbyName("VVjecUp_jetFakes"+sub,strVar), -1 );

    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("WjecDown_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("TTjecDown_jetFakes"+sub,strVar),  -1 );
    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("ZjecDown_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("EWKZjecDown_jetFakes"+sub,strVar), -1 );
    this->GetHistbyName("jetFakesjecDown"+sub,strVar)->Add( this->GetHistbyName("VVjecDown_jetFakes"+sub,strVar), -1 );
  }

  this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_Up"+sub,strVar), -1 );  
  this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_Down"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_Up"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_Down"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  

  if(resetZero){
    for(int i = 1; i < this->GetHistbyName("jetFakes"+sub,strVar)->GetNbinsX(); i++){
      if( this->GetHistbyName("jetFakes"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes"+sub,strVar)->SetBinContent(i,0.);
      if(jecShift) {
        if( this->GetHistbyName("jetFakesjecUp"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakesjecUp"+sub,strVar)->SetBinContent(i,0.);
        if( this->GetHistbyName("jetFakesjecDown"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakesjecDown"+sub,strVar)->SetBinContent(i,0.);
      }
      if( this->GetHistbyName("jetFakes"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes"+sub,strVar)->SetBinContent(i,0.);
      if( this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes_Zvar_systUp"+sub,strVar)->SetBinContent(i,0.);
      if( this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes_Zvar_systDown"+sub,strVar)->SetBinContent(i,0.);
      if( this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes_TTvar_systUp"+sub,strVar)->SetBinContent(i,0.);
      if( this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes_TTvar_systDown"+sub,strVar)->SetBinContent(i,0.);
    }
  }
  
  for( auto syst : FFsyst[channel] ){
    TString tmp=syst; tmp.ReplaceAll("_down","Down"); tmp.ReplaceAll("_up","Up"); 
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("data_jetFakes_"+tmp+sub,strVar)   );
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("W_jetFakes_"+tmp+sub,strVar), -1 );
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_"+tmp+sub,strVar),  -1 );
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_"+tmp+sub,strVar), -1 );
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("EWKZ_jetFakes_"+tmp+sub,strVar), -1 );
    this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes_"+tmp+sub,strVar), -1 );
    if(resetZero){
      for(int i = 1; i < this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->GetNbinsX(); i++){
        if( this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("jetFakes_"+tmp+sub,strVar)->SetBinContent(i,0.);
      }
    }
    //FIXME: if combination of uncertainties and normalization to std histo is required this is a first attempt 
    this->GetHistbyName("jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("jetFakes_"+tmp+sub,strVar)  );
    double ratio = ( this->GetHistbyName("jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_jetFakes_syst = TMath::Sqrt( TMath::Power(normUp_jetFakes_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_jetFakes_syst = TMath::Sqrt( TMath::Power(normDown_jetFakes_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_jetFakes = TMath::Sqrt( TMath::Power(normUp_jetFakes,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_jetFakes = TMath::Sqrt( TMath::Power(normDown_jetFakes,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName("data_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("data_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("data_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("data_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("data_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("data_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("data_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_data_syst = TMath::Sqrt( TMath::Power(normUp_data_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_data_syst = TMath::Sqrt( TMath::Power(normDown_data_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_data = TMath::Sqrt( TMath::Power(normUp_data,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_data = TMath::Sqrt( TMath::Power(normDown_data,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName("W_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("W_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("W_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("W_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("W_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("W_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("W_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_W_syst = TMath::Sqrt( TMath::Power(normUp_W_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_W_syst = TMath::Sqrt( TMath::Power(normDown_W_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_W = TMath::Sqrt( TMath::Power(normUp_W,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_W = TMath::Sqrt( TMath::Power(normDown_W,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName("TT_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("TT_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("TT_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("TT_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("TT_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("TT_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_TT_syst = TMath::Sqrt( TMath::Power(normUp_TT_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_TT_syst = TMath::Sqrt( TMath::Power(normDown_TT_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_TT = TMath::Sqrt( TMath::Power(normUp_TT,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_TT = TMath::Sqrt( TMath::Power(normDown_TT,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName("Z_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("Z_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("Z_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("Z_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("Z_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("Z_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_Z_syst = TMath::Sqrt( TMath::Power(normUp_Z_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_Z_syst = TMath::Sqrt( TMath::Power(normDown_Z_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_Z = TMath::Sqrt( TMath::Power(normUp_Z,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_Z = TMath::Sqrt( TMath::Power(normDown_Z,2) + TMath::Power( 1-ratio,2 ) );
    }

    this->GetHistbyName("EWKZ_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("EWKZ_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("EWKZ_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("EWKZ_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("EWKZ_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("EWKZ_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("EWKZ_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_EWKZ_syst = TMath::Sqrt( TMath::Power(normUp_EWKZ_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_EWKZ_syst = TMath::Sqrt( TMath::Power(normDown_EWKZ_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_EWKZ = TMath::Sqrt( TMath::Power(normUp_EWKZ,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_EWKZ = TMath::Sqrt( TMath::Power(normDown_EWKZ,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName("VV_jetFakes_norm_"+tmp+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName("VV_jetFakes"+sub,strVar)->Integral( 0, this->GetHistbyName("VV_jetFakes"+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName("VV_jetFakes_norm_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName("VV_jetFakes_norm_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName("VV_jetFakes_norm_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio > 1) normUp_VV_syst = TMath::Sqrt( TMath::Power(normUp_VV_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_VV_syst = TMath::Sqrt( TMath::Power(normDown_VV_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio > 1) normUp_VV = TMath::Sqrt( TMath::Power(normUp_VV,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_VV = TMath::Sqrt( TMath::Power(normDown_VV,2) + TMath::Power( 1-ratio,2 ) );
    }
  }
  this->GetHistbyName("jetFakes_norm"+sub,"norm");
  this->GetHistbyName("data_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("W_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("TT_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("Z_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("EWKZ_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("VV_jetFakes_norm"+sub,"norm");
  this->GetHistbyName("jetFakes_norm"+sub)->SetBinContent(1,normUp_jetFakes);
  this->GetHistbyName("data_jetFakes_norm"+sub)->SetBinContent(1,normUp_data);
  this->GetHistbyName("W_jetFakes_norm"+sub)->SetBinContent(1,normUp_W);
  this->GetHistbyName("TT_jetFakes_norm"+sub)->SetBinContent(1,normUp_TT);
  this->GetHistbyName("Z_jetFakes_norm"+sub)->SetBinContent(1,normUp_Z);
  this->GetHistbyName("EWKZ_jetFakes_norm"+sub)->SetBinContent(1,normUp_EWKZ);
  this->GetHistbyName("VV_jetFakes_norm"+sub)->SetBinContent(1,normUp_VV);
  this->GetHistbyName("jetFakes_norm"+sub)->SetBinContent(2,normDown_jetFakes);
  this->GetHistbyName("data_jetFakes_norm"+sub)->SetBinContent(2,normDown_data);
  this->GetHistbyName("W_jetFakes_norm"+sub)->SetBinContent(2,normDown_W);
  this->GetHistbyName("TT_jetFakes_norm"+sub)->SetBinContent(2,normDown_TT);
  this->GetHistbyName("Z_jetFakes_norm"+sub)->SetBinContent(2,normDown_Z);
  this->GetHistbyName("EWKZ_jetFakes_norm"+sub)->SetBinContent(2,normDown_EWKZ);
  this->GetHistbyName("VV_jetFakes_norm"+sub)->SetBinContent(2,normDown_VV);
  this->GetHistbyName("jetFakes_norm"+sub)->SetBinContent(3,normUp_jetFakes_syst);
  this->GetHistbyName("data_jetFakes_norm"+sub)->SetBinContent(3,normUp_data_syst);
  this->GetHistbyName("W_jetFakes_norm"+sub)->SetBinContent(3,normUp_W_syst);
  this->GetHistbyName("TT_jetFakes_norm"+sub)->SetBinContent(3,normUp_TT_syst);
  this->GetHistbyName("Z_jetFakes_norm"+sub)->SetBinContent(3,normUp_Z_syst);
  this->GetHistbyName("EWKZ_jetFakes_norm"+sub)->SetBinContent(3,normUp_EWKZ_syst);
  this->GetHistbyName("VV_jetFakes_norm"+sub)->SetBinContent(3,normUp_VV_syst);
  this->GetHistbyName("jetFakes_norm"+sub)->SetBinContent(4,normDown_jetFakes_syst);
  this->GetHistbyName("data_jetFakes_norm"+sub)->SetBinContent(4,normDown_data_syst);
  this->GetHistbyName("W_jetFakes_norm"+sub)->SetBinContent(4,normDown_W_syst);
  this->GetHistbyName("TT_jetFakes_norm"+sub)->SetBinContent(4,normDown_TT_syst);
  this->GetHistbyName("Z_jetFakes_norm"+sub)->SetBinContent(4,normDown_Z_syst);
  this->GetHistbyName("EWKZ_jetFakes_norm"+sub)->SetBinContent(4,normDown_EWKZ_syst);
  this->GetHistbyName("VV_jetFakes_norm"+sub)->SetBinContent(4,normDown_VV_syst);

  if(channel == "tt"){
    this->GetHistbyName("FF_rest"+sub,strVar)->Add( this->GetHistbyName("ZJ_rest"+sub,strVar)   );
    this->GetHistbyName("FF_rest"+sub,strVar)->Add( this->GetHistbyName("EWKZJ_rest"+sub,strVar)   );
    this->GetHistbyName("FF_rest"+sub,strVar)->Add( this->GetHistbyName("W_rest"+sub,strVar)   );
    this->GetHistbyName("FF_rest"+sub,strVar)->Add( this->GetHistbyName("TTJ_rest"+sub,strVar)   );
    this->GetHistbyName("FF_rest"+sub,strVar)->Add( this->GetHistbyName("VVJ_rest"+sub,strVar)   );

    if(jecShift){
      this->GetHistbyName("FF_restjecUp"+sub,strVar)->Add( this->GetHistbyName("ZJ_restjecUp"+sub,strVar)   );
      this->GetHistbyName("FF_restjecUp"+sub,strVar)->Add( this->GetHistbyName("EWKZJ_restjecUp"+sub,strVar)   );
      this->GetHistbyName("FF_restjecUp"+sub,strVar)->Add( this->GetHistbyName("W_restjecUp"+sub,strVar)   );
      this->GetHistbyName("FF_restjecUp"+sub,strVar)->Add( this->GetHistbyName("TTJ_restjecUp"+sub,strVar)   );
      this->GetHistbyName("FF_restjecUp"+sub,strVar)->Add( this->GetHistbyName("VVJ_restjecUp"+sub,strVar)   );

      this->GetHistbyName("FF_restjecDown"+sub,strVar)->Add( this->GetHistbyName("ZJ_restjecDown"+sub,strVar)   );
      this->GetHistbyName("FF_restjecDown"+sub,strVar)->Add( this->GetHistbyName("EWKZJ_restjecDown"+sub,strVar)   );
      this->GetHistbyName("FF_restjecDown"+sub,strVar)->Add( this->GetHistbyName("W_restjecDown"+sub,strVar)   );
      this->GetHistbyName("FF_restjecDown"+sub,strVar)->Add( this->GetHistbyName("TTJ_restjecDown"+sub,strVar)   );
      this->GetHistbyName("FF_restjecDown"+sub,strVar)->Add( this->GetHistbyName("VVJ_restjecDown"+sub,strVar)   );
    }
  }
  //cout << "NormUp: " << normUp << endl;
  //cout << "NormDown: " << normDown << endl;
  //cout << "//////////////////////////////////////" << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CreateHistos::OS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::SS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::relaxed_W(TString cat, TString mt){

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

int CreateHistos::SS_Low(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("base")
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::SS_Low_relaxed(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("relaxed")
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat, "wo")
     && this->Vetos()) return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CreateHistos::QCD_OSSS(TString cat){
  if(channel == "mt"){
    if(cat == "0Jet"
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
    if(cat == "0Jet"
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CreateHistos::writeHistos( TString channel, vector<TString> cats, vector<TString> vars ){

  stringstream outfile_name;
  TString sub;
  TString tmp;
  TString D2="";
  if(do2DFit) D2+="2D";

  
  for(auto var : vars){
    TString suffix = (doOfficialNaming) ? "" : var+"-";
    outfile_name << "histos/"  << "htt_" << channel << ".inputs-sm-13TeV-"<<suffix<<D2<<".root";
    outfile = new TFile(outfile_name.str().c_str(), "RECREATE") ;

    for(auto cat : cats){
      outfile->mkdir(channel +"_"+ cat );
      outfile->cd(channel +"_"+ cat); 
      sub = "+" + var +"_" + cat + "+";

      for( int i=0;i<histo_names.size();i++ ){

        if(histo_names.at(i).Contains(sub) ){
          if(!keepDebugHistos
             && ( histo_names.at(i).Contains("SS_Low")
                  || histo_names.at(i).Contains("relaxed")
                  || histo_names.at(i).Contains("SS_W")
                  || histo_names.at(i).Contains("OS_W")
                  || histo_names.at(i).Contains("QCD_OSW")
                  || histo_names.at(i).Contains("W_OSW")
                  || histo_names.at(i).Contains("WjecUp_OSW")
                  || histo_names.at(i).Contains("WjecDown_OSW")
                  )
             ) continue;
          if(!keepFFDebugHistos
             && ( histo_names.at(i).Contains("data_jetFakes")
                  || histo_names.at(i).Contains("W_jetFakes")
                  || histo_names.at(i).Contains("WjecUp_jetFakes")
                  || histo_names.at(i).Contains("WjecDown_jetFakes")
                  || histo_names.at(i).Contains("Z_jetFakes")
                  || histo_names.at(i).Contains("ZjecUp_jetFakes")
                  || histo_names.at(i).Contains("ZjecDown_jetFakes")
                  || histo_names.at(i).Contains("TT_jetFakes")
                  || histo_names.at(i).Contains("TTjecUp_jetFakes")
                  || histo_names.at(i).Contains("TTjecDown_jetFakes")
                  || histo_names.at(i).Contains("VV_jetFakes")
                  || histo_names.at(i).Contains("VVjecUp_jetFakes")
                  || histo_names.at(i).Contains("VVjecDown_jetFakes")
                  || histo_names.at(i).Contains("EWKZ_jetFakes")
                  || histo_names.at(i).Contains("EWKZjecUp_jetFakes")
                  || histo_names.at(i).Contains("EWKZjecDown_jetFakes")                  
                  )
             ) continue;
          if(!keepZGenJetsSplitting
             && ( histo_names.at(i).Contains("Z_0Jets")
                  || histo_names.at(i).Contains("Z_1Jets")
                  || histo_names.at(i).Contains("Z_ge2Jets")
                  )
             ) continue;
          tmp = histo_names.at(i);
          tmp.ReplaceAll(sub, "");
          tmp.ReplaceAll("jecUp","_CMS_scale_j_13TeVUp");
          tmp.ReplaceAll("jecDown","_CMS_scale_j_13TeVDown");
          histos.at(i)->SetName(tmp);
          if(do2DFit ){
            if( is2DCategories(cat) ){
              tmp.ReplaceAll("2D","");
              histos.at(i)->Write(tmp, TObject::kWriteDelete);
            }
            else if( cat == "inclusive") histos.at(i)->Write(tmp, TObject::kWriteDelete);
            else continue;
          }
          else histos.at(i)->Write(tmp, TObject::kWriteDelete);
        }
      }
    }
    outfile->Close();
    outfile_name.str("");
  }
    
}

