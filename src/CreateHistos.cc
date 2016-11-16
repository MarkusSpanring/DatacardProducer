#include "interface/CreateHistos.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TObject.h>

using namespace std;

CreateHistos::CreateHistos(){

  histos.clear();
  files.clear();
  histo_names.clear();
  
  files.push_back({Parameter.dataset.Z,"Z"});
  files.push_back({Parameter.dataset.W,"W"});
  files.push_back({Parameter.dataset.TT,"TT"});
  files.push_back({Parameter.dataset.VV,"VV"});
  files.push_back({Parameter.dataset.data,"data"});
  files.push_back({Parameter.dataset.ggH,"ggH"});
  files.push_back({Parameter.dataset.qqH,"qqH"});
  if(ptShift){
    files.push_back({Parameter.dataset.ZtauUp,"ZtauUp"});
    files.push_back({Parameter.dataset.ZtauDown,"ZtauDown"});
    files.push_back({Parameter.dataset.TTtauUp,"TTtauUp"});
    files.push_back({Parameter.dataset.TTtauDown,"TTtauDown"});
    files.push_back({Parameter.dataset.VVtauUp,"VVtauUp"});
    files.push_back({Parameter.dataset.VVtauDown,"VVtauDown"});
  }

  for(int i=0; i<variables.size(); i++) vars.push_back(variables.at(i));
  for(int i=0; i<categories.size(); i++) cats.push_back(categories.at(i));
  
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

void CreateHistos::initFakeFactors(){
  for(auto cat : cats){
    if( !this->is1DCategories(cat) ) continue;
    FFfile[cat] = TFile::Open("HTTutilities/Jet2TauFakes/data/"+channel+"/"+cat+"/fakeFactors_20161023.root");
    FFObj[cat] = (FakeFactor*)FFfile[cat]->Get("ff_comb");
  }
}

void CreateHistos::run(){


  //clearHistos();
  float weight = 1;
  float var = -999;
  
  initFakeFactors();
  

  for(int i =0;i < files.size();i++){

    this->loadFile(files[i][0]);
    Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
    //Int_t nentries = min( Int_t(NtupleView->fChain->GetEntries()), Int_t( 1000000 ) );
    cout<<"The input chain contains: "<<nentries<<" events."<<endl;
    float perc;
    for (Int_t jentry=0; jentry<100000;jentry++){
    //for (Int_t jentry=0; jentry<nentries;jentry++){       

      if(jentry % 200000 == 0){
        if(nentries > 0){
          perc =  ( (float)jentry / (float)nentries ) * 100;
        }
        cout.precision(3);
        cout<< "                                                             \r"<< flush;
        cout<< jentry << "/" << nentries <<"\t\t" << perc << "%\r"<< flush;
      }

      NtupleView->GetEntry(jentry);    

      weight = NtupleView->stitchedWeight*NtupleView->puweight*NtupleView->effweight*usedLuminosity;

      if(NtupleView->idisoweight_2 != 1) weight = weight * (0.9/0.83);
      weight = weight * this->getAntiLep_tauscaling();

      for(auto cat : cats){

        for(auto strVar : vars){

          var = -999;

          if(strVar == "m_vis")                              var = NtupleView->m_vis;
          else if(strVar == "mt_1")                          var = this->getMT();
          else if(strVar == "jpt_1")                         var = NtupleView->jpt_1;
          else if(strVar == "jpt_2")                         var = NtupleView->jpt_2;

          else if(strVar == "pt_1")                          var = NtupleView->pt_1;
          else if(strVar == "pt_2")                          var = NtupleView->pt_2;
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

          else if(strVar == "mjj"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = NtupleView->mjj;

          else if(strVar == "jeta1eta2"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = NtupleView->jeta1eta2;

          else if(strVar == "jdeta"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = this->CalcJdeta();

          else continue;



          if(files[i][1] == "Z" 
             || files[i][1] == "ZtauUp"
             || files[i][1] == "ZtauDown")       this->DYSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "TT"
             || files[i][1] == "TTtauUp"
             || files[i][1] == "TTtauDown")      this->TSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "VV"
             || files[i][1] == "VVtauUp"
             || files[i][1] == "VVtauDown")      this->VVSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "W" )           this->WSelections(var, weight, cat, strVar, files[i][1]);

          else if(files[i][1] == "data" )        this->dataSelections(var, 1., cat, strVar, files[i][1]);

          else if(files[i][1] == "qqH"
             || files[i][1] == "ggH"
             || files[i][1] == "qqHtauUp"
             || files[i][1] == "qqHtauDown"
             || files[i][1] == "ggHtauUp"
             || files[i][1] == "ggHtauDown")     this->signalSelections(var, weight, cat, strVar, files[i][1]);
        }
      }
      
    }
  }
  for(auto cat : cats){
    for(auto strVar : vars){
      this->Estimate_W_QCD(strVar, cat);
    }
  }
  if(calcFF){
    for(auto cat : cats){
      for(auto strVar : vars){
        this->EstimateFF(strVar, cat);
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
    return 1.0;

}

float CreateHistos::CalcJdeta(){
    if(NtupleView->jeta_1 != -999 && NtupleView->jeta_2 != -999 ){
        return fabs( NtupleView->jeta_1 - NtupleView->jeta_2 );
    }
    else return -999;

}
float CreateHistos::CalcHPt(){
  TLorentzVector tau;
  TLorentzVector mu;
  TLorentzVector MVAmet;

  mu.SetPtEtaPhiM(NtupleView->pt_1, NtupleView->eta_1, NtupleView->phi_1, NtupleView->m_1);
  tau.SetPtEtaPhiM(NtupleView->pt_2, NtupleView->eta_2, NtupleView->phi_2, NtupleView->m_2);
  MVAmet.SetPtEtaPhiM(NtupleView->mvamet,0.,NtupleView->mvametphi,0.);

  return (tau + mu + MVAmet).Pt();
}

double CreateHistos::getMT(){
  if(useMVAMET) return NtupleView->mt_1;
  return NtupleView->pfmt_1;
}


void CreateHistos::getFFInputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( NtupleView->njets );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( this->getMT() );
  inputs.push_back( NtupleView->iso_1 );

}

void CreateHistos::applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData){

  TString sub = "+" + strVar +"_" + cat + "+";
  
  if( this->Baseline("FF",cat) &&  is1DCategories(cat) ){
    if( isData || NtupleView->gen_match_2 < 6 ){
      FFinputs.clear();
      this->getFFInputs(FFinputs);
      this->GetHistbyName( fname+"_jetFakes"+sub,strVar)->Fill(var, weight*FFObj[cat]->value(FFinputs) );
      if(channel=="mt"){
        for( auto syst : Parameter.FFsystematics.mt.syst ){
          this->GetHistbyName( fname+"_jetFakes_"+syst+sub,strVar)->Fill(var, weight*FFObj[cat]->value(FFinputs, syst) );
        }
      }else if(channel=="et"){
        for( auto syst : Parameter.FFsystematics.et.syst ){
          this->GetHistbyName( fname+"_jetFakes_"+syst+sub,strVar)->Fill(var, weight*FFObj[cat]->value(FFinputs, syst) );
        }
      }else if(channel=="tt"){
        for( auto syst : Parameter.FFsystematics.tt.syst ){
          this->GetHistbyName( fname+"_jetFakes_"+syst+sub,strVar)->Fill(var, weight*FFObj[cat]->value(FFinputs, syst) );
        }
      }
    }
  }

}

void CreateHistos::DYSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";

    if(fname == "Z"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("Z"+sub,strVar)->Fill(var, weight);
        ////////////////////////////////////////////////////////////////
        if(NtupleView->gen_match_2 < 5){
          this->GetHistbyName("ZLL"+sub,strVar)->Fill(var, weight);
          this->GetHistbyName("ZL"+sub,strVar)->Fill(var, weight);
        }
        else if(NtupleView->gen_match_2 == 5){
          this->GetHistbyName("ZTT"+sub,strVar)->Fill(var, weight);
          this->GetHistbyName("ZTT_CMS_htt_dyShape_13TeVUp"+sub,strVar)->Fill(var, NtupleView->ZWeight * NtupleView->ZWeight);
          this->GetHistbyName("ZTT_CMS_htt_dyShape_13TeVDown"+sub,strVar)->Fill(var, 1);
        }
        else if(NtupleView->gen_match_2 == 6){
          this->GetHistbyName("ZLL"+sub,strVar)->Fill(var, weight);
          this->GetHistbyName("ZJ"+sub,strVar)->Fill(var, weight);
        }
        ////////////////////////////////////////////////////////////////
        if( NtupleView->NUP == 0)          this->GetHistbyName("Z_0Jets"+sub,strVar)->Fill(var, weight);
        else if( NtupleView->NUP == 1)     this->GetHistbyName("Z_1Jets"+sub,strVar)->Fill(var, weight);  
        else if(  NtupleView->NUP > 1)     this->GetHistbyName("Z_ge2Jets"+sub,strVar)->Fill(var, weight);
        ////////////////////////////////////////////////////////////////
      }
      if(calcFF) this->applyFF(var,weight,cat,strVar,fname,0);
      
      if( this->OS_W(cat) )                        this->GetHistbyName("OS_W_Z"+sub,strVar)->Fill(var, weight);
      if( this->SS_W(cat) )                        this->GetHistbyName("SS_W_Z"+sub,strVar)->Fill(var, weight);
      if( this->SS_Low(cat) )                      this->GetHistbyName("SS_Low_Z"+sub,strVar)->Fill(var, weight);
      if( this->SS_Low_relaxed(cat) )              this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Fill(var, weight);
      
    }
    else if(fname == "ZtauUp"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )       this->GetHistbyName("ZTT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ZtauDown"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )       this->GetHistbyName("ZTT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::TSelections(float var, float weight, TString cat, TString strVar, TString fname){
    TString sub = "+" + strVar +"_" + cat + "+";

    if(fname == "TT"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("TT"+sub,strVar)->Fill(var, weight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(var, NtupleView->topWeight * NtupleView->topWeight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(var, 1);
        if(NtupleView->gen_match_2 < 5)       this->GetHistbyName("TTL"+sub,strVar)->Fill(var, weight);
        else if(NtupleView->gen_match_2 == 5)  this->GetHistbyName("TTT"+sub,strVar)->Fill(var, weight);
        else if(NtupleView->gen_match_2 == 6)  this->GetHistbyName("TTJ"+sub,strVar)->Fill(var, weight);
      }

      if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_TT"+sub,strVar)->Fill(var, weight);
      if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_TT"+sub,strVar)->Fill(var, weight);
      if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_TT"+sub,strVar)->Fill(var, weight);

      if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Fill(var, weight);

      if(calcFF) this->applyFF(var,weight,cat,strVar,fname,0);
      
    }
    else if(fname == "TTtauUp"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("TTT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "TTtauDown"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("TTT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::WSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "W"){
      //if( this->Baseline("OS",cat) )         this->GetHistbyName("W"+sub,strVar)->Fill(var, weight);

      if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_W" + sub, strVar)->Fill(var, weight);
      if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_W" + sub, strVar)->Fill(var, weight);
      if( this->relaxed_W(cat, "low") )     this->GetHistbyName("relaxed_W_low_W" + sub, strVar)->Fill(var, weight);
      if( this->relaxed_W(cat, "high") )    this->GetHistbyName("relaxed_W_high_W" + sub, strVar)->Fill(var, weight);
      if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_W" + sub, strVar)->Fill(var, weight);
      
      if( this->SS_Low_relaxed(cat) )       this->GetHistbyName("SS_Low_relaxed_W" + sub, strVar)->Fill(var, weight);
      
      if(calcFF) this->applyFF(var,weight,cat,strVar,fname,0);
      
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::VVSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "VV"){
      if( this->Baseline("OS",cat) ){
        this->GetHistbyName("VV"+sub,strVar)->Fill(var, weight);
        if( NtupleView->gen_match_2 < 5)       this->GetHistbyName("VVL"+sub,strVar)->Fill(var, weight);
        else if( NtupleView->gen_match_2 == 5)  this->GetHistbyName("VVT"+sub,strVar)->Fill(var, weight);
        else if( NtupleView->gen_match_2 == 6)  this->GetHistbyName("VVJ"+sub,strVar)->Fill(var, weight);
      }

      if( this->OS_W(cat) )                   this->GetHistbyName("OS_W_VV"+sub,strVar)->Fill(var, weight);
      if( this->SS_W(cat) )                   this->GetHistbyName("SS_W_VV"+sub,strVar)->Fill(var, weight);
      if( this->SS_Low(cat) )                 this->GetHistbyName("SS_Low_VV"+sub,strVar)->Fill(var, weight);
      
      if( this->SS_Low_relaxed(cat) )         this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Fill(var, weight);

      if(calcFF) this->applyFF(var,weight,cat,strVar,fname,0);
      
    }
    else if(fname == "VVtauUp"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("VVT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "VVtauDown"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("VVT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::signalSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "ggH"
       || fname == "qqH"){

      if( this->Baseline("OS",cat) )                   this->GetHistbyName(fname+"125"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ggHtauUp"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("ggH_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ggHtauDown"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("ggH_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "qqHtauUp"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("qqH_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "qqHtauDown"){
      if( this->Baseline("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("qqH_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::dataSelections(float var, float weight, TString cat, TString strVar, TString fname){
 
    TString sub = "+" + strVar +"_" + cat + "+";

    if( this->Baseline("OS",cat) )        this->GetHistbyName("data_obs"+sub,strVar)->Fill(var, weight);

    if( this->OS_W(cat) )                 this->GetHistbyName("OS_W_data"+sub,strVar)->Fill(var, weight);
    if( this->SS_W(cat) )                 this->GetHistbyName("SS_W_data"+sub,strVar)->Fill(var, weight);
    if( this->SS_Low(cat) )               this->GetHistbyName("SS_Low_data"+sub,strVar)->Fill(var, weight);

    if(calcFF) this->applyFF(var,weight,cat,strVar,fname,1);
    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD_osw(TString strVar, TString cat){
  TString sub = "+" + strVar +"_" + cat + "+";

  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_data"+sub,strVar)   );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_VV"+sub,strVar), -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_Z"+sub,strVar),  -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_TT"+sub,strVar), -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Add( this->GetHistbyName("SS_W_W"+sub,strVar),  -1 );
  this->GetHistbyName("QCD_OSW"+sub,strVar)->Scale( this->QCD_OSSS(cat) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateW(TString strVar, TString cat){
  TString sub = "+" + strVar +"_" + cat + "+";

  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_data"+sub,strVar)   );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_VV"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_Z"+sub,strVar),  -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_TT"+sub,strVar), -1 );
  this->GetHistbyName("W_OSW"+sub,strVar)->Add( this->GetHistbyName("QCD_OSW"+sub,strVar), -1 );

  this->GetHistbyName("HL"+sub, strVar)->Add( this->GetHistbyName("relaxed_W_low_W"+sub,strVar) );
  this->GetHistbyName("HL"+sub, strVar)->Divide( this->GetHistbyName("relaxed_W_high_W"+sub,strVar) );

  this->GetHistbyName("W"+sub,strVar)->Add( this->GetHistbyName("W_OSW"+sub,strVar) );
  this->GetHistbyName("W"+sub,strVar)->Multiply( this->GetHistbyName("HL"+sub, strVar) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD(TString strVar, TString cat){
  TString sub = "+" + strVar +"_" + cat + "+";
  
  double W_SF = this->GetHistbyName("W_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_W"+sub,strVar)->Integral();

  double W_rel =  this->GetHistbyName("SS_Low_W"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Integral();
  double Z_rel =  this->GetHistbyName("SS_Low_Z"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Integral();
  double VV_rel = this->GetHistbyName("SS_Low_VV"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Integral();
  double TT_rel = this->GetHistbyName("SS_Low_TT"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Integral();

  this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Scale(  W_SF   );
  this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Scale(  W_rel  );
  this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Scale( VV_rel );
  this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Scale(  Z_rel  );
  this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Scale( TT_rel );

  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_data"+sub,strVar)  );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar), -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar),  -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar), -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar),  -1 );
  this->GetHistbyName("QCD"+sub,strVar)->Scale( this->QCD_OSSS(cat) );

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::Estimate_W_QCD(TString strVar, TString cat){
  this->CreateQCD_osw(strVar, cat);
  this->CreateW(strVar, cat);
  this->CreateQCD(strVar, cat);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::EstimateFF(TString strVar, TString cat){
  TString sub = "+" + strVar +"_" + cat + "+";

  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("data_jetFakes"+sub,strVar)   );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("W_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes"+sub,strVar),  -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes"+sub,strVar), -1 );
  this->GetHistbyName("jetFakes"+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes"+sub,strVar), -1 );

  if(channel=="mt"){
    for( auto syst : Parameter.FFsystematics.mt.syst ){
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("data_jetFakes_"+syst+sub,strVar)   );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("W_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_"+syst+sub,strVar),  -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes_"+syst+sub,strVar), -1 );
    }
  }
  if(channel=="et"){
    for( auto syst : Parameter.FFsystematics.et.syst ){
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("data_jetFakes_"+syst+sub,strVar)   );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("W_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_"+syst+sub,strVar),  -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes_"+syst+sub,strVar), -1 );
    }
  }
  if(channel=="tt"){
    for( auto syst : Parameter.FFsystematics.tt.syst ){
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("data_jetFakes_"+syst+sub,strVar)   );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("W_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("TT_jetFakes_"+syst+sub,strVar),  -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("Z_jetFakes_"+syst+sub,strVar), -1 );
      this->GetHistbyName("jetFakes_"+syst+sub,strVar)->Add( this->GetHistbyName("VV_jetFakes_"+syst+sub,strVar), -1 );
    }
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CreateHistos::Baseline(TString sign, TString cat){

    
    if( NtupleView->q_1 * NtupleView->q_2 < 0
       && this->passMTCut()
       && this->passIso("base")
       && this->Vetos()
       && this->CategorySelection(cat)
        ){
      if( sign=="OS"
          && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
          ) return 1;
      if( sign=="FF"
          && !NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
          && NtupleView->byVLooseIsolationMVArun2v1DBoldDMwLT_2
          ) return 1;
    } 
    
   
  return 0;
}


int CreateHistos::OS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat)
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::SS_W(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("base")
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat)
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::relaxed_W(TString cat, TString mt){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && this->passIso("relaxed")
     && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
     && this->Vetos()
     && this->CategorySelection(cat)){

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
     && this->CategorySelection(cat)
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::SS_Low_relaxed(TString cat){

  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && this->passIso("relaxed")
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
     && this->CategorySelection(cat)
     && this->Vetos()) return 1;
  return 0;
}


int CreateHistos::Vetos(){
  if(NtupleView->passesTauLepVetos
     && NtupleView->passesThirdLepVeto){
    if(channel == "et" && NtupleView->passesDiElectronVeto) return 1;
    if(channel == "mt" && NtupleView->passesDiMuonVeto) return 1;
  }
  return 0;
}

int CreateHistos::passMTCut(){
  if( !applyMTCut ) return 1;
  if( this->getMT() < Parameter.analysisCut.mTLow ) return 1;
  return 0;
}

int CreateHistos::passIso(TString type){
  if(type == "base"){
    if(channel == "et" && NtupleView->iso_1 < Parameter.analysisCut.elIso_base) return 1;
    if(channel == "mt" && NtupleView->iso_1 < Parameter.analysisCut.muIso_base) return 1;
  }
  else if(type == "relaxed"){
    if(channel == "et" && NtupleView->iso_1 < Parameter.analysisCut.elIso_relaxed) return 1;
    if(channel == "mt" && NtupleView->iso_1 < Parameter.analysisCut.muIso_relaxed) return 1;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::CategorySelection(TString cat){


//////////////////  Vertex dependendcy ////////////////////////////////////
  if(cat == "splitlowv"
     || cat == "splithighv"
     || cat == "lowv"
     || cat == "highv"){

    if( NtupleView->njets > 1
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
  if(cat == "VBF_low") return this->VBF_low();
  if(cat == "PUId_lo_VBF_low") return ( this->PUJetIdSelection("loose") & this->VBF_low() );
  if(cat == "PUId_me_VBF_low") return ( this->PUJetIdSelection("medium") & this->VBF_low() );
  if(cat == "PUId_ti_VBF_low") return ( this->PUJetIdSelection("tight") & this->VBF_low() );

/////////////////////////  VBF high category ///////////////////////////////
  if(cat == "VBF_high") return this->VBF_high();
  if(cat == "PUId_VBF_high") return ( this->PUJetIdSelection("tight") & this->VBF_high() );

/////////////////////////  1Jet low category ///////////////////////////////
  if(cat == "1Jet_low") return this->Jet1_low();
  if(cat == "PUId_1Jet_low") return ( this->PUJetIdSelection("tight") & this->Jet1_low() );

/////////////////////////  1Jet high category ///////////////////////////////
  if(cat == "1Jet_high") return this->Jet1_high();
  if(cat == "PUId_1Jet_high") return ( this->PUJetIdSelection("tight") & this->Jet1_high() );

/////////////////////////  2Jet low category ///////////////////////////////
  if(cat == "0Jet_low") return this->Jet0_low();
  if(cat == "PUId_0Jet_low") return ( this->PUJetIdSelection("tight") & this->Jet0_low() );

/////////////////////////  Jet high category ///////////////////////////////
  if(cat == "0Jet_high") return this->Jet0_high();  
  if(cat == "PUId_0Jet_high") return ( this->PUJetIdSelection("tight") & this->Jet0_high() );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::jet2_mvis(){

  if( NtupleView->njets > 1
      && NtupleView->jpt_1 > 30
      && NtupleView->jpt_2 > 30
      && this->getMT() < Parameter.analysisCut.mTLow
      && NtupleView->m_vis > 50
      && NtupleView->m_vis < 80)  return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::VBF_low(){


  if(this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->njets == 2
     && NtupleView->pt_2 > 20
     && NtupleView->mjj > 500
     && (NtupleView->mjj < 800
         || this->CalcHPt() < 100
        )
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::VBF_high(){

  if(this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->njets == 2
     && NtupleView->pt_2 > 20
     && NtupleView->mjj > 800
     && this->CalcHPt() > 100
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet1_low(){

  if(this->getMT() < Parameter.analysisCut.mTLow
     && (NtupleView->njets == 1 
         || (NtupleView->njets == 2
             && NtupleView->mjj < 500
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
int CreateHistos::Jet1_high(){

  if(this->getMT() < Parameter.analysisCut.mTLow
     && (NtupleView->njets == 1 
         || (NtupleView->njets == 2
             && NtupleView->mjj < 500
            ) 
        )
     && NtupleView->pt_2 > 40
     && this->CalcHPt() > 140             
        
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet0_low(){

  if(NtupleView->njets == 0
     && NtupleView->pt_2 > 20
     && NtupleView->pt_2 < 50
     && this->getMT() < Parameter.analysisCut.mTLow
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet0_high(){

  if(NtupleView->njets == 0
     && NtupleView->pt_2 > 50
     && this->getMT() < Parameter.analysisCut.mTLow
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::PUJetIdSelection(TString wp){
  if( wp == "tight" && NtupleView->jmva_1 > this->PUIdCutParamsTight(NtupleView->jeta_1) )    return 1;
  if( wp == "medium" && NtupleView->jmva_1 > this->PUIdCutParamsMedium(NtupleView->jeta_1) )  return 1;
  if( wp == "loose" && NtupleView->jmva_1 > this->PUIdCutParamsLoose(NtupleView->jeta_1) )    return 1;
  return 0;
}

float CreateHistos::PUIdCutParamsTight(float eta){
  if( fabs(eta) < 2.5) return 0.62;
  else if( fabs(eta) < 2.75 ) return -0.21;
  else if( fabs(eta) < 3.0 )  return -0.07;
  else if( fabs(eta) < 5.0 ) return  -0.03;

  return -1.;
}

float CreateHistos::PUIdCutParamsMedium(float eta){
  if( fabs(eta) < 2.5) return -0.06;
  else if( fabs(eta) < 2.75 ) return -0.42;
  else if( fabs(eta) < 3.0 )  return -0.3;
  else if( fabs(eta) < 5.0 ) return  -0.23;

  return -1.;
}

float CreateHistos::PUIdCutParamsLoose(float eta){
  if( fabs(eta) < 2.5) return -0.92;
  else if( fabs(eta) < 2.75 ) return -0.56;
  else if( fabs(eta) < 3.0 )  return -0.44;
  else if( fabs(eta) < 5.0 ) return  -0.39;

  return -1.;
}

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

TH1D* CreateHistos::GetHistbyName(TString name, TString strVar){

   for(int i = 0;i<histo_names.size();i++){
      if(histo_names.at(i) == name) return histos.at(i);
   }

   return this->JITHistoCreator(name, strVar);

}

void CreateHistos::returnBinning(double* returnHisto, vector<double> input){
  for(int i=0; i<input.size(); i++) returnHisto[i]=input.at(i);
}
int CreateHistos::returnBins(vector<double> input){
  return input.size();
}
TH1D* CreateHistos::getBinnedHisto(TString name,vector<double> input){
  double binning[this->returnBins(input)];
  this->returnBinning(binning,input);
  TH1D* tmp=new TH1D(name,"",this->returnBins(input)-1,&binning[0]);
  return tmp;
}

TH1D* CreateHistos::JITHistoCreator(TString name, TString strVar){

  int nbins = 1;
  double nmin = 0;
  double nmax = 1;

  int usingVarBins = 0;
  
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
void CreateHistos::writeHistos( TString channel, vector<TString> cats, vector<TString> vars ){

  stringstream outfile_name;
  TString sub;
  TString tmp;


  
  for(auto var : vars){
    outfile_name << "histos/"  << "htt_mt.inputs-sm-13TeV-"<<var<<".root";
    outfile = new TFile(outfile_name.str().c_str(), "RECREATE") ;

    for(auto cat : cats){
      outfile->mkdir(channel +"_"+ cat );
      outfile->cd(channel +"_"+ cat); 
      sub = "+" + var +"_" + cat + "+";

      for( int i=0;i<histo_names.size();i++ ){

        if(histo_names.at(i).Contains(sub) ){
          tmp = histo_names.at(i);
          tmp.ReplaceAll(sub, "");
          histos.at(i)->SetName(tmp);
          histos.at(i)->Write(tmp, TObject::kWriteDelete);
        }
      }
    }
    outfile->Close();
    outfile_name.str("");
  }
    
}

