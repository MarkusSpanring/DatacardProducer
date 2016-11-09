#include "interface/CreateHistos.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TObject.h>

using namespace std;

CreateHistos::CreateHistos(){
  
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
}

void CreateHistos::loadFile(TString filename){

  TChain *tchain = new TChain("TauCheck");
  tchain->Add(filename);
  
  NtupleView = new ntuple(tchain);
  cout<<"File: "<<filename<<" loaded"<<endl;
 
}

void CreateHistos::run(){


  //clearHistos();
  float weight = 1;
  float var =0;
  float W_ratio = 1;
  float SS_OS_ratio = 1;


  for(int i =0;i < files.size();i++){
    this->loadFile(files[i][0]);
    Int_t nentries = Int_t(NtupleView->fChain->GetEntries());
    cout<<"The input chain contains: "<<nentries<<" events."<<endl;
    float perc;
    //for (Int_t jentry=0; jentry<nentries;jentry++){
    for (Int_t jentry=0; jentry<100000;jentry++){       

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


      for(auto cat : cats){

        for(auto strVar : vars){

          var = -999;

          if(strVar == "m_vis")                              var = NtupleView->m_vis;
          else if(strVar == "mt_1")                          var = this->getMT();
          else if(strVar == "jpt_1")                         var = NtupleView->jpt_1;
          else if(strVar == "jpt_2")                         var = NtupleView->jpt_2;



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
             && NtupleView->jpt_1 > 50)     var = NtupleView->jeta_1;

          else if(strVar == "jeta_2"
             && NtupleView->jpt_2 > 50)     var = NtupleView->jeta_2;

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
    if( cat == "inclusive" ){
      W_ratio = this->GetWRatio();
      cout << "---------------" <<endl;
      cout << "W ratio: " << W_ratio << endl;

      SS_OS_ratio = this->GetQCDRatio();
      cout << "QCD ratio: " << SS_OS_ratio << endl;
      cout << "---------------" <<endl;
    }
  }

  for(auto cat : cats){
    for(auto var : vars){
      TString sub = "+" + var +"_" + cat + "+";
      this->ExtractQCD(1.06, cat,  var);
      this->GetHistbyName("W"+sub)->Scale(0.92);
    }
  }
  cout << "Done running over events." << endl;
  writeHistos( channel, cats, vars );
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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



void CreateHistos::DYSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";

    if(fname == "Z"){
      if( this->Common("OS",cat) ){
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
      else if( this->Common("SS",cat) )          this->GetHistbyName("Z"+sub+"_SS",strVar)->Fill(var, weight);


      if(cat == "inclusive"){

        if( this->CR_W() )                        this->GetHistbyName("CR_W_Z","CR_mt_1")->Fill(this->getMT(), weight);
        if( this->CR_QCD() )                      this->GetHistbyName("CR_QCD_Z","iso_1")->Fill(NtupleView->iso_1, weight);
        if( this->SR_QCD() )                      this->GetHistbyName("SR_QCD_Z","iso_1")->Fill(NtupleView->iso_1, weight);
      }
    }
    else if(fname == "ZtauUp"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )       this->GetHistbyName("ZTT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ZtauDown"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )       this->GetHistbyName("ZTT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::TSelections(float var, float weight, TString cat, TString strVar, TString fname){
    TString sub = "+" + strVar +"_" + cat + "+";

    if(fname == "TT"){
      if( this->Common("OS",cat) ){
        this->GetHistbyName("TT"+sub,strVar)->Fill(var, weight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(var, NtupleView->topWeight * NtupleView->topWeight);
        this->GetHistbyName("TT_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(var, 1);
        if(NtupleView->gen_match_2 == 5)       this->GetHistbyName("TTT"+sub,strVar)->Fill(var, weight);
        else if(NtupleView->gen_match_2 != 5)  this->GetHistbyName("TTJ"+sub,strVar)->Fill(var, weight);
      }
      else if( this->Common("SS",cat) )     this->GetHistbyName("TT"+sub+"_SS",strVar)->Fill(var, weight);

      if(cat == "inclusive"){
        if( this->CR_W() )                 this->GetHistbyName("CR_W_TT","CR_mt_1")->Fill(this->getMT(), weight);
        if( this->CR_QCD() )               this->GetHistbyName("CR_QCD_TT","iso_1")->Fill(NtupleView->iso_1, weight);
        if( this->SR_QCD() )               this->GetHistbyName("SR_QCD_TT","iso_1")->Fill(NtupleView->iso_1, weight);
      }
    }
    else if(fname == "TTtauUp"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("TTT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "TTtauDown"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("TTT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::WSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "W"){
      if( this->Common("OS",cat) )         this->GetHistbyName("W"+sub,strVar)->Fill(var, weight);
      else if( this->Common("SS",cat) )    this->GetHistbyName("W"+sub+"_SS",strVar)->Fill(var, weight);

      if(cat == "inclusive"){
        if( this->CR_W() )                 this->GetHistbyName("CR_W_W","CR_mt_1")->Fill(this->getMT(), weight);
        if( this->CR_QCD() )               this->GetHistbyName("CR_QCD_W","iso_1")->Fill(NtupleView->iso_1, weight);
        if( this->SR_QCD() )               this->GetHistbyName("SR_QCD_W","iso_1")->Fill(NtupleView->iso_1, weight);
      }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::VVSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "VV"){
      if( this->Common("OS",cat) ){
        this->GetHistbyName("VV"+sub,strVar)->Fill(var, weight);
        if( NtupleView->gen_match_2 == 5)       this->GetHistbyName("VVT"+sub,strVar)->Fill(var, weight);
        else if( NtupleView->gen_match_2 != 5)  this->GetHistbyName("VVJ"+sub,strVar)->Fill(var, weight);
      }
      else if( this->Common("SS",cat) )         this->GetHistbyName("VV"+sub+"_SS",strVar)->Fill(var, weight);

      if(cat == "inclusive"){
        if( this->CR_W() )                   this->GetHistbyName("CR_W_VV","CR_mt_1")->Fill(this->getMT(), weight);
        if( this->CR_QCD() )                 this->GetHistbyName("CR_QCD_VV","iso_1")->Fill(NtupleView->iso_1, weight);
        if( this->SR_QCD() )                 this->GetHistbyName("SR_QCD_VV","iso_1")->Fill(NtupleView->iso_1, weight);
      }
    }
    else if(fname == "VVtauUp"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("VVT_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "VVtauDown"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("VVT_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::signalSelections(float var, float weight, TString cat, TString strVar, TString fname){

    TString sub = "+" + strVar +"_" + cat + "+";
    if(fname == "ggH"
       || fname == "qqH"){

      if( this->Common("OS",cat) )                   this->GetHistbyName(fname+"125"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ggHtauUp"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("ggH_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "ggHtauDown"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("ggH_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "qqHtauUp"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("qqH_CMS_scale_t_mt_13TeVUp"+sub,strVar)->Fill(var, weight);
    }
    else if(fname == "qqHtauDown"){
      if( this->Common("OS",cat)
          && NtupleView->gen_match_2 == 5 )          this->GetHistbyName("qqH_CMS_scale_t_mt_13TeVDown"+sub,strVar)->Fill(var, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::dataSelections(float var, float weight, TString cat, TString strVar, TString fname){
 
    TString sub = "+" + strVar +"_" + cat + "+";

    if( this->Common("OS",cat) )         this->GetHistbyName("data_obs"+sub,strVar)->Fill(var, weight);
    else if( this->Common("SS",cat) )    this->GetHistbyName("data_obs"+sub+"_SS",strVar)->Fill(var, weight);
    if(cat == "inclusive"){

      if( this->CR_W() )                 this->GetHistbyName("CR_W_data","CR_mt_1")->Fill(this->getMT(), weight);
      if( this->CR_QCD() )               this->GetHistbyName("CR_QCD_data","iso_1")->Fill(NtupleView->iso_1, weight);
      if( this->SR_QCD() )               this->GetHistbyName("SR_QCD_data","iso_1")->Fill(NtupleView->iso_1, weight);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float CreateHistos::GetWRatio(){

  this->GetHistbyName("CR_W_data")->Add(this->GetHistbyName("CR_W_TT"), -1);
  this->GetHistbyName("CR_W_data")->Add(this->GetHistbyName("CR_W_Z"), -1);
  this->GetHistbyName("CR_W_data")->Add(this->GetHistbyName("CR_W_VV"), -1);

  return ( this->GetHistbyName("CR_W_data")->Integral() /  this->GetHistbyName("CR_W_W")->Integral() ) ;
}

float CreateHistos::GetQCDRatio(){
  this->GetHistbyName("CR_QCD_data")->Add(this->GetHistbyName("CR_QCD_TT"), -1);
  this->GetHistbyName("CR_QCD_data")->Add(this->GetHistbyName("CR_QCD_Z"), -1);
  this->GetHistbyName("CR_QCD_data")->Add(this->GetHistbyName("CR_QCD_VV"), -1);
  this->GetHistbyName("CR_QCD_data")->Add(this->GetHistbyName("CR_QCD_W"), -1);

  this->GetHistbyName("SR_QCD_data")->Add(this->GetHistbyName("SR_QCD_TT"), -1);
  this->GetHistbyName("SR_QCD_data")->Add(this->GetHistbyName("SR_QCD_Z"), -1);
  this->GetHistbyName("SR_QCD_data")->Add(this->GetHistbyName("SR_QCD_VV"), -1);
  this->GetHistbyName("SR_QCD_data")->Add(this->GetHistbyName("SR_QCD_W"), -1);


  return ( this->GetHistbyName("SR_QCD_data")->Integral() /  this->GetHistbyName("CR_QCD_data")->Integral() );
}

void CreateHistos::ExtractQCD(float SS_OS_ratio, TString cat, TString strVar ){

  TString sub = "+" + strVar +"_" + cat + "+";

  this->GetHistbyName("QCD"+sub,strVar)->Add(this->GetHistbyName("data_obs"+sub+"_SS",strVar)  );
  this->GetHistbyName("QCD"+sub,strVar)->Add(this->GetHistbyName("TT"+sub+"_SS",strVar) , -1);
  this->GetHistbyName("QCD"+sub,strVar)->Add(this->GetHistbyName("W"+sub+"_SS",strVar) , -1);
  this->GetHistbyName("QCD"+sub,strVar)->Add(this->GetHistbyName("Z"+sub+"_SS",strVar) , -1);
  this->GetHistbyName("QCD"+sub,strVar)->Add(this->GetHistbyName("VV"+sub+"_SS",strVar) , -1);
  this->GetHistbyName("QCD"+sub,strVar)->Scale(SS_OS_ratio);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CreateHistos::Common(TString sign, TString cat){

    if(sign == "OS"
       && NtupleView->q_1 * NtupleView->q_2 < 0
       && NtupleView->iso_1 < 0.15
       && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
       && this->Vetos()
       && this->CategorySelection(cat)
      ) return 1;

    if(sign == "SS"
       && NtupleView->q_1 * NtupleView->q_2 > 0
       && !(cat == "VBF_low"
           || cat == "VBF_high"
           || cat == "1Jet_high"
          )
       && NtupleView->iso_1 < 0.15
       && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
       && this->Vetos()
       && this->CategorySelection(cat)
      ) return 1;

    if(sign == "SS"
       && NtupleView->q_1 * NtupleView->q_2 > 0
       && (cat == "VBF_low"
           || cat == "VBF_high"
           || cat == "1Jet_high"
          )
       && NtupleView->iso_1 < 0.3
       && NtupleView->byMediumIsolationMVArun2v1DBoldDMwLT_2
       && this->Vetos()
       && this->CategorySelection(cat)
      ) return 1;
   
  return 0;
}

int CreateHistos::CR_W(){

  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && NtupleView->iso_1 < 0.15
     && this->getMT() > Parameter.analysisCut.mTHigh
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->Vetos()) return 1;
  return 0;
}
int CreateHistos::CR_QCD(){
  if(NtupleView->q_1 * NtupleView->q_2 > 0
     && (NtupleView->iso_1 > 0.15
         && NtupleView->iso_1 < 0.4)
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->Vetos()) return 1;
  return 0;
}
int CreateHistos::SR_QCD(){
  if(NtupleView->q_1 * NtupleView->q_2 < 0
     && (NtupleView->iso_1 > 0.15
         && NtupleView->iso_1 < 0.4)
     && this->getMT() < Parameter.analysisCut.mTLow
     && NtupleView->byTightIsolationMVArun2v1DBoldDMwLT_2
     && this->Vetos()) return 1;
  return 0;
}

int CreateHistos::Vetos(){
  if(NtupleView->passesTauLepVetos
     && NtupleView->passesThirdLepVeto
     && NtupleView->passesDiMuonVeto) return 1;
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
      nmin = Parameter.variable.m_vis.nmin;;
      nmax = Parameter.variable.m_vis.nmax;
    }
  }

  if(strVar == "jpt_1"
     || strVar == "jpt_2"
     || strVar == "jpt_1_2"
     || strVar == "jpt_1_2p5"
     || strVar == "jpt_1_3"
     || strVar == "jpt_2_2"
     || strVar == "jpt_2_2p5"
     || strVar == "jpt_2_3"){
    nbins = 70;
    nmin = 0;
    nmax = 350;
  }
  else if(strVar == "jeta_1" || strVar == "jeta_2"){
    nbins = 40;
    nmin = -5;
    nmax = 5;
  }
  else if(strVar == "jdeta"){
    nbins = 40;
    nmin = 0;
    nmax = 10;
  }
  else if(strVar == "mt_1"){
    nbins = 40;
    nmin = 0;
    nmax = 40;
  }
  else if(strVar == "CR_mt_1"){
    nbins = 40;
    nmin = 70;
    nmax = 1000;
  }
  else if(strVar == "iso_1"){
    nbins = 25;
    nmin = 0;
    nmax = 0.5;
  }
  else if(strVar == "mjj"){
    nbins = 70;
    nmin = 0;
    nmax = 1400;
  }
  else if(strVar == "jeta1eta2"){
    nbins = 40;
    nmin = -10;
    nmax = 10;
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

double CreateHistos::getMT(){
  if(useMVAMET) return NtupleView->mt_1;
  else return NtupleView->pfmt_1;
}
