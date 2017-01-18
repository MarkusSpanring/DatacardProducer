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
      if(files[i][1] == "Z" || files[i][1] == "EWKZ") weight *= NtupleView->ZWeight;
      if(files[i][1] == "TT") weight *= NtupleView->topWeight;

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

          else if(strVar == "mjj")          var = NtupleView->mjj;

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

          if(files[i][1] == "EWKZ" 
             || files[i][1] == "EWKZtauUp"
             || files[i][1] == "EWKZtauDown")       this->EWKZSelections(var, weight, cat, strVar, files[i][1]);

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
          

          if(do2DFit){

            if ( !this->is2DCategories(cat) ) continue;
            
            if(files[i][1] == "Z" 
               || files[i][1] == "ZtauUp"
               || files[i][1] == "ZtauDown")       this->DYSelections(var, weight, cat, strVar, files[i][1], "2D");

            else if(files[i][1] == "EWKZ" 
               || files[i][1] == "EWKZtauUp"
               || files[i][1] == "EWKZtauDown")       this->EWKZSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if(files[i][1] == "TT"
                    || files[i][1] == "TTtauUp"
                    || files[i][1] == "TTtauDown")      this->TSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if(files[i][1] == "VV"
                    || files[i][1] == "VVtauUp"
                    || files[i][1] == "VVtauDown")      this->VVSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "W" )           this->WSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "data" )        this->dataSelections(var, 1., cat, strVar, files[i][1], "2D");
              
            else if(files[i][1] == "qqH"
                    || files[i][1] == "ggH"
                    || files[i][1] == "qqHtauUp"
                    || files[i][1] == "qqHtauDown"
                    || files[i][1] == "ggHtauUp"
                    || files[i][1] == "ggHtauDown")     this->signalSelections(var, weight, cat, strVar, files[i][1], "2D");
            
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

float CreateHistos::CalcJdeta(){
    if(NtupleView->jeta_1 != -999 && NtupleView->jeta_2 != -999 ){
        return fabs( NtupleView->jeta_1 - NtupleView->jeta_2 );
    }
    else return -999;

}
float CreateHistos::CalcHPt(){

  TLorentzVector tau;
  TLorentzVector mu;
  TLorentzVector met;

  mu.SetPtEtaPhiM(NtupleView->pt_1, NtupleView->eta_1, NtupleView->phi_1, NtupleView->m_1);
  tau.SetPtEtaPhiM(NtupleView->pt_2, NtupleView->eta_2, NtupleView->phi_2, NtupleView->m_2);
  if(useMVAMET)met.SetPtEtaPhiM(NtupleView->mvamet,0.,NtupleView->mvametphi,0.);
  else met.SetPtEtaPhiM(NtupleView->met,0.,NtupleView->metphi,0.);

  return (tau + mu + met).Pt();
}

double CreateHistos::getMT(){
  if(useMVAMET) return NtupleView->mt_1;
  return NtupleView->pfmt_1;
}

double CreateHistos::getMT2(){
  if(useMVAMET) return NtupleView->mt_2;
  return NtupleView->pfmt_2;
}

double CreateHistos::getMTTOT(){
  return TMath::Sqrt( TMath::Power(this->getMT(),2) + TMath::Power(this->getMT2(),2) + 2*NtupleView->pt_1*NtupleView->pt_2*( 1-TMath::Cos( TVector2::Phi_mpi_pi( NtupleView->phi_1-NtupleView->phi_2 ) ) ) );
}

void CreateHistos::getFFInputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( NtupleView->njets );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( this->getMT() );
  inputs.push_back( NtupleView->iso_1 );
}
void CreateHistos::getFF1Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->decayMode_1 );
  inputs.push_back( NtupleView->njets );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}
void CreateHistos::getFF2Inputs(vector<double>&inputs){
  inputs.push_back( NtupleView->pt_2 );
  inputs.push_back( NtupleView->pt_1 );
  inputs.push_back( NtupleView->decayMode_2 );
  inputs.push_back( NtupleView->njets );
  inputs.push_back( NtupleView->m_vis );
  inputs.push_back( 0 );
}


void CreateHistos::applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend){

  TString sub = extend + "+" + strVar +"_" + cat + "+";
  float usedVar=var;
  if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

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
      int binX=this->returnBin(Parameter.variable2D.D2_vbf.binsX,NtupleView->mjj);
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
      int binX=this->returnBin(Parameter.variable2D_tt.D2_vbf.binsX,NtupleView->mjj);
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
          this->GetHistbyName("TTL"+sub,strVar)->Fill(usedVar, weight);
          this->GetHistbyName("TTL_CMS_htt_ttbarShape_13TeVUp"+sub,strVar)->Fill(usedVar, weight*NtupleView->topWeight);
          this->GetHistbyName("TTL_CMS_htt_ttbarShape_13TeVDown"+sub,strVar)->Fill(usedVar, weight/NtupleView->topWeight);
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
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend){

    TString sub = extend + "+" + strVar +"_" + cat + "+";
    float usedVar=var;
    if(extend=="2D") usedVar = this->get2DVar(sub)+0.1;
    
    if(fname == "W"){
      if( this->Baseline("OS",cat) ){ 
        this->GetHistbyName("W_MC"+sub,strVar)->Fill(usedVar, weight);
        
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
          this->GetHistbyName("VVL"+sub,strVar)->Fill(usedVar, weight);
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
    for(int i = 1; i < entries; i++){
      if( this->GetHistbyName("QCD_OSW"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("QCD_OSW"+sub,strVar)->SetBinContent(i,0.);
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
    for(int i = 1; i < entries; i++){
      if( this->GetHistbyName("W"+sub,strVar)->GetBinContent(i) < 0 ) this->GetHistbyName("W"+sub,strVar)->SetBinContent(i,0.);
    }
  }

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";
  
  double W_SF = this->GetHistbyName("W_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_W"+sub,strVar)->Integral();

  double W_rel =  this->GetHistbyName("SS_Low_W"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_W"+sub,strVar)->Integral();
  double Z_rel =  this->GetHistbyName("SS_Low_Z"+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_Z"+sub,strVar)->Integral();
  double VV_rel = this->GetHistbyName("SS_Low_VV"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_VV"+sub,strVar)->Integral();
  double TT_rel = this->GetHistbyName("SS_Low_TT"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_TT"+sub,strVar)->Integral();
  double EWKZ_rel = this->GetHistbyName("SS_Low_EWKZ"+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_EWKZ"+sub,strVar)->Integral();

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
  }
  //cout << "NormUp: " << normUp << endl;
  //cout << "NormDown: " << normDown << endl;
  //cout << "//////////////////////////////////////" << endl;

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

//FIXME: temporary fix, has to be undone as soon as tauLepVeto is correct for etau channel
int CreateHistos::Vetos(){
  if(NtupleView->passesThirdLepVeto){
    if(channel == "et" && NtupleView->againstMuonLoose3_2 && NtupleView->againstElectronTightMVA6_2 && NtupleView->passesDiElectronVeto) return 1;
    if(channel == "mt" && NtupleView->passesTauLepVetos && NtupleView->passesDiMuonVeto) return 1;
    if(channel == "tt" && NtupleView->passesTauLepVetos) return 1; 
  }
  return 0;
}

/*int CreateHistos::Vetos(){
  if(NtupleView->passesTauLepVetos
     && NtupleView->passesThirdLepVeto){
    if(channel == "et" && NtupleView->passesDiElectronVeto) return 1;
    if(channel == "mt" && NtupleView->passesDiMuonVeto) return 1;
  }
  return 0;
  }*/

int CreateHistos::passMTCut(){
  if( !applyMTCut ) return 1;
  if( this->getMT() < Parameter.analysisCut.mTLow ) return 1;
  return 0;
}

int CreateHistos::passIso(TString type){
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::CategorySelection(TString cat, TString mtcut){


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
  if(cat == "0Jet") return this->Jet0(mtcut);
  
  ///////////////////////  boosted category     ///////////////////////////////
  if(cat == "boosted") return this->Boosted(mtcut);

  ///////////////////////  vbf category     ///////////////////////////////
  if(cat == "vbf") return this->VBF(mtcut);
  
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
int CreateHistos::VBF_low(TString mtcut){


  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
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
int CreateHistos::VBF_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
     && NtupleView->njets == 2
     && NtupleView->pt_2 > 20
     && NtupleView->mjj > 800
     && this->CalcHPt() > 100
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet1_low(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
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
int CreateHistos::Jet1_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
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
int CreateHistos::Jet0_low(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
     && NtupleView->njets == 0
     && NtupleView->pt_2 > 20
     && NtupleView->pt_2 < 50
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet0_high(TString mtcut){

  if( (this->getMT() < Parameter.analysisCut.mTLow
       || mtcut == "wo"
      )
     && NtupleView->njets == 0
     && NtupleView->pt_2 > 50
     )return 1;
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Jet0(TString mtcut){

  if( channel != "tt" ){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && NtupleView->njets == 0
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if( NtupleView->njets == 0
        && NtupleView->pt_1 > 50
        )return 1;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::Boosted(TString mtcut){

  if( channel != "tt" ){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && (NtupleView->njets == 1
            || (NtupleView->njets==2 && NtupleView->mjj<300)
            || NtupleView->njets>2
            )
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if(NtupleView->pt_1 > 50 
       && ( NtupleView->njets == 1
            || (NtupleView->njets>=2 && !( abs(NtupleView->jdeta) > 2.5 && this->CalcHPt() > 100 ) )
            )
       )return 1;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CreateHistos::VBF(TString mtcut){

  if( channel != "tt"){
    if( (this->getMT() < Parameter.analysisCut.mTLow
         || mtcut == "wo"
         )
        && NtupleView->njets == 2
        && NtupleView->mjj > 300
        && NtupleView->pt_2 > 30
        )return 1;
  }
  else{
    if(NtupleView->pt_1 > 50
       && NtupleView->njets >= 2
       && abs(NtupleView->jdeta) > 2.5
       && this->CalcHPt() > 100
       )return 1;
  }
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
void CreateHistos::writeHistos( TString channel, vector<TString> cats, vector<TString> vars ){

  stringstream outfile_name;
  TString sub;
  TString tmp;
  TString D2="";
  if(do2DFit) D2+="_2D";

  
  for(auto var : vars){
    outfile_name << "histos/"  << "htt_" << channel << ".inputs-sm-13TeV-"<<var<<D2<<".root";
    outfile = new TFile(outfile_name.str().c_str(), "RECREATE") ;

    for(auto cat : cats){
      outfile->mkdir(channel +"_"+ cat );
      outfile->cd(channel +"_"+ cat); 
      sub = "+" + var +"_" + cat + "+";

      for( int i=0;i<histo_names.size();i++ ){

        if(histo_names.at(i).Contains(sub) ){
          if(!keepFFDebugHistos
             && ( histo_names.at(i).Contains("data_jetFakes")
                  || histo_names.at(i).Contains("W_jetFakes")
                  || histo_names.at(i).Contains("Z_jetFakes")
                  || histo_names.at(i).Contains("TT_jetFakes")
                  || histo_names.at(i).Contains("VV_jetFakes")
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

