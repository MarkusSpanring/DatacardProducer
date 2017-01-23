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
  
  files.push_back({Parameter.dataset.Z,s_Z});
  files.push_back({Parameter.dataset.EWKZ,s_EWKZ});
  files.push_back({Parameter.dataset.W,s_W});
  files.push_back({Parameter.dataset.TT,s_TT});
  files.push_back({Parameter.dataset.VV,s_VV});
  if(channel=="mt")files.push_back({Parameter.dataset.data_mt,s_data});
  if(channel=="et")files.push_back({Parameter.dataset.data_et,s_data});
  if(channel=="tt")files.push_back({Parameter.dataset.data_tt,s_data});
  files.push_back({Parameter.dataset.ggH,s_ggH});
  files.push_back({Parameter.dataset.qqH,s_qqH});
  if(ptShift){
    files.push_back({Parameter.dataset.ZtauUp,s_ZtauUp});
    files.push_back({Parameter.dataset.ZtauDown,s_ZtauDown});
    files.push_back({Parameter.dataset.EWKZtauUp,s_EWKZtauUp});
    files.push_back({Parameter.dataset.EWKZtauDown,s_EWKZtauDown});
    files.push_back({Parameter.dataset.TTtauUp,s_TTtauUp});
    files.push_back({Parameter.dataset.TTtauDown,s_TTtauDown});
    files.push_back({Parameter.dataset.VVtauUp,s_VVtauUp});
    files.push_back({Parameter.dataset.VVtauDown,s_VVtauDown});
    files.push_back({Parameter.dataset.ggHtauUp,s_ggHtauUp});
    files.push_back({Parameter.dataset.ggHtauDown,s_ggHtauDown});
    files.push_back({Parameter.dataset.qqHtauUp,s_qqHtauUp});
    files.push_back({Parameter.dataset.qqHtauDown,s_qqHtauDown});
  }
  if(jecShift){
    files.push_back({Parameter.dataset.Z,s_ZjecUp});
    files.push_back({Parameter.dataset.Z,s_ZjecDown});
    files.push_back({Parameter.dataset.EWKZ,s_EWKZjecUp});
    files.push_back({Parameter.dataset.EWKZ,s_EWKZjecDown});
    files.push_back({Parameter.dataset.W,s_WjecUp});
    files.push_back({Parameter.dataset.W,s_WjecDown});
    files.push_back({Parameter.dataset.TT,s_TTjecUp});
    files.push_back({Parameter.dataset.TT,s_TTjecDown});
    files.push_back({Parameter.dataset.VV,s_VVjecUp});
    files.push_back({Parameter.dataset.VV,s_VVjecDown});
    files.push_back({Parameter.dataset.ggH,s_ggHjecUp});
    files.push_back({Parameter.dataset.ggH,s_ggHjecDown});
    files.push_back({Parameter.dataset.qqH,s_qqHjecUp});
    files.push_back({Parameter.dataset.qqH,s_qqHjecDown});
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
      if( isZFile(files[i][1]) || isEWKZFile(files[i][1]) ) weight *= NtupleView->ZWeight;
      if( isTTFile(files[i][1]) ) weight *= NtupleView->topWeight;
      
      if(channel == "et") weight = weight * this->getAntiLep_tauscaling();

      for(auto cat : cats){

        for(auto strVar : vars){

          var = -999;

          if(strVar == s_mvis)                                var = NtupleView->m_vis;
          else if(strVar == s_msv)                            var = NtupleView->m_sv;
          else if(strVar == s_ptsv)                           var = NtupleView->pt_sv;
          else if(strVar == s_mt1)                            var = this->getMT();
          else if(strVar == s_jpt1)                           var = NtupleView->jpt_1;
          else if(strVar == s_jpt2)                           var = NtupleView->jpt_2;
          
          else if(strVar == s_pt1)                            var = NtupleView->pt_1;
          else if(strVar == s_pt2)                            var = NtupleView->pt_2;
          else if(strVar == s_eta1)                           var = NtupleView->eta_1;
          else if(strVar == s_eta2)                           var = NtupleView->eta_2;
          else if(strVar == s_met)                            var = NtupleView->met;
          else if(strVar == s_mttot)                          var = this->getMTTOT();
          else if(strVar == s_Hpt)                            var = this->CalcHPt(); 


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

          else if(strVar == "mjj")          var = this->getMjj();

          else if(strVar == "jeta1eta2"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = NtupleView->jeta1eta2;

          else if(strVar == "jdeta"
             && NtupleView->jpt_1 > 30
             && NtupleView->jpt_2 > 30)     var = this->getJdeta();

          else continue;


          if( isZFile(files[i][1]) )                     this->DYSelections(var, weight, cat, strVar, files[i][1]);
          
          else if( isEWKZFile(files[i][1]) )             this->EWKZSelections(var, weight, cat, strVar, files[i][1]);
          
          else if( isTTFile(files[i][1]) )               this->TSelections(var, weight, cat, strVar, files[i][1]);
          
          else if( isVVFile(files[i][1]) )               this->VVSelections(var, weight, cat, strVar, files[i][1]);

          else if( isWFile(files[i][1]) )                this->WSelections(var, weight, cat, strVar, files[i][1]);

          else if( files[i][1] == "data" )               this->dataSelections(var, 1., cat, strVar, files[i][1]);

          else if( isSignalFile(files[i][1]) )           this->signalSelections(var, weight, cat, strVar, files[i][1]);
          

          if(do2DFit){

            if ( !this->is2DCategories(cat) ) continue;
            
            if( isZFile(files[i][1]) )                   this->DYSelections(var, weight, cat, strVar, files[i][1], "2D");

            else if( isEWKZFile(files[i][1]) )           this->EWKZSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if( isTTFile(files[i][1]) )             this->TSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if( isVVFile(files[i][1]) )             this->VVSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if( isWFile(files[i][1]) )              this->WSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if( files[i][1] == "data" )             this->dataSelections(var, 1., cat, strVar, files[i][1], "2D");
              
            else if( isSignalFile(files[i][1]) )         this->signalSelections(var, weight, cat, strVar, files[i][1], "2D");
            
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
//UGLIEST FUNCTION EVER WRITTEN -> PLEASE SIMPLIFY!!!!!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
          tmp.ReplaceAll("jecUp",s_CMSjecScale+"13TeVUp");
          tmp.ReplaceAll("jecDown",s_CMSjecScale+"13TeVDown");
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

int CreateHistos::isZFile(TString fileName){
  if(fileName == s_Z) return 1;
  if(fileName == s_ZtauUp) return 1;
  if(fileName == s_ZtauDown) return 1;
  if(fileName == s_ZjecUp) return 1;
  if(fileName == s_ZjecDown) return 1;

  return 0;
}

int CreateHistos::isEWKZFile(TString fileName){
  if(fileName == s_EWKZ) return 1;
  if(fileName == s_EWKZtauUp) return 1;
  if(fileName == s_EWKZtauDown) return 1;
  if(fileName == s_EWKZjecUp) return 1;
  if(fileName == s_EWKZjecDown) return 1;

  return 0;
}

int CreateHistos::isSignalFile(TString fileName){
  if(fileName == s_ggH) return 1;
  if(fileName == s_ggHtauUp) return 1;
  if(fileName == s_ggHtauDown) return 1;
  if(fileName == s_ggHjecUp) return 1;
  if(fileName == s_ggHjecDown) return 1;
  if(fileName == s_qqH) return 1;
  if(fileName == s_qqHtauUp) return 1;
  if(fileName == s_qqHtauDown) return 1;
  if(fileName == s_qqHjecUp) return 1;
  if(fileName == s_qqHjecDown) return 1;

  return 0;
}

int CreateHistos::isTTFile(TString fileName){
  if(fileName == s_TT) return 1;
  if(fileName == s_TTtauUp) return 1;
  if(fileName == s_TTtauDown) return 1;
  if(fileName == s_TTjecUp) return 1;
  if(fileName == s_TTjecDown) return 1;

  return 0;  
}

int CreateHistos::isWFile(TString fileName){
  if(fileName == s_W) return 1;
  if(fileName == s_WjecUp) return 1;
  if(fileName == s_WjecDown) return 1;

  return 0;
}

int CreateHistos::isVVFile(TString fileName){
  if(fileName == s_VV) return 1;
  if(fileName == s_VVtauUp) return 1;
  if(fileName == s_VVtauDown) return 1;
  if(fileName == s_VVjecUp) return 1;
  if(fileName == s_VVjecDown) return 1;

  return 0;
}
