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
  if(channel=="mt" && !doMC )files.push_back({Parameter.dataset.data_mt,s_data});
  else if(channel=="mt" && doMC)files.push_back({Parameter.dataset.MCsum_mt,s_data});
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
  float weight_data = 1;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(doInitialize){

    for(auto cat : cats){
      for(auto strVar : vars){
        if( !do2DFit || cat == s_inclusive ){
          initDYSelections(cat,strVar);
          initTSelections(cat,strVar);
          initVVSelections(cat,strVar);
          initEWKZSelections(cat,strVar);
          initSignalSelections(cat,strVar);
        }
        else if( do2DFit ){
          if( channel == "tt"
              && (cat.Contains(s_wjets)
                  || cat.Contains(s_antiiso) )
              ) continue;
          initDYSelections(cat,strVar,"2D");
          initTSelections(cat,strVar,"2D");
          initVVSelections(cat,strVar,"2D");
          initEWKZSelections(cat,strVar,"2D");
          initSignalSelections(cat,strVar,"2D");
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int i =0;i < files.size();i++){

    if(files[i][1].Contains(s_jecUp)) this->isJEC=1;
    else if(files[i][1].Contains(s_jecDown)) this->isJEC=-1;
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
      if(!doMC){
        if( isZFile(files[i][1]) || isEWKZFile(files[i][1]) ) weight *= NtupleView->ZWeight;
        if( isTTFile(files[i][1]) ) weight *= NtupleView->topWeight;
      }
      else if(doMC) weight_data=weight;

      //FIXME: delete this function if not needed anymore
      //if(channel == "et") weight = weight * this->getAntiLep_tauscaling();

      for(auto cat : cats){

        if( !doMC && isZFile(files[i][1]) ) weight *= this->getZmumuWeight( cat );

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

          if( !do2DFit || cat == s_inclusive ){

            if( isZFile(files[i][1]) )                     this->DYSelections(var, weight, cat, strVar, files[i][1]);
            
            else if( isEWKZFile(files[i][1]) )             this->EWKZSelections(var, weight, cat, strVar, files[i][1]);
          
            else if( isTTFile(files[i][1]) )               this->TSelections(var, weight, cat, strVar, files[i][1]);
            
            else if( isVVFile(files[i][1]) )               this->VVSelections(var, weight, cat, strVar, files[i][1]);

            else if( isWFile(files[i][1]) )                this->WSelections(var, weight, cat, strVar, files[i][1]);

            else if( files[i][1] == s_data )               this->dataSelections(var, weight_data, cat, strVar, files[i][1]);

            else if( isSignalFile(files[i][1]) )           this->signalSelections(var, weight, cat, strVar, files[i][1]);

          }
          

          if(do2DFit){

            if ( !this->is2DCategories(cat) ) continue;
            if( channel == "tt"
                && (cat.Contains(s_wjets)
                    || cat.Contains(s_antiiso) )
                ) continue;
            
            if( isZFile(files[i][1]) )                   this->DYSelections(var, weight, cat, strVar, files[i][1], "2D");

            else if( isEWKZFile(files[i][1]) )           this->EWKZSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if( isTTFile(files[i][1]) )             this->TSelections(var, weight, cat, strVar, files[i][1], "2D");
            
            else if( isVVFile(files[i][1]) )             this->VVSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if( isWFile(files[i][1]) )              this->WSelections(var, weight, cat, strVar, files[i][1], "2D");
              
            else if( files[i][1] == s_data )             this->dataSelections(var, weight_data, cat, strVar, files[i][1], "2D");
              
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
        if( do2DFit && this->is2DCategories(cat) )      this->Estimate_W_QCD(strVar, cat, "2D"); 
      }
    }
  }
  if(calcFF){
    for(auto cat : cats){
      for(auto strVar : vars){
        if( !do2DFit || cat == s_inclusive )            this->EstimateFF(strVar, cat);
        else if( do2DFit && this->is2DCategories(cat) ) this->EstimateFF(strVar, cat, "2D");
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

  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_data+sub,strVar)   );
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_VV+sub,strVar), -1 );
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_Z+sub,strVar),  -1 );
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_TT+sub,strVar), -1 );
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_EWKZ+sub,strVar),  -1 );
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+sub,strVar) );
  this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVUp+sub,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+sub,strVar) );
  this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVDown+sub,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+"_"+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar) );
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+"_"+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar) );
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+sub,strVar),  -1 );
  this->GetHistbyName("OS_W_"+s_QCD+sub,strVar)->Scale( this->QCD_OSSS(cat) );

  this->resetZeroBins("OS_W_"+s_QCD+sub,strVar);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //WSFUncertainties
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVUp+sub,strVar), -1  );
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVUp+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  this->resetZeroBins("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVUp+sub,strVar);
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+s_WSFUncert+channel+"_"+s_13TeVDown+sub,strVar), -1  );
  this->GetHistbyName("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVDown+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  this->resetZeroBins("OS_W_"+s_QCD+s_WSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVDown+sub,strVar);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //QCDSFUncertainties
  this->GetHistbyName("OS_W_"+s_QCD+s_QCDSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar) );
  this->GetHistbyName("OS_W_"+s_QCD+s_QCDSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVUp+sub,strVar)->Scale(1+this->getQCDSFUncertainty(cat) );
  this->GetHistbyName("OS_W_"+s_QCD+s_QCDSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar) );
  this->GetHistbyName("OS_W_"+s_QCD+s_QCDSFUncert+channel+this->return2DString(cat)+"_"+s_13TeVDown+sub,strVar)->Scale(1-this->getQCDSFUncertainty(cat) );
  
  if(jecShift){
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_data+sub,strVar)   );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_VV+s_jecUp+sub,strVar), -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_Z+s_jecUp+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_TT+s_jecUp+sub,strVar), -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+s_jecUp+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_EWKZ+s_jecUp+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar)->Scale( this->QCD_OSSS(cat) );

   this->resetZeroBins("OS_W_"+s_QCDjecUp+sub,strVar);
   
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_data+sub,strVar)   );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_VV+s_jecDown+sub,strVar), -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_Z+s_jecDown+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_TT+s_jecDown+sub,strVar), -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_W+s_jecDown+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName("SS_W_"+s_EWKZ+s_jecDown+sub,strVar),  -1 );
   this->GetHistbyName("OS_W_"+s_QCDjecDown+sub,strVar)->Scale( this->QCD_OSSS(cat) );

   this->resetZeroBins("OS_W_"+s_QCDjecDown+sub,strVar);
   
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateW(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";

  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_data+sub,strVar)   );
  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_VV+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_Z+sub,strVar),  -1 );
  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_TT+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_EWKZ+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar), -1 );

  double w_normFactor = this->GetHistbyName(s_W+"_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_W+sub,strVar)->Integral();
  double w_normFactorSR = this->GetHistbyName("SR_MC"+s_W+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_low_"+s_W+sub,strVar)->Integral();
  double w_normFactorCR = this->GetHistbyName("OS_W_"+s_W+sub,strVar)->Integral() / this->GetHistbyName("relaxed_W_high_"+s_W+sub,strVar)->Integral();

  this->GetHistbyName(s_West+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+sub,strVar) );
  this->GetHistbyName(s_West+sub,strVar)->Scale( w_normFactor );
  this->GetHistbyName("SR_"+s_W+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+sub,strVar)->Scale( w_normFactorSR );
  this->GetHistbyName("WCR_"+s_W+sub,strVar)->Add( this->GetHistbyName("relaxed_W_high_"+s_W+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+sub,strVar)->Scale( w_normFactorCR );

  this->resetZeroBins(s_West+sub,strVar);
  this->resetZeroBins("SR_"+s_W+sub,strVar);
  this->resetZeroBins("WCR_"+s_W+sub,strVar);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //WSFUncertainties
  this->GetHistbyName(s_West+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_West+sub,strVar) );
  this->GetHistbyName(s_West+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName(s_West+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_West+sub,strVar) );
  this->GetHistbyName(s_West+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );
  this->GetHistbyName("SR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SR_"+s_W+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName("SR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SR_"+s_W+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );
  this->GetHistbyName("WCR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("WCR_"+s_W+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName("WCR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("WCR_"+s_W+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //QCDSFUncertainties
  this->GetHistbyName(s_West+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_West+sub,strVar) );
  this->GetHistbyName(s_West+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_West+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SR_"+s_W+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SR_"+s_W+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("WCR_"+s_W+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("WCR_"+s_W+sub,strVar) );
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //jetToTauFakeShift Up
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_data+sub,strVar)   );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_VV+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeUp+sub,strVar),  -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeUp+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_EWKZ+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar), -1 );

  w_normFactor = this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeUp+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeUp+sub,strVar)->Integral();
  w_normFactorSR = this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeUp+sub,strVar)->Integral() / this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeUp+sub,strVar)->Integral();
  w_normFactorCR = this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeUp+sub,strVar)->Integral() / this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeUp+sub,strVar)->Integral();

  this->GetHistbyName(s_West+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeUp+sub,strVar) );
  this->GetHistbyName(s_West+s_jetToTauFakeUp+sub,strVar)->Scale( w_normFactor );
  this->GetHistbyName("SR_"+s_W+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeUp+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_jetToTauFakeUp+sub,strVar)->Scale( w_normFactorSR );
  this->GetHistbyName("WCR_"+s_W+s_jetToTauFakeUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeUp+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_jetToTauFakeUp+sub,strVar)->Scale( w_normFactorCR );

  this->resetZeroBins(s_West+s_jetToTauFakeUp+sub,strVar);
  this->resetZeroBins("SR_"+s_W+s_jetToTauFakeUp+sub,strVar);
  this->resetZeroBins("WCR_"+s_W+s_jetToTauFakeUp+sub,strVar);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //jetToTauFakeShift Down
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_data+sub,strVar)   );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_VV+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_Z+s_jetToTauFakeDown+sub,strVar),  -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_TT+s_jetToTauFakeDown+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_EWKZ+sub,strVar), -1 );
  this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCD+sub,strVar), -1 );
  
  w_normFactor = this->GetHistbyName(s_W+"_OSW"+s_jetToTauFakeDown+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeDown+sub,strVar)->Integral();
  w_normFactorSR = this->GetHistbyName("SR_MC"+s_W+s_jetToTauFakeDown+sub,strVar)->Integral() / this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeDown+sub,strVar)->Integral();
  w_normFactorCR = this->GetHistbyName("OS_W_"+s_W+s_jetToTauFakeDown+sub,strVar)->Integral() / this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeDown+sub,strVar)->Integral();

  this->GetHistbyName(s_West+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeDown+sub,strVar) );
  this->GetHistbyName(s_West+s_jetToTauFakeDown+sub,strVar)->Scale( w_normFactor );
  this->GetHistbyName("SR_"+s_W+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_W+s_jetToTauFakeDown+sub,strVar) );
  this->GetHistbyName("SR_"+s_W+s_jetToTauFakeDown+sub,strVar)->Scale( w_normFactorSR );
  this->GetHistbyName("WCR_"+s_W+s_jetToTauFakeDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_high_"+s_W+s_jetToTauFakeDown+sub,strVar) );
  this->GetHistbyName("WCR_"+s_W+s_jetToTauFakeDown+sub,strVar)->Scale( w_normFactorCR );

  this->resetZeroBins(s_West+s_jetToTauFakeDown+sub,strVar);
  this->resetZeroBins("SR_"+s_W+s_jetToTauFakeDown+sub,strVar);
  this->resetZeroBins("WCR_"+s_W+s_jetToTauFakeDown+sub,strVar);

  if(jecShift){
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_data+sub,strVar)   );
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_VVjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_ZjecUp+sub,strVar),  -1 );
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_TTjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_EWKZjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar), -1 );
    
    w_normFactor = this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_WjecUp+sub,strVar)->Integral();
    w_normFactorSR = this->GetHistbyName("SR_MC"+s_WjecUp+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_low_"+s_WjecUp+sub,strVar)->Integral();
    w_normFactorCR = this->GetHistbyName("OS_W_"+s_WjecUp+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_WjecUp+sub,strVar)->Integral();
    
    this->GetHistbyName(s_WestjecUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_WjecUp+sub,strVar) );
    this->GetHistbyName(s_WestjecUp+sub,strVar)->Scale( w_normFactor );
    this->GetHistbyName("SR_"+s_WjecUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_WjecUp+sub,strVar) );
    this->GetHistbyName("SR_"+s_WjecUp+sub,strVar)->Scale( w_normFactorSR );
    this->GetHistbyName("WCR_"+s_WjecUp+sub,strVar)->Add( this->GetHistbyName("relaxed_W_high_"+s_WjecUp+sub,strVar) );
    this->GetHistbyName("WCR_"+s_WjecUp+sub,strVar)->Scale( w_normFactorCR );

    this->resetZeroBins(s_WestjecUp+sub,strVar);
    this->resetZeroBins("SR_"+s_WjecUp+sub,strVar);
    this->resetZeroBins("WCR_"+s_WjecUp+sub,strVar);
    
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_data+sub,strVar)   );
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_VVjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_ZjecUp+sub,strVar),  -1 );
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_TTjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_EWKZjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Add( this->GetHistbyName("OS_W_"+s_QCDjecUp+sub,strVar), -1 );
    
    w_normFactor = this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_WjecDown+sub,strVar)->Integral();
    w_normFactorSR = this->GetHistbyName("SR_MC"+s_WjecDown+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_low_"+s_WjecDown+sub,strVar)->Integral();
    w_normFactorCR = this->GetHistbyName("OS_W_"+s_WjecDown+sub,strVar)->Integral()  / this->GetHistbyName("relaxed_W_high_"+s_WjecDown+sub,strVar)->Integral();
    
    this->GetHistbyName(s_WestjecDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_WjecDown+sub,strVar) );
    this->GetHistbyName(s_WestjecDown+sub,strVar)->Scale( w_normFactor );
    this->GetHistbyName("SR_"+s_WjecDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_low_"+s_WjecDown+sub,strVar) );
    this->GetHistbyName("SR_"+s_WjecDown+sub,strVar)->Scale( w_normFactorSR );
    this->GetHistbyName("WCR_"+s_WjecDown+sub,strVar)->Add( this->GetHistbyName("relaxed_W_high_"+s_WjecDown+sub,strVar) );
    this->GetHistbyName("WCR_"+s_WjecDown+sub,strVar)->Scale( w_normFactorCR );
    
    this->resetZeroBins(s_WestjecDown+sub,strVar);
    this->resetZeroBins("SR_"+s_WjecDown+sub,strVar);
    this->resetZeroBins("WCR_"+s_WjecDown+sub,strVar);
    
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateHistos::CreateQCD(TString strVar, TString cat, TString extend){
  TString sub = extend + "+" + strVar +"_" + cat + "+";
  
  double W_SF = ( this->GetHistbyName("W_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("W_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_"+s_W+sub,strVar)->Integral() : 0;
  double W_rel = ( this->GetHistbyName("SS_Low_"+s_W+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_W+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar)->Integral() : 0;
  double Z_rel = ( this->GetHistbyName("SS_Low_"+s_Z+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_Z+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_Z+sub,strVar)->Integral() : 0;
  double VV_rel = ( this->GetHistbyName("SS_Low_"+s_VV+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_VV+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_VV+sub,strVar)->Integral() : 0;
  double TT_rel = ( this->GetHistbyName("SS_Low_"+s_TT+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_TT+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_TT+sub,strVar)->Integral() : 0;
  double EWKZ_rel = ( this->GetHistbyName("SS_Low_"+s_EWKZ+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_EWKZ+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_EWKZ+sub,strVar)->Integral() : 0;

  this->GetHistbyName("AICR_"+s_W+sub,strVar)->Add( this->GetHistbyName( "SR_MC"+s_W+sub,strVar ) );
  this->GetHistbyName("AICR_"+s_W+sub,strVar)->Scale(W_SF);

  if(cat.Contains(s_antiiso)){
    W_rel=1;
    Z_rel=1;
    VV_rel=1;
    TT_rel=1;
    EWKZ_rel=1;
  }
  
  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar)->Scale(  W_SF   );
  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar)->Scale(  W_rel  );
  this->GetHistbyName("SS_Low_relaxed_"+s_VV+sub,strVar)->Scale( VV_rel );
  this->GetHistbyName("SS_Low_relaxed_"+s_Z+sub,strVar)->Scale(  Z_rel  );
  this->GetHistbyName("SS_Low_relaxed_"+s_TT+sub,strVar)->Scale( TT_rel );
  this->GetHistbyName("SS_Low_relaxed_"+s_EWKZ+sub,strVar)->Scale( EWKZ_rel );

  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVUp,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar) );
  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVUp,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVDown,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar) );
  this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVDown,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );

  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_"+s_data+sub,strVar)  );
  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_VV+sub,strVar), -1 );
  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_Z+sub,strVar),  -1 );
  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_TT+sub,strVar), -1 );
  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_EWKZ+sub,strVar),  -1 );
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  this->GetHistbyName(s_QCDest+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub,strVar),  -1 );
  this->GetHistbyName(s_QCD+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCDest+sub,strVar)->Scale( this->QCD_OSSS(cat) );

  this->resetZeroBins(s_QCDest+sub,strVar);
  this->resetZeroBins(s_QCD+sub,strVar);
  this->resetZeroBins("AICR_"+s_W+sub,strVar);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //WSFUncertainties
  this->GetHistbyName("AICR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("AICR_"+s_W+sub,strVar) );
  this->GetHistbyName("AICR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1.+this->getWSFUncertainty(cat) );
  this->GetHistbyName("AICR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("AICR_"+s_W+sub,strVar) );
  this->GetHistbyName("AICR_"+s_W+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1.-this->getWSFUncertainty(cat) );
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVUp,strVar), -1 );
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVDown,strVar), -1 );
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  this->GetHistbyName(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( this->QCD_OSSS(cat) );
  this->resetZeroBins(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar);
  this->resetZeroBins(s_QCDest+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar);
  this->GetHistbyName(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVUp,strVar), -1 );
  this->GetHistbyName(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_W+sub+s_WSFUncert+channel+"_"+s_13TeVDown,strVar), -1 );
  this->resetZeroBins(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar);
  this->resetZeroBins(s_QCD+s_WSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //QCDSFUncertainties
  this->GetHistbyName("AICR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName("AICR_"+s_W+sub,strVar) );
  this->GetHistbyName("AICR_"+s_W+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName("AICR_"+s_W+sub,strVar) );
  this->GetHistbyName(s_QCDest+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCDest+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1+this->getQCDSFUncertainty(cat) );
  this->GetHistbyName(s_QCDest+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCDest+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1-this->getQCDSFUncertainty(cat) );
  this->GetHistbyName(s_QCD+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCD+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVUp+sub,strVar)->Scale( 1+this->getQCDSFUncertainty(cat) );
  this->GetHistbyName(s_QCD+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Add( this->GetHistbyName(s_QCDest+sub,strVar) );
  this->GetHistbyName(s_QCD+s_QCDSFUncert+channel+"_"+this->return2DString(cat)+s_13TeVDown+sub,strVar)->Scale( 1-this->getQCDSFUncertainty(cat) );

  if(jecShift){
    W_SF = ( this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName(s_WjecUp+"_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_"+s_WjecUp+sub,strVar)->Integral() : 0;
    W_rel = ( this->GetHistbyName("SS_Low_"+s_WjecUp+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_WjecUp+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_WjecUp+sub,strVar)->Integral() : 0;
    Z_rel = ( this->GetHistbyName("SS_Low_"+s_ZjecUp+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_ZjecUp+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_ZjecUp+sub,strVar)->Integral() : 0;
    VV_rel = ( this->GetHistbyName("SS_Low_"+s_VVjecUp+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_VVjecUp+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_VVjecUp+sub,strVar)->Integral() : 0;
    TT_rel = ( this->GetHistbyName("SS_Low_"+s_TTjecUp+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_TTjecUp+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_TTjecUp+sub,strVar)->Integral() : 0;
    EWKZ_rel = ( this->GetHistbyName("SS_Low_"+s_EWKZjecUp+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_EWKZjecUp+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecUp+sub,strVar)->Integral() : 0;

    this->GetHistbyName("AICR_"+s_WjecUp+sub,strVar)->Add( this->GetHistbyName( "SR_MC"+s_WjecUp+sub,strVar ) );
    this->GetHistbyName("AICR_"+s_WjecUp+sub,strVar)->Scale(W_SF);

    if(cat.Contains(s_antiiso)){
      W_rel=1;
      Z_rel=1;
      VV_rel=1;
      TT_rel=1;
      EWKZ_rel=1;
    }
    
    this->GetHistbyName("SS_Low_relaxed_"+s_WjecUp+sub,strVar)->Scale(  W_SF   );
    this->GetHistbyName("SS_Low_relaxed_"+s_WjecUp+sub,strVar)->Scale(  W_rel  );
    this->GetHistbyName("SS_Low_relaxed_"+s_VVjecUp+sub,strVar)->Scale( VV_rel );
    this->GetHistbyName("SS_Low_relaxed_"+s_ZjecUp+sub,strVar)->Scale(  Z_rel  );
    this->GetHistbyName("SS_Low_relaxed_"+s_TTjecUp+sub,strVar)->Scale( TT_rel );
    this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecUp+sub,strVar)->Scale( EWKZ_rel );

    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_"+s_data+sub,strVar)  );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_VVjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_ZjecUp+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_TTjecUp+sub,strVar), -1 );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_WjecUp+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecUp+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDjecUp+sub,strVar)->Add( this->GetHistbyName(s_QCDestjecUp+sub,strVar) );
    this->GetHistbyName(s_QCDestjecUp+sub,strVar)->Scale( this->QCD_OSSS(cat) );

    this->resetZeroBins(s_QCDestjecUp+sub,strVar);
    this->resetZeroBins(s_QCDjecUp+sub,strVar);
    this->resetZeroBins("AICR_"+s_WjecUp+sub,strVar);
    
    W_SF = ( this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName(s_WjecDown+"_OSW"+sub,strVar)->Integral() / this->GetHistbyName("OS_W_"+s_WjecDown+sub,strVar)->Integral() : 0;
    W_rel = ( this->GetHistbyName("SS_Low_"+s_WjecDown+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_WjecDown+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_WjecDown+sub,strVar)->Integral() : 0;
    Z_rel = ( this->GetHistbyName("SS_Low_"+s_ZjecDown+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_ZjecDown+sub,strVar)->Integral()  / this->GetHistbyName("SS_Low_relaxed_"+s_ZjecDown+sub,strVar)->Integral() : 0;
    VV_rel = ( this->GetHistbyName("SS_Low_"+s_VVjecDown+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_VVjecDown+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_VVjecDown+sub,strVar)->Integral() : 0;
    TT_rel = ( this->GetHistbyName("SS_Low_"+s_TTjecDown+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_TTjecDown+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_TTjecDown+sub,strVar)->Integral() : 0;
    EWKZ_rel = ( this->GetHistbyName("SS_Low_"+s_EWKZjecDown+sub,strVar)->Integral() > 0 )
    ? this->GetHistbyName("SS_Low_"+s_EWKZjecDown+sub,strVar)->Integral() / this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecDown+sub,strVar)->Integral() : 0;

    this->GetHistbyName("AICR_"+s_WjecDown+sub,strVar)->Add( this->GetHistbyName( "SR_MC"+s_WjecDown+sub,strVar ) );
    this->GetHistbyName("AICR_"+s_WjecDown+sub,strVar)->Scale(W_SF);

    if(cat.Contains(s_antiiso)){
      W_rel=1;
      Z_rel=1;
      VV_rel=1;
      TT_rel=1;
      EWKZ_rel=1;
    }
    
    this->GetHistbyName("SS_Low_relaxed_"+s_WjecDown+sub,strVar)->Scale(  W_SF   );
    this->GetHistbyName("SS_Low_relaxed_"+s_WjecDown+sub,strVar)->Scale(  W_rel  );
    this->GetHistbyName("SS_Low_relaxed_"+s_VVjecDown+sub,strVar)->Scale( VV_rel );
    this->GetHistbyName("SS_Low_relaxed_"+s_ZjecDown+sub,strVar)->Scale(  Z_rel  );
    this->GetHistbyName("SS_Low_relaxed_"+s_TTjecDown+sub,strVar)->Scale( TT_rel );
    this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecDown+sub,strVar)->Scale( EWKZ_rel );

    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_"+s_data+sub,strVar)  );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_VVjecDown+sub,strVar), -1 );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_ZjecDown+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_TTjecDown+sub,strVar), -1 );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_WjecDown+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Add( this->GetHistbyName("SS_Low_relaxed_"+s_EWKZjecDown+sub,strVar),  -1 );
    this->GetHistbyName(s_QCDjecDown+sub,strVar)->Add( this->GetHistbyName(s_QCDestjecDown+sub,strVar) );
    this->GetHistbyName(s_QCDestjecDown+sub,strVar)->Scale( this->QCD_OSSS(cat) );

    this->resetZeroBins(s_QCDestjecDown+sub,strVar);
    this->resetZeroBins(s_QCDjecDown+sub,strVar);
    this->resetZeroBins("AICR_"+s_WjecDown+sub,strVar);
    
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
  
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+sub,strVar),  -1 );
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar), -1 );

  if(jecShift){
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_WjecUp+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_TTjecUp+"_"+s_jetFakes+sub,strVar),  -1 );
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_ZjecUp+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_EWKZjecUp+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecUp+sub,strVar)->Add( this->GetHistbyName(s_VVjecUp+"_"+s_jetFakes+sub,strVar), -1 );

    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_WjecDown+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_TTjecDown+"_"+s_jetFakes+sub,strVar),  -1 );
    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_ZjecDown+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_EWKZjecDown+"_"+s_jetFakes+sub,strVar), -1 );
    this->GetHistbyName(s_jetFakesjecDown+sub,strVar)->Add( this->GetHistbyName(s_VVjecDown+"_"+s_jetFakes+sub,strVar), -1 );
  }

  this->GetHistbyName(s_jetFakes+"_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
  this->GetHistbyName(s_jetFakes+"_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+sub,strVar),  -1 );
  this->GetHistbyName(s_jetFakes+"_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_Up"+sub,strVar), -1 );  
  this->GetHistbyName(s_jetFakes+"_Zvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar), -1 );

  this->GetHistbyName(s_jetFakes+"_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
  this->GetHistbyName(s_jetFakes+"_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+sub,strVar),  -1 );
  this->GetHistbyName(s_jetFakes+"_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_Down"+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_Zvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar), -1 );

  this->GetHistbyName(s_jetFakes+"_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
  this->GetHistbyName(s_jetFakes+"_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_Up"+sub,strVar),  -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systUp"+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar), -1 );

  this->GetHistbyName(s_jetFakes+"_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)   );
  this->GetHistbyName(s_jetFakes+"_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_Down"+sub,strVar),  -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+sub,strVar), -1 );
  this->GetHistbyName(s_jetFakes+"_TTvar_systDown"+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar), -1 );

  this->resetZeroBins(s_jetFakes+sub,strVar);
  this->resetZeroBins(s_jetFakes+"_Zvar_systUp"+sub,strVar);
  this->resetZeroBins(s_jetFakes+"_Zvar_systDown"+sub,strVar);
  this->resetZeroBins(s_jetFakes+"_TTvar_systUp"+sub,strVar);
  this->resetZeroBins(s_jetFakes+"_TTvar_systDown"+sub,strVar);
  if(jecShift){
    this->resetZeroBins(s_jetFakesjecUp+sub,strVar);
    this->resetZeroBins(s_jetFakesjecDown+sub,strVar);
  }

  double ratio = 0;
  for( auto syst : FFsyst[channel] ){
    TString tmp=syst;
    this->getCorrectUncertaintyString( tmp );
    this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+tmp+sub,strVar)   );
    if( this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) > 0
        && this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) < 10e10 )       this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar), -1 );
    if( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) > 0
        && this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) < 10e10 )     this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar),  -1 );
    if( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) > 0
        && this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) < 10e10 )       this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar), -1 );
    if( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) > 0
        && this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) < 10e10 ) this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar), -1 );
    if( this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) > 0
        && this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar)->Integral(0,this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar)->GetNbinsX()+1 ) < 10e10 )     this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar), -1 );

    this->resetZeroBins(s_jetFakes+"_"+tmp+sub,strVar);
    
    this->GetHistbyName(s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      cout << ratio << endl;
      if( ratio < 1) normUp_jetFakes_syst = TMath::Sqrt( TMath::Power(normUp_jetFakes_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_jetFakes_syst = TMath::Sqrt( TMath::Power(normDown_jetFakes_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_jetFakes = TMath::Sqrt( TMath::Power(normUp_jetFakes,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_jetFakes = TMath::Sqrt( TMath::Power(normDown_jetFakes,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_data+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_data_syst = TMath::Sqrt( TMath::Power(normUp_data_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_data_syst = TMath::Sqrt( TMath::Power(normDown_data_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_data = TMath::Sqrt( TMath::Power(normUp_data,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_data = TMath::Sqrt( TMath::Power(normDown_data,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_W+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_W_syst = TMath::Sqrt( TMath::Power(normUp_W_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_W_syst = TMath::Sqrt( TMath::Power(normDown_W_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_W = TMath::Sqrt( TMath::Power(normUp_W,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_W = TMath::Sqrt( TMath::Power(normDown_W,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_TT+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_TT+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_TT_syst = TMath::Sqrt( TMath::Power(normUp_TT_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_TT_syst = TMath::Sqrt( TMath::Power(normDown_TT_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_TT = TMath::Sqrt( TMath::Power(normUp_TT,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_TT = TMath::Sqrt( TMath::Power(normDown_TT,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_Z+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_Z+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_Z_syst = TMath::Sqrt( TMath::Power(normUp_Z_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_Z_syst = TMath::Sqrt( TMath::Power(normDown_Z_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_Z = TMath::Sqrt( TMath::Power(normUp_Z,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_Z = TMath::Sqrt( TMath::Power(normDown_Z,2) + TMath::Power( 1-ratio,2 ) );
    }

    this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_EWKZ_syst = TMath::Sqrt( TMath::Power(normUp_EWKZ_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_EWKZ_syst = TMath::Sqrt( TMath::Power(normDown_EWKZ_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_EWKZ = TMath::Sqrt( TMath::Power(normUp_EWKZ,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_EWKZ = TMath::Sqrt( TMath::Power(normDown_EWKZ,2) + TMath::Power( 1-ratio,2 ) );
    }
    
    this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+tmp+sub,strVar)  );
    ratio = ( this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar)->Integral( 0, this->GetHistbyName(s_VV+"_"+s_jetFakes+sub,strVar)->GetNbinsX()+1  ) )/( this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Integral( 0, this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->GetNbinsX()+1  ) );
    this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+"_"+tmp+sub,strVar)->Scale( ratio );
    if(tmp.Contains("syst")){
      if( ratio < 1) normUp_VV_syst = TMath::Sqrt( TMath::Power(normUp_VV_syst,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_VV_syst = TMath::Sqrt( TMath::Power(normDown_VV_syst,2) + TMath::Power( 1-ratio,2 ) );
    } else{
      if( ratio < 1) normUp_VV = TMath::Sqrt( TMath::Power(normUp_VV,2) + TMath::Power( 1-ratio,2 ) );
      else normDown_VV = TMath::Sqrt( TMath::Power(normDown_VV,2) + TMath::Power( 1-ratio,2 ) );
    }
  }
  this->GetHistbyName(s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+sub,s_norm);
  this->GetHistbyName(s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_jetFakes);
  this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_data);
  this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_W);
  this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_TT);
  this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_Z);
  this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_EWKZ);
  this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(1,normUp_VV);
  this->GetHistbyName(s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_jetFakes);
  this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_data);
  this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_W);
  this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_TT);
  this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_Z);
  this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_EWKZ);
  this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(2,normDown_VV);
  this->GetHistbyName(s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_jetFakes_syst);
  this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_data_syst);
  this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_W_syst);
  this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_TT_syst);
  this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_Z_syst);
  this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_EWKZ_syst);
  this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(3,normUp_VV_syst);
  this->GetHistbyName(s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_jetFakes_syst);
  this->GetHistbyName(s_data+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_data_syst);
  this->GetHistbyName(s_W+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_W_syst);
  this->GetHistbyName(s_TT+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_TT_syst);
  this->GetHistbyName(s_Z+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_Z_syst);
  this->GetHistbyName(s_EWKZ+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_EWKZ_syst);
  this->GetHistbyName(s_VV+"_"+s_jetFakes+"_"+s_norm+sub)->SetBinContent(4,normDown_VV_syst);

  if(channel == "tt"){
    this->GetHistbyName(s_FFrest+sub,strVar)->Add( this->GetHistbyName(s_ZJ+"_"+s_rest+sub,strVar)   );
    this->GetHistbyName(s_FFrest+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_rest+sub,strVar)   );
    this->GetHistbyName(s_FFrest+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_rest+sub,strVar)   );
    this->GetHistbyName(s_FFrest+sub,strVar)->Add( this->GetHistbyName(s_TTJ+"_"+s_rest+sub,strVar)   );
    this->GetHistbyName(s_FFrest+sub,strVar)->Add( this->GetHistbyName(s_VVJ+"_"+s_rest+sub,strVar)   );

    this->GetHistbyName(s_ZJ+"_"+s_rest+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_rest+sub,strVar)   );

    if(jecShift){
      this->GetHistbyName(s_FFrestjecUp+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_restjecUp+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecUp+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_restjecUp+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecUp+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_restjecUp+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecUp+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_restjecUp+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecUp+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_restjecUp+sub,strVar)   );

      this->GetHistbyName(s_ZJ+"_"+s_restjecUp+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_restjecUp+sub,strVar)   );

      this->GetHistbyName(s_FFrestjecDown+sub,strVar)->Add( this->GetHistbyName(s_Z+"_"+s_restjecDown+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecDown+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_restjecDown+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecDown+sub,strVar)->Add( this->GetHistbyName(s_W+"_"+s_restjecDown+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecDown+sub,strVar)->Add( this->GetHistbyName(s_TT+"_"+s_restjecDown+sub,strVar)   );
      this->GetHistbyName(s_FFrestjecDown+sub,strVar)->Add( this->GetHistbyName(s_VV+"_"+s_restjecDown+sub,strVar)   );

      this->GetHistbyName(s_ZJ+"_"+s_restjecDown+sub,strVar)->Add( this->GetHistbyName(s_EWKZ+"_"+s_restjecDown+sub,strVar)   );
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
  if(do2DFit) D2+="-2D";
  if(doMC) D2+="-MCsum";

  
  for(auto var : vars){
    outfile_name << "histos/"  << "htt_" << channel << ".inputs-sm-13TeV-"<<var<<D2<<".root";
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
                  || histo_names.at(i).Contains("SS_"+s_W)
                  //|| histo_names.at(i).Contains("OS_"+s_W)
                  || histo_names.at(i).Contains(s_W+"_OSW")
                  || histo_names.at(i).Contains(s_WjecUp+"_OSW")
                  || histo_names.at(i).Contains(s_WjecDown+"_OSW")
                  )
             ) continue;
          if(!keepFFDebugHistos
             && ( histo_names.at(i).Contains(s_data+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_W+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_WjecUp+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_WjecDown+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_Z+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_ZjecUp+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_ZjecDown+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_TT+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_TTjecUp+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_TTjecDown+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_VV+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_VVjecUp+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_VVjecDown+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_EWKZ+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_EWKZjecUp+"_"+s_jetFakes)
                  || histo_names.at(i).Contains(s_EWKZjecDown+"_"+s_jetFakes)                  
                  )
             ) continue;
          if(!keepZGenJetsSplitting
             && ( histo_names.at(i).Contains("Z_0Jets")
                  || histo_names.at(i).Contains("Z_1Jets")
                  || histo_names.at(i).Contains("Z_ge2Jets")
                  )
             ) continue;
          tmp = histo_names.at(i);
          if( (cat.Contains(s_wjets)
               ||cat.Contains(s_antiiso) )
              && tmp.Contains(s_West) ) continue;
          if( (cat.Contains(s_wjets)
               ||cat.Contains(s_antiiso) )
              && tmp.Contains(s_QCDest) ) continue;
          if( tmp.Contains("SR_MC") ) continue;
          if( ( cat.Contains(s_wjets)
                || cat.Contains(s_antiiso) )
              && tmp.Contains("SR_"+s_W) ) continue;
          else if( !cat.Contains(s_wjets)
                   && !cat.Contains(s_antiiso)
                   && tmp.Contains("SR_"+s_W ) ) tmp.ReplaceAll("SR_","");
          else if( cat.Contains(s_wjets)
                   && tmp.Contains("WCR_"+s_W ) ) tmp.ReplaceAll("WCR_","");
          else if( cat.Contains(s_antiiso)
                   && tmp.Contains("AICR_"+s_W ) ) tmp.ReplaceAll("AICR_","");
          else if( tmp.Contains("WCR_"+s_W)
                   || tmp.Contains("AICR_"+s_W) ) continue;
          if( cat.Contains(s_wjets)
              && tmp.Contains("OS_W_"+s_QCD ) ) tmp.ReplaceAll("OS_W_","");
          else if( cat.Contains(s_wjets)
                   && !tmp.Contains(s_QCDSFUncert)
                   && tmp.Contains(s_QCD) ) continue;

          if(tmp.Contains("OS_"+s_W)) continue;
          tmp.ReplaceAll(sub, "");
          tmp.ReplaceAll(s_jecUp,s_CMSjecScale+s_13TeVUp);
          tmp.ReplaceAll(s_jecDown,s_CMSjecScale+s_13TeVDown);
          tmp.ReplaceAll(s_jetToTauFakeUp,s_CMSjetToTauFake+s_13TeVUp);
          tmp.ReplaceAll(s_jetToTauFakeDown,s_CMSjetToTauFake+s_13TeVDown);
          histos.at(i)->SetName(tmp);
          if(do2DFit ){
            if( is2DCategories(cat) ){
              tmp.ReplaceAll("2D","");
              this->resetZeroBins(histos.at(i));
              histos.at(i)->Write(tmp, TObject::kWriteDelete);
            }
            else if( cat == s_inclusive){
              this->resetZeroBins(histos.at(i));
              histos.at(i)->Write(tmp, TObject::kWriteDelete);
            }
            else continue;
          }
          else {
            this->resetZeroBins(histos.at(i));
            histos.at(i)->Write(tmp, TObject::kWriteDelete);
          }
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
