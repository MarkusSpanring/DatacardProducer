#ifndef __CREATEHISTOS__
#define __CREATEHISTOS__

#include "interface/SelectionAnalyzer.h"

class CreateHistos : public SelectionAnalyzer{
 public:

  vector<TString> cats;
  vector<TString> vars;
  vector< vector<TString> > files; 

  CreateHistos();
  ~CreateHistos();


  void initFakeFactors();
  int is1DCategories(TString category);
  int is2DCategories(TString category);
  void DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void VVSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void signalSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend="");

  float getAntiLep_tauscaling();
  void getFFInputs(vector<double>&inputs);
  void getFF1Inputs(vector<double>&inputs);
  void getFF2Inputs(vector<double>&inputs);

  double QCD_OSSS(TString cat);
  void CreateQCD_osw(TString strVar, TString cat, TString extend="");
  void CreateW(TString strVar, TString cat, TString extend="");
  void CreateQCD(TString strVar, TString cat, TString extend="");
  void Estimate_W_QCD(TString strVar, TString cat, TString extend="");
  void EstimateFF(TString strVar, TString cat, TString extend="");


  int OS_W(TString cat);
  int SS_W(TString cat);
  int relaxed_W(TString cat, TString mt);
  int SS_Low(TString cat);
  int SS_Low_relaxed(TString cat);

  void loadFile(TString filename);
  void run(TString isTest);
  void clearHistos();
  void writeHistos(TString channel, vector<TString> cats, vector<TString> vars);
  
  double get2DVar(TString sub);
  int returnBin(vector<double> bins, double value);

  TFile *outfile;

  vector<Double_t> FFinputs;
  map< TString, TFile*> FFfile;
  map< TString, FakeFactor*> FFObj;
  map< TString, vector<string> > FFsyst;

  
};



#endif
