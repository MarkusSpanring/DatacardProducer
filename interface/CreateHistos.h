#ifndef __CREATEHISTOS__
#define __CREATEHISTOS__

#include "interface/SelectionAnalyzer.h"

class CreateHistos : public SelectionAnalyzer{
 public:

  CreateHistos();
  ~CreateHistos();

  float getAntiLep_tauscaling();

  int isZFile(TString fileName);
  int isEWKZFile(TString fileName);
  int isSignalFile(TString fileName);
  int isTTFile(TString fileName);
  int isWFile(TString fileName);
  int isVVFile(TString fileName); 

  void CreateQCD_osw(TString strVar, TString cat, TString extend="");
  void CreateW(TString strVar, TString cat, TString extend="");
  void CreateQCD(TString strVar, TString cat, TString extend="");
  void Estimate_W_QCD(TString strVar, TString cat, TString extend="");
  void EstimateFF(TString strVar, TString cat, TString extend="");

  void loadFile(TString filename);
  void run(TString isTest);
  void clearHistos();
  void writeHistos(TString channel, vector<TString> cats, vector<TString> vars);
  
  TFile *outfile;

  
};



#endif
