#ifndef __FFCalculator__
#define __FFCalculator__

#include "interface/GlobalClass.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTGraph.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH2F.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH3D.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTFormula.h"
#include "HTTutilities/Jet2TauFakes/interface/IFunctionWrapper.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"

class FFCalculator : public GlobalClass{
 public:
  FFCalculator();
  ~FFCalculator();
  
  void doUpDownReplace( TString &replaceString );
  void doStatUncertaintyReplace( TString &replaceString );
  void doSystUncertaintyReplace( TString &replaceString );
  void getCorrectUncertaintyString( TString &replaceString );

  void getFFInputs(vector<double>&inputs);
  void getFF1Inputs(vector<double>&inputs);
  void getFF2Inputs(vector<double>&inputs);
  void initFakeFactors();
  void applyFF(float var, float weight, TString cat, TString strVar, TString fname, int isData, TString extend="");

  vector<Double_t> FFinputs;
  map< TString, TFile*> FFfile;
  map< TString, FakeFactor*> FFObj;
  map< TString, vector<string> > FFsyst;
  
};



#endif
