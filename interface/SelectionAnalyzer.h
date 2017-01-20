#ifndef __SELECTIONANALYZER__
#define __SELECTIONANALYZER__

#include "interface/FFCalculator.h"

class SelectionAnalyzer : public FFCalculator{
 public:
  SelectionAnalyzer();
  ~SelectionAnalyzer();
  
  void DYSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void EWKZSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void VVSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void TSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void WSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void dataSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  void signalSelections(float var, float weight, TString cat, TString strVar, TString fname, TString extend="");
  

};



#endif
