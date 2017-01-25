#include "interface/NameStrings.h"

using namespace std;

const TString channel="et";
const TString version="v6";
const double usedLuminosity=27.9; //20.2 wo runG //27.9 
const TString doSvfit="woSVFIT";
const TString reduced="";//"";//"_reduced";
const TString FFversion="fakeFactors_20170111.root";
const vector<TString> variables={s_mvis};
const vector<TString> categories={s_inclusive,s_0jet,s_boosted,s_vbf};
///////////////////////////////////////////////////////////////////////////
const int applyMTCut=1; //MTcut in inclusive selection, default==1 for mt&et, default==0 for tt
const int resetZero=1;
const int do2DFit=1;
const int ptShift=1;
const int jecShift=1;
const int useMVAMET=0;
const int calcFF=1;
const int keepDebugHistos=0;
const int keepFFDebugHistos=0;
const int keepZGenJetsSplitting=0;

