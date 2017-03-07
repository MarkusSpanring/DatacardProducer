#include "interface/NameStrings.h"

using namespace std;

const TString channel="mt";
const TString version="v10";
const int doMC=0;
const double usedLuminosity=35.9; 
const TString doSvfit="woSVFIT";
const TString reduced="";//"";//"_reduced";
const TString FFversion="fakeFactors_20170228.root";
const vector<TString> variables={s_mvis}; //s_mvis,s_met,s_pt1,s_pt2
const vector<TString> categories={s_inclusive,s_0jet,s_wjets_0jet_cr,s_antiiso_0jet_cr,s_boosted,s_wjets_boosted_cr,s_antiiso_boosted_cr,s_vbf,s_wjets_vbf_cr,s_antiiso_vbf_cr};
///////////////////////////////////////////////////////////////////////////
const int applyMTCut=1;//MTcut in inclusive selection, default==1 for mt&et, default==0 for tt
const int resetZero=1;
const int do2DFit=1;
const int ptShift=1;
const int jecShift=1;
const int doInitialize=1;
const int useMVAMET=0;
const int calcFF=1;
const int keepDebugHistos=0;
const int keepFFDebugHistos=0;
const int keepZGenJetsSplitting=0;

