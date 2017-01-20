const TString channel="et";
const TString version="v5";
const double usedLuminosity=27.9; //20.2 wo runG //27.9 
const int ptShift=0;
const int jecShift=0;
const TString doSvfit="woSVFIT";
const TString reduced="";//"";//"_reduced";
const int useMVAMET=0;
const int calcFF=1;
const int keepDebugHistos=0;
const int keepFFDebugHistos=0;
const int keepZGenJetsSplitting=0;
const TString FFversion="fakeFactors_20170111.root";
const int applyMTCut=1; //MTcut in inclusive selection, default==1 for mt&et, default==0 for tt
const int resetZero=1;
const int do2DFit=1;
const int doOfficialNaming=0;
//const vector<TString> variables={"m_sv","pt_sv","m_vis","pt_1","pt_2","mttot","met","Hpt","mjj"};
const vector<TString> variables={"m_vis"};
//const vector<TString> categories={"inclusive","0Jet_low","0Jet_high","1Jet_low","1Jet_high","VBF_low","VBF_high"};
const vector<TString> categories={"inclusive","0jet","boosted","vbf"};


