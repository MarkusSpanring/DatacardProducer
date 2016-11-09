const TString channel="mt";
const TString version="v3";
const double usedLuminosity=12.9;
const int ptShift=1;
const TString svfit="woSVFIT";
const int useMVAMET=0;
const int applyMTCut=1; //MTcut in inclusive selection, default==1
const vector<TString> variables={"m_vis"}; //{"m_vis","mt_1","jpt_1","jpt_2","jeta_1","jeta_2","jdeta","mjj","jeta1eta2"};
const vector<TString> categories={"inclusive","Jet0_low","Jet0_high","Jet1_low","Jet1_high","VBF_low","VBF_high"};
