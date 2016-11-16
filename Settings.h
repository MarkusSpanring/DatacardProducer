const TString channel="mt";
const TString version="v3";
const double usedLuminosity=12.9;
const int ptShift=1;
const TString svfit="woSVFIT";
const int useMVAMET=0;
const int calcFF=1;
const int applyMTCut=1; //MTcut in inclusive selection, default==1
<<<<<<< HEAD
const vector<TString> variables={"m_vis"}; //{"m_vis","mt_1","jpt_1","jpt_2","jeta_1","jeta_2","jdeta","mjj","jeta1eta2"};
const vector<TString> categories={"inclusive"}; //{"0Jet_low,0Jet_high,1Jet_low,1Jet_high,VBF_low,VBF_high"};
=======
const vector<TString> variables={"m_vis","pt_2", "jeta_1", "jeta_2", "jdeta","mjj"}; //{"m_vis","mt_1","jpt_1","jpt_2","jeta_1","jeta_2","jdeta","mjj","jeta1eta2"};
const vector<TString> categories={"inclusive","1Jet_low","1Jet_high","0jet_low","0Jet_high","VBF_low","VBF_high"};
>>>>>>> datacard-remote/master
