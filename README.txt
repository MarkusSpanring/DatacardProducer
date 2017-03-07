cmsrel CMSSW_XXX
git clone https://github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes
mkdir HTTutilities/Jet2TauFakes/data/
cp -r /afs/cern.ch/user/j/jbrandst/public/Htautau/FakeRate/2016/20170228/* HTTutilities/Jet2TauFakes/data/
scram b -j 8
git clone git@github.com:MarkusSpanring/DatacardProducer.git
