CXX=g++  -std=c++11
INCDIR=.
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

all: CreateHistos.o ntuple.o runFile

%.o: src/%.cc
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -fpic -c $<

CreateHistos.o: src/CreateHistos.cc
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c src/CreateHistos.cc

ntuple.o: ntuple.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c ntuple.C

runFile: runFile.cc
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -o $@ ntuple.o CreateHistos.o ../tmp/slc6_amd64_gcc530/src/HTTutilities/Jet2TauFakes/src/HTTutilitiesJet2TauFakes/libHTTutilitiesJet2TauFakes.so runFile.cc


clean:
	rm *.o runFile 
