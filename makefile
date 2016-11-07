CXX=g++  -std=c++11
INCDIR=.
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

all: CreateHistos.o ntuple.o runFile 

CreateHistos.o: src/CreateHistos.cc
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c src/CreateHistos.cc

ntuple.o: ntuple.C
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -c ntuple.C

runFile: runFile.cc
	$(CXX) -I$(INCDIR) -I$(ROOTINC) $(ROOTLIB) -o $@ ntuple.o CreateHistos.o runFile.cc


clean:
	rm *.o runFile 
