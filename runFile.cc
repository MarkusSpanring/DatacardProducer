#include "interface/CreateHistos.h"

#include <iostream>
#include <string>

void runFile(TString isTest) {

	
    CreateHistos *Analyzer = new CreateHistos();

    Analyzer->run(isTest);

    delete Analyzer;
  

}
#ifndef __CINT__
int main(int argc, char* argv[]) {
  TString isTest = argc > 1 ? "test" : ""; 
  runFile(isTest);    
}
#endif
