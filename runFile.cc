#include "interface/CreateHistos.h"

#include <iostream>
#include <string>

void runFile() {

	
    CreateHistos *Analyzer = new CreateHistos();

    Analyzer->run();

    delete Analyzer;
  

}
#ifndef __CINT__
int main(int argc, char* argv[]) {

  runFile();    
}
#endif