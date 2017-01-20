#ifndef __FFCalculator__
#define __FFCalculator__

#include "interface/GlobalClass.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTGraph.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH2F.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH3D.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTFormula.h"
#include "HTTutilities/Jet2TauFakes/interface/IFunctionWrapper.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"

class FFCalculator : public GlobalClass{
 public:
  FFCalculator();
  ~FFCalculator();
  
};



#endif
