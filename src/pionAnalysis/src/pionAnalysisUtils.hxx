#ifndef pionAnalysisUtils_h
#define pionAnalysisUtils_h

#include "PionAnaDataClasses.hxx"

namespace pionAnaUtils{

  // Add custom categories
  void AddCustomCategories();
  
  // Fill Custom categories
  void FillCustomCategories(AnaEventB* event, AnaParticle* part, PionAnaCounters& counters);
  void FillDaupionanaCategory(AnaEventB* event, AnaTrueParticlePionAna* part, AnaTrueParticlePionAna*  trueBeamPart);
  
}

#endif
