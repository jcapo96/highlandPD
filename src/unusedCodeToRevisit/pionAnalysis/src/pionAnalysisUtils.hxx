#ifndef pionAnalysisUtils_h
#define pionAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace pionAnaUtils{

  // Add custom categories
  void AddCustomCategories();
  
  // Fill Custom categories
  void FillCustomCategories(AnaEventB* event, AnaParticle* part, PDCounters& counters);
  void FillDaupionanaCategory(AnaEventB* event, AnaTrueParticlePD* part, AnaTrueParticlePD*  trueBeamPart);
  
}

#endif
