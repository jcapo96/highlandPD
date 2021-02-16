#ifndef kaonAnalysisUtils_h
#define kaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace kaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();
  
  // Fill Custom categories
  void FillCustomCategories(AnaEventB* event, AnaParticle* part);
  void FillDaughterCategory(AnaEventB* event, AnaTrueParticlePD* part);
  
}

#endif
