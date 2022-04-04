#ifndef kaonAnalysisUtils_h
#define kaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace kaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();
  
  void AddCandidateDaughterMuonCategory();
  
  // Fill Custom categories
  void FillCandidateDaughterMuonCategory(AnaParticlePD* parent, AnaParticlePD* daughter);
}

#endif
