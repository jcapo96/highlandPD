#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddSignalCandidateCategory();

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
}

#endif
