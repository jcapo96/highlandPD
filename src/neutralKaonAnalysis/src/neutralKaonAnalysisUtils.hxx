#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddSignalCandidateCategory();
  void AddLooseSignalCandidateCategory();

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillLooseSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
}

#endif
