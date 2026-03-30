#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddSignalCandidateCategory();

  // Signal definition used for categories and event-display filtering.
  bool IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
}

#endif
