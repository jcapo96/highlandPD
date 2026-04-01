#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddSignalCandidateCategory();

  // Legit annihilation vertex: both reco daughters share the same true parent.
  bool IsLegitVertexCandidate(AnaNeutralParticlePD* neutralParticle);

  // Legit-vertex subtype tags by true-parent daughter multiplicity.
  bool IsLegitVertexFromTwoBodyDecay(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  bool IsLegitVertexFromMultiBodyDecay(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Signal definition used for categories and event-display filtering.
  bool IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // When IsSignalCandidate is satisfied, returns the true K0 parent; otherwise nullptr.
  AnaTrueParticlePD* GetSignalTrueParent(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
}

#endif
