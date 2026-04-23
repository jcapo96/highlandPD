#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddSignalCandidateCategory();
  void AddVertexCandidateCategory();

  // Legit annihilation vertex: both reco daughters share the same true parent.
  bool IsLegitVertexCandidate(AnaNeutralParticlePD* neutralParticle);

  // Signal definition used for categories and event-display filtering.
  bool IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // When IsSignalCandidate is satisfied, returns the true K0 parent; otherwise nullptr.
  AnaTrueParticlePD* GetSignalTrueParent(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillVertexCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  /// Same integer codes as category "signal" in FillSignalCandidateCategory, using only the annihilation vertex and
  /// event (truth-matched evaluation). For vertex creation before AnaNeutralParticlePD exists.
  int GetSignalCategoryCodeForAnnihilationVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
}

#endif
