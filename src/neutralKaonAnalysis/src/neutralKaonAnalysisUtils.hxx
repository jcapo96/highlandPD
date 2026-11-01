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

  // Signal definition used for categories and event-display filtering: the two reco particles of the decay vertex are
  // the pi+ and pi- of one true K0S -> pi+ pi- decay, whatever the reco parent (the definition of the Lambda signal,
  // from the daughters only). A K0S missing from the stored truth (kept only when it descends from a K+) is
  // identified by the true pi+ pi- mass of a pair with a common true start.
  bool IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // The true K0S of a signal candidate when it is stored in the truth; otherwise nullptr (also for signal whose K0S
  // was not stored).
  AnaTrueParticlePD* GetSignalTrueParent(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Signal (IsSignalCandidate) whose neutral is also truth-matched through its reco parent: the true K0S is a daughter
  // of the reco parent's true particle (the signal definition until 2026-11-01; tree variable k0parenttruematch).
  bool IsParentMatchedSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  // Fill Custom categories
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillVertexCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);

  /// Same integer codes as category "signal" in FillSignalCandidateCategory, using only the annihilation vertex and
  /// event (truth-matched evaluation). For vertex creation before AnaNeutralParticlePD exists.
  int GetSignalCategoryCodeForAnnihilationVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
}

#endif
