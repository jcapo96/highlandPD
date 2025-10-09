#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddNeutralParticleSignalBackgroundCategory();
  void AddNeutralParticlePDGCategory();
  void AddNeutralParticleChargeCategory();
  void AddNeutralParticleNoTruthCategory();
  void AddSignalCandidateCategory();

  // Fill Custom categories
  void FillNeutralParticleSignalBackgroundCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillNeutralParticlePDGCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillNeutralParticleChargeCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillNeutralParticleNoTruthCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle);
}

#endif
