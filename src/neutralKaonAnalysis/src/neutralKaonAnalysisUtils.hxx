#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddNeutralParticleSignalBackgroundCategory();
  void AddNeutralParticleDetailedBackgroundCategory();

  // Fill Custom categories
  void FillNeutralParticleSignalBackgroundCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
  void FillNeutralParticleDetailedBackgroundCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event);
}

#endif
