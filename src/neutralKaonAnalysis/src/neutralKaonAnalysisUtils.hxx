#ifndef neutralKaonAnalysisUtils_h
#define neutralKaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace neutralKaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddNeutralParticleSignalBackgroundCategory();

  // Fill Custom categories
  void FillNeutralParticleSignalBackgroundCategory(AnaNeutralParticlePD* neutralParticle);
}

#endif
