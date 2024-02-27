#ifndef kaonAnalysisUtils_h
#define kaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace kaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();

  void AddBeamParticleReducedCategory();
  void AddBeamParticleOriginCategory();
  void AddBestCandidateParticleReducedCategory();
  void AddCandidateParticleReducedCategory();
  void AddCandidateDaughterParticleReducedCategory();
  void AddCandidateDaughterMuonCategory();
  void AddCandidateDaughterMuonReducedCategory();
  
  // Fill Custom categories
  void FillBeamParticleReducedCategory(AnaParticlePD* beampart);
  void FillBeamParticleOriginCategory(AnaParticlePD* beampart);
  void FillBestCandidateParticleReducedCategory(AnaParticlePD* part);
  void FillCandidateParticleReducedCategory(AnaParticlePD* part);
  void FillCandidateDaughterParticleReducedCategory(AnaParticlePD* part);
  void FillCandidateDaughterMuonCategory(AnaParticlePD* parent, AnaParticlePD* daughter);
  void FillCandidateDaughterMuonReducedCategory(AnaParticlePD* parent, AnaParticlePD* daughter);
}

#endif
