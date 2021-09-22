#ifndef kaonAnalysisUtils_h
#define kaonAnalysisUtils_h

#include "pdDataClasses.hxx"

namespace kaonAnaUtils{

  // Add custom categories
  void AddCustomCategories();
  
  // Fill Custom categories
  void FillCustomCategories(AnaEventB* event, AnaParticle* part);
 
  void FillDaughterCategory(AnaEventB* event, AnaParticlePD* part);
  void FillDaughterTSCategory(AnaEventB* event, AnaParticlePD* part1, AnaParticlePD* part2);
 
  void FillGDaughterCategory(AnaEventB* event, AnaParticlePD* part, Int_t indx1, Int_t indx2);
  void FillGDaughterKaonCategory(std::string catname, AnaEventB* event, AnaParticlePD* part1, AnaParticlePD* part2, Int_t indx1, Int_t indx2);

  void FillGGDaughterCategory(AnaEventB* event, AnaParticlePD* part, Int_t indx1, Int_t indx2, Int_t indx3);
  void FillGGDaughterKaonCategory(AnaEventB* event, AnaParticlePD* part1, AnaParticlePD* part2, AnaParticlePD* part3, Int_t indx1, Int_t indx2, Int_t indx3);

  std::pair<Int_t, Int_t> GetKaonDecayMode (AnaTrueParticleB* trueKaon, std::vector<AnaTrueParticleB*> trueParts);
  std::pair<bool,  Int_t> MuonFromKaonChain(AnaTrueParticleB* trueKaon, std::vector<AnaTrueParticleB*> trueParts);
  
}

#endif
