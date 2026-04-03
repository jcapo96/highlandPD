#ifndef pdUtilsTruthBeam_h
#define pdUtilsTruthBeam_h

#include "pdDataClasses.hxx"
#include <vector>

namespace pdAnaUtils {

  AnaTrueParticle* GetBeamTrueParticle(const AnaSpillB& spill);

  AnaTrueParticlePD* GetTrueParticle(AnaEventB* event, Int_t ID);
  AnaTrueParticlePD* GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID);

  AnaParticlePD* GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID);

  bool isBeamLike(AnaParticlePD* part, AnaBeamPD* beam);

  AnaParticlePD* GetBeamParticle(const AnaEventC& event);
  AnaTrueParticlePD* GetTrueBeamParticle(const AnaEventC& event);

}

#endif
