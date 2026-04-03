#ifndef pdUtilsTopology_h
#define pdUtilsTopology_h

#include "pdDataClasses.hxx"
#include <vector>

namespace pdAnaUtils {

  void ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance);

  Float_t ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);

  Float_t ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);

  Double_t ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts);

}

#endif
