#ifndef pdUtilsEventKinematics_h
#define pdUtilsEventKinematics_h

#include "pdDataClasses.hxx"

namespace pdAnaUtils {

  Float_t ComputeCalorimetricMomentum(AnaParticlePD* part, int pdg, bool includeDecayProducts = false);

  void GetBeamQualityCuts(AnaEventPD* event,
                          double &mean_x, double &mean_y, double &mean_z,
                          double &sigma_x, double &sigma_y, double &sigma_z,
                          double &cos);

  Float_t ComputeTrueInvariantMass(const AnaTrueParticlePD& part1, const AnaTrueParticlePD& part2,
                                   Float_t mass1, Float_t mass2);

}

#endif
