#ifndef pdUtilsTrack_h
#define pdUtilsTrack_h

#include "pdDataClasses.hxx"
#include <vector>

namespace pdAnaUtils {

  Float_t ComputeTrackLengthFromHitPosition(const AnaParticlePD* part);
  Float_t ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part);

  void ComputeParticlePositionAndDirection(AnaParticlePD* part);

  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> dEdx);
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD> hits);

  Float_t ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange = 9999);

  bool IsStoppingInFV(AnaParticlePD *part);

  int GetHitTPCid(AnaHitPD& hit);
  int GetPosTPCid(TVector3 pos);

  void EstimateHitsDirection(AnaParticlePD* part);

  void ComputeResidualRange(AnaParticlePD* part);

  Double_t ComputeDepositedEnergy(AnaParticlePD* part);

  Double_t EstimateTrueMomAtAPABorder(AnaParticlePD* part);

}

#endif
