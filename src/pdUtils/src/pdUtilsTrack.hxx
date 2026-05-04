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

  /// Strict-inequality TPC box used by SpaceCharge::IsInsideBoundaries
  /// (|X|>0 && |X|<360 cm; Y>5.2 && Y<604 cm; Z>-0.5 && Z<695.3 cm).
  /// Keep in sync with SpaceCharge::IsInsideBoundaries (SpaceCharge.cxx).
  bool IsInTPCSCEBox(const TVector3& p);

  void EstimateHitsDirection(AnaParticlePD* part);

  void ComputeResidualRange(AnaParticlePD* part);

  Double_t ComputeDepositedEnergy(AnaParticlePD* part);

  Double_t EstimateTrueMomAtAPABorder(AnaParticlePD* part);

}

#endif
