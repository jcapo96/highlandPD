#ifndef pdUtilsRangePID_h
#define pdUtilsRangePID_h

#include "pdDataClasses.hxx"
#include <utility>

namespace pdAnaUtils {

  Float_t ComputeRangeMomentum(double trkrange, int pdg);
  Float_t ComputeCSDARange(double beammom, int pdg);
  Float_t ComputeKineticEnergy(const AnaParticlePD &part);

  std::pair<double, int> Chi2PID(const AnaParticlePD& part, const int pdg);
  std::pair<double, int> Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR);

  /// Chi2PID first / second (mean χ² per contributing hit); -999 if null or Chi2PID invalid. Same as neutral-kaon truth tree Reco*Chi2Ndf.
  Float_t Chi2PIDChi2PerHit(const AnaParticlePD* part, int pdg);

  /// Pion (211): same per-hit χ² as Chi2PID (dedx_range_pi template + bin error + dedx_res), but only collection hits with
  /// 0 < RR < maxResidualRangeCm and the same skip / optional dE/dx window as ComputePionBraggWindowChi2PiEq61.
  bool ComputePionBraggWindowChi2PionRangeTemplate(AnaParticlePD* part, double maxResidualRangeCm, int minHits,
                                                   int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                                   double dedxMaxMeVcm, double& meanChi2PerHit, int& nHitsUsed);

}

#endif
