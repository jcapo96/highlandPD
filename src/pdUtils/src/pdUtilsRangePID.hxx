#ifndef pdUtilsRangePID_h
#define pdUtilsRangePID_h

#include "pdDataClasses.hxx"
#include <utility>

namespace pdAnaUtils {

  Float_t ComputeRangeMomentum(double trkrange, int pdg);
  Float_t ComputeCSDARange(double beammom, int pdg);
  Float_t ComputeKineticEnergy(const AnaParticlePD &part);

  std::pair<double, int> Chi2PID(const AnaParticlePD& part, const int pdg);
  /// Chi2PID on the hits of one wire plane (Chi2PID is plane 2, the collection plane); the dE/dx resolution model is
  /// the collection-plane one for every plane.
  std::pair<double, int> Chi2PIDPlane(const AnaParticlePD& part, const int pdg, const int plane);
  std::pair<double, int> Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR);

  /// Chi2PID first / second (mean χ² per contributing hit); -999 if null or Chi2PID invalid. Same as neutral-kaon truth tree Reco*Chi2Ndf.
  Float_t Chi2PIDChi2PerHit(const AnaParticlePD* part, int pdg);
  /// As Chi2PIDChi2PerHit on the induction plane with more stored hits (filled only from PDSPAnalyzer ntuples with
  /// pdIO.PDSPAnalyzerFillAllPlanes = 1); -999 if there are none.
  Float_t Chi2PIDChi2PerHitInduction(const AnaParticlePD* part, int pdg);

  /// Pion (211): same per-hit χ² as Chi2PID (dedx_range_pi template + bin error + dedx_res), but only collection hits with
  /// 0 < RR < maxResidualRangeCm and the same skip / optional dE/dx window as ComputePionBraggWindowChi2PiEq61.
  bool ComputePionBraggWindowChi2PionRangeTemplate(AnaParticlePD* part, double maxResidualRangeCm, int minHits,
                                                   int skipHitsFirst, int skipHitsLast, double dedxMinMeVcm,
                                                   double dedxMaxMeVcm, double& meanChi2PerHit, int& nHitsUsed);

}

#endif
