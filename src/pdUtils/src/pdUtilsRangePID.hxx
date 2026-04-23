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

}

#endif
