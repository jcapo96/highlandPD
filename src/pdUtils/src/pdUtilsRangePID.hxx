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

}

#endif
