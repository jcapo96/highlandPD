#ifndef pdMomReconstructionTLEFit_h
#define pdMomReconstructionTLEFit_h

#include "pdDataClasses.hxx"
#include <vector>

namespace pdMomReconstruction {

/// Configuration for pion TLE (free-range dE/dx) momentum scan — mirrors pionMomentumAnalysis free-range parameters.
struct PionTLEFitConfig {
  int minInteriorHits = 15;
  int skipHitsFirst = 3;
  int skipHitsLast = 3;
  double dedxMinMeVcm = 0.5;
  double dedxMaxMeVcm = 5.0;
  double scanLmaxCm = 450.0;
  double scanStepCm = 1.0;
  double scanStepFineCm = 0.25;
  double lowPMomentumRefineGeV = 0.4;
  double dedxPdfPathCm = 0.65;
};

struct PionTLEFitResult {
  double bestMomentumGeV = -999.0;
  double bestLogLikelihood = -999.0;
  bool valid = false;
};

/// Best momentum from the TLE / free-range dE/dx likelihood scan (collection plane hits).
/// If \a pAxisGeVOut and \a logLOut are non-null, fills them with the scan curve (for diagnostics).
bool EstimatePionMomentumTLEFit(AnaParticlePD* track, const PionTLEFitConfig& cfg, PionTLEFitResult& out,
                                 std::vector<double>* pAxisGeVOut = nullptr, std::vector<double>* logLOut = nullptr);

}  // namespace pdMomReconstruction

#endif
