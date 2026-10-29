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
  /// Multiplies the measured dE/dx before the fit (calibration to the likelihood model; 1 = off).
  double dedxScale = 1.0;
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

/// Recommended pion momentum for K0S daughters (pionMomentum/interacting study): TLE with the CSDA momentum of
/// the reco length as a lower bound (an interacting or stopping pion has at least that momentum), CSDA alone if
/// the TLE fit fails. "measured" flags momenta below the minimum-ionising plateau, where dE/dx still constrains
/// the momentum; above it the value is only a lower bound (TLE saturates near 0.5 GeV/c).
struct PionMomentumRecoConfig {
  PionTLEFitConfig tle;
  double measuredMaxGeV = 0.45;
};

struct PionMomentumRecoResult {
  double momentumGeV = -999.0;
  double tleMomentumGeV = -999.0;
  double rangeMomentumGeV = -999.0;
  bool tleValid = false;
  bool measured = false;
  bool valid = false;
};

bool EstimatePionMomentumRecommended(AnaParticlePD* track, const PionMomentumRecoConfig& cfg,
                                     PionMomentumRecoResult& out);

}  // namespace pdMomReconstruction

#endif
