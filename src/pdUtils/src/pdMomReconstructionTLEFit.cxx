#include "pdMomReconstructionTLEFit.hxx"
#include "pdUtilsDEdx.hxx"

#include <cmath>
#include <limits>

namespace pdMomReconstruction {

bool EstimatePionMomentumTLEFit(AnaParticlePD* track, const PionTLEFitConfig& cfg, PionTLEFitResult& out,
                                std::vector<double>* pAxisGeVOut, std::vector<double>* logLOut) {
  out.bestMomentumGeV = -999.0;
  out.bestLogLikelihood = -999.0;
  out.valid = false;
  if (pAxisGeVOut) pAxisGeVOut->clear();
  if (logLOut) logLOut->clear();
  if (!track) return false;

  std::vector<double> pAxisGeV;
  std::vector<double> logL;
  if (!pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
          track, cfg.scanLmaxCm, cfg.scanStepCm, cfg.minInteriorHits, cfg.skipHitsFirst, cfg.skipHitsLast,
          cfg.dedxMinMeVcm, cfg.dedxMaxMeVcm, cfg.dedxPdfPathCm, pAxisGeV, logL, cfg.scanStepFineCm,
          cfg.lowPMomentumRefineGeV) ||
      logL.size() != pAxisGeV.size())
    return false;

  double bestLogL = -std::numeric_limits<double>::infinity();
  double bestP = -999.0;
  for (size_t i = 0; i < pAxisGeV.size(); ++i) {
    if (!std::isfinite(logL[i])) continue;
    if (logL[i] > bestLogL) {
      bestLogL = logL[i];
      bestP = pAxisGeV[i];
    }
  }
  if (!std::isfinite(bestP) || bestP <= 0.0) return false;

  out.bestMomentumGeV = bestP;
  out.bestLogLikelihood = std::isfinite(bestLogL) ? bestLogL : -999.0;
  out.valid = true;
  if (pAxisGeVOut) *pAxisGeVOut = pAxisGeV;
  if (logLOut) *logLOut = logL;
  return true;
}

}  // namespace pdMomReconstruction
