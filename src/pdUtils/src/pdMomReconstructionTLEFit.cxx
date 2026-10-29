#include "pdMomReconstructionTLEFit.hxx"
#include "pdMomReconstruction.hxx"
#include "pdUtilsDEdx.hxx"

#include <algorithm>
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
          cfg.lowPMomentumRefineGeV, cfg.dedxScale) ||
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

bool EstimatePionMomentumRecommended(AnaParticlePD* track, const PionMomentumRecoConfig& cfg,
                                     PionMomentumRecoResult& out) {
  out = PionMomentumRecoResult();
  if (!track) return false;
  const Float_t pRange = EstimatePionMomentumFromCSDA(track);
  out.rangeMomentumGeV = (pRange > 0.f) ? static_cast<double>(pRange) : -999.0;
  PionTLEFitResult tle;
  out.tleValid = EstimatePionMomentumTLEFit(track, cfg.tle, tle) && tle.valid;
  if (out.tleValid) out.tleMomentumGeV = tle.bestMomentumGeV;
  if (out.tleValid)
    out.momentumGeV = std::max(out.tleMomentumGeV, out.rangeMomentumGeV);
  else
    out.momentumGeV = out.rangeMomentumGeV;
  out.valid = out.momentumGeV > 0.0;
  out.measured = out.valid && out.momentumGeV < cfg.measuredMaxGeV;
  return out.valid;
}

}  // namespace pdMomReconstruction
