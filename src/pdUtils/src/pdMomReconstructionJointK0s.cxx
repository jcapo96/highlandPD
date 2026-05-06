#include "pdMomReconstructionJointK0s.hxx"
#include "pdMomReconstructionTLEFit.hxx"
#include "pdUtilsDEdx.hxx"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pdMomReconstruction {
namespace {

double ArgmaxMomentumFromCurve(const std::vector<double>& pAxis, const std::vector<double>& logL) {
  if (pAxis.empty() || pAxis.size() != logL.size()) return std::numeric_limits<double>::quiet_NaN();
  size_t im = 0;
  for (size_t i = 1; i < logL.size(); ++i) {
    if (logL[i] > logL[im]) im = i;
  }
  return pAxis[im];
}

}  // namespace

bool BuildNeutralKaonJointDiagnosticsCurvesForDaughter(AnaParticlePD* reco, const PionTLEFitConfig& tle,
                                                       const MCSLikelihoodConfig& mcsLikelihood,
                                                       double diagPMinGeV, std::vector<double>& pAxisDiagOut,
                                                       std::vector<double>& rawLogLTleOut,
                                                       std::vector<double>& rawLogLMcsOut, double& pBestTleOut,
                                                       double& pBestMcsOut) {
  pAxisDiagOut.clear();
  rawLogLTleOut.clear();
  rawLogLMcsOut.clear();
  pBestTleOut = -1.;
  pBestMcsOut = -1.;
  if (!reco) return false;

  std::vector<double> pAxis;
  std::vector<double> logL;
  if (!pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
          reco, tle.scanLmaxCm, tle.scanStepCm, tle.minInteriorHits, tle.skipHitsFirst, tle.skipHitsLast,
          tle.dedxMinMeVcm, tle.dedxMaxMeVcm, tle.dedxPdfPathCm, pAxis, logL, tle.scanStepFineCm,
          tle.lowPMomentumRefineGeV))
    return false;
  if (pAxis.empty() || pAxis.size() != logL.size()) return false;

  pAxisDiagOut.reserve(pAxis.size() + 64u);
  const double diagPMin = (std::isfinite(diagPMinGeV) && diagPMinGeV > 0.0) ? diagPMinGeV : 0.01;
  double pStepDiag = 0.05;
  if (pAxis.size() >= 2u) {
    std::vector<double> dP;
    dP.reserve(pAxis.size() - 1u);
    for (size_t ip = 1; ip < pAxis.size(); ++ip) {
      const double dp = pAxis[ip] - pAxis[ip - 1u];
      if (std::isfinite(dp) && dp > 1e-9) dP.push_back(dp);
    }
    if (!dP.empty()) {
      const size_t mid = dP.size() / 2u;
      std::nth_element(dP.begin(), dP.begin() + static_cast<std::ptrdiff_t>(mid), dP.end());
      pStepDiag = std::max(1e-3, dP[mid]);
    }
  }
  const double pStartDiag = (std::isfinite(diagPMin) && diagPMin > 0.0) ? diagPMin : 0.01;
  if (pStartDiag < pAxis.front() - 1e-9) {
    for (double p = pStartDiag; p < pAxis.front() - 0.5 * pStepDiag; p += pStepDiag) pAxisDiagOut.push_back(p);
  }
  pAxisDiagOut.insert(pAxisDiagOut.end(), pAxis.begin(), pAxis.end());
  if (pAxisDiagOut.empty()) return false;

  rawLogLTleOut.reserve(pAxisDiagOut.size());
  for (double p : pAxisDiagOut)
    rawLogLTleOut.push_back(pdJointK0sPionMomentum::InterpolateLogLikelihoodClamped(pAxis, logL, p));
  pBestTleOut = ArgmaxMomentumFromCurve(pAxisDiagOut, rawLogLTleOut);
  if (!std::isfinite(pBestTleOut) || pBestTleOut <= 0.) pBestTleOut = -1.;

  std::vector<double> logLMcs;
  const bool hasMcs = BuildPionMCSLogLikelihoodVsMomentumCurve(*reco, pAxisDiagOut, mcsLikelihood, logLMcs) &&
                      logLMcs.size() == pAxisDiagOut.size();
  if (hasMcs) {
    rawLogLMcsOut = std::move(logLMcs);
    pBestMcsOut = ArgmaxMomentumFromCurve(pAxisDiagOut, rawLogLMcsOut);
    if (!std::isfinite(pBestMcsOut) || pBestMcsOut <= 0.) pBestMcsOut = -1.;
  }
  return true;
}

bool FitJointK0sTwoPionMomentaOnGrid(AnaParticlePD* daughter1, AnaParticlePD* daughter2, const TVector3& dir1,
                                    const TVector3& dir2, const JointK0sTwoPionGridFitConfig& cfg,
                                    JointK0sPionMomentumGridResult& out, JointK0sSigmaMEventDiagnostics* sigmaDiagOut,
                                    JointK0sTwoPionFitAuxiliary* auxOut) {
  out = JointK0sPionMomentumGridResult{};
  if (sigmaDiagOut) *sigmaDiagOut = JointK0sSigmaMEventDiagnostics{};
  if (auxOut) *auxOut = JointK0sTwoPionFitAuxiliary{};
  if (!daughter1 || !daughter2) return false;

  PionTLEFitResult r1;
  PionTLEFitResult r2;
  std::vector<double> p1v;
  std::vector<double> logL1;
  std::vector<double> p2v;
  std::vector<double> logL2;
  if (!EstimatePionMomentumTLEFit(daughter1, cfg.tle, r1, &p1v, &logL1)) return false;
  if (!EstimatePionMomentumTLEFit(daughter2, cfg.tle, r2, &p2v, &logL2)) return false;
  if (p1v.size() != logL1.size() || p2v.size() != logL2.size() || p1v.empty() || p2v.empty()) return false;

  if (auxOut) {
    auxOut->p1AxisGeV = p1v;
    auxOut->logL1Tle = logL1;
    auxOut->p2AxisGeV = p2v;
    auxOut->logL2Tle = logL2;
    auxOut->p1TleArgmaxGeV = ArgmaxMomentumFromCurve(p1v, logL1);
    auxOut->p2TleArgmaxGeV = ArgmaxMomentumFromCurve(p2v, logL2);
    if (!std::isfinite(auxOut->p1TleArgmaxGeV)) auxOut->p1TleArgmaxGeV = -999.;
    if (!std::isfinite(auxOut->p2TleArgmaxGeV)) auxOut->p2TleArgmaxGeV = -999.;
  }

  std::vector<double> logL1Joint = logL1;
  std::vector<double> logL2Joint = logL2;
  if (cfg.useMCS && std::isfinite(cfg.mcsWeight) && cfg.mcsWeight > 0.) {
    std::vector<double> logLMcs1;
    std::vector<double> logLMcs2;
    const bool okMcs1 = BuildPionMCSLogLikelihoodVsMomentumCurve(*daughter1, p1v, cfg.mcs.likelihood, logLMcs1);
    const bool okMcs2 = BuildPionMCSLogLikelihoodVsMomentumCurve(*daughter2, p2v, cfg.mcs.likelihood, logLMcs2);
    if (okMcs1 && okMcs2 && logLMcs1.size() == logL1Joint.size() && logLMcs2.size() == logL2Joint.size()) {
      if (auxOut) {
        auxOut->p1McsArgmaxGeV = ArgmaxMomentumFromCurve(p1v, logLMcs1);
        auxOut->p2McsArgmaxGeV = ArgmaxMomentumFromCurve(p2v, logLMcs2);
        if (!std::isfinite(auxOut->p1McsArgmaxGeV)) auxOut->p1McsArgmaxGeV = -999.;
        if (!std::isfinite(auxOut->p2McsArgmaxGeV)) auxOut->p2McsArgmaxGeV = -999.;
      }
      for (size_t i = 0; i < logL1Joint.size(); ++i) logL1Joint[i] += cfg.mcsWeight * logLMcs1[i];
      for (size_t i = 0; i < logL2Joint.size(); ++i) logL2Joint[i] += cfg.mcsWeight * logLMcs2[i];
    }
  }

  const double sigmaLo = std::min(cfg.sigmaMassMinGeV, cfg.sigmaMassMaxGeV);
  const double sigmaHi = std::max(cfg.sigmaMassMinGeV, cfg.sigmaMassMaxGeV);
  double sigma_m_for_grid = cfg.sigmaMassGeV;
  JointK0sSigmaMEventDiagnostics diag{};
  if (cfg.useEventSigmaM &&
      pdJointK0sPionMomentum::ComputeSigmaMEventGeV(p1v, logL1, p2v, logL2, dir1, dir2, cfg.sigmaMassGeV, diag,
                                                    sigmaLo, sigmaHi) &&
      std::isfinite(diag.sigma_m_event_gev) && diag.sigma_m_event_gev > 0.) {
    sigma_m_for_grid = diag.sigma_m_event_gev;
    if (sigmaDiagOut) *sigmaDiagOut = diag;
  } else if (sigmaDiagOut) {
    *sigmaDiagOut = JointK0sSigmaMEventDiagnostics{};
  }
  if (std::isfinite(sigma_m_for_grid) && sigma_m_for_grid > 0.) {
    sigma_m_for_grid = std::max(sigma_m_for_grid, sigmaLo);
    sigma_m_for_grid = std::min(sigma_m_for_grid, sigmaHi);
  }

  out = pdJointK0sPionMomentum::FitJointMomentaOnGrid(p1v, logL1Joint, p2v, logL2Joint, dir1, dir2, cfg.pMinGeV,
                                                     cfg.pMaxGeV, cfg.pStepGeV, cfg.mK0sMassGeV, sigma_m_for_grid,
                                                     cfg.massPenaltyScale, cfg.refineFactor);
  return true;
}

}  // namespace pdMomReconstruction
