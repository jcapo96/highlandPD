#include "pdUtilsNeutralKaonSignalDiagnostics.hxx"

#include "CategoryManager.hxx"
#include "Parameters.hxx"
#include "TGraph.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TTree.h"
#include "neutralKaonAnalysisUtils.hxx"
#include "pdJointK0sPionMomentum.hxx"
#include "pdMomReconstructionMCS.hxx"
#include "pdMomReconstructionJointK0s.hxx"
#include "pdMomReconstructionFromParams.hxx"
#include "pdUtilsDEdx.hxx"
#include "pdUtilsLineFit.hxx"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <tuple>
#include <unordered_set>
#include <limits>
#include <TH2F.h>
#include <TGaxis.h>
#include <TVector3.h>

namespace {
  Int_t sSignalDiagRun = -999999999;
  Int_t sSignalDiagSub = -999999999;
  Int_t sSignalDiagEvt = -999999999;
  std::unordered_set<Int_t> sSignalDiagTrueK0Ids;
  std::unordered_set<std::string> sSigDedxAcceptedDauKeys;
  std::unordered_map<std::string, int> sSigDedxPairIndexByDauKey;
  std::unordered_map<int, int> sK0SigDedxSerialByCode;

  std::vector<TMultiGraph*> gK0SignalDedxMultiGraphs;
  std::vector<TH1F*> gK0SignalDedxBiasHistograms;
  std::vector<TMultiGraph*> gK0JointMomentumLogLMultiGraphs;
  std::vector<TMultiGraph*> gK0JointMomentumMCSLogLMultiGraphs;
  std::vector<TGraph*> gK0JointMomentumMCSAngleGraphs;
  std::vector<TGraph*> gK0JointMomentumMCSDedxVsRRGraphs;
  struct McsAngleStats {
    int entries = 0;
    double meanObs = std::numeric_limits<double>::quiet_NaN();
    double rmsObs = std::numeric_limits<double>::quiet_NaN();
    double rmsTrue = std::numeric_limits<double>::quiet_NaN();
    bool detectorSigmaEnabled = false;
    bool detectorSigmaApplied = false;
    double rmsDetUsed = std::numeric_limits<double>::quiet_NaN();
    double rmsMcsFromObsSubDet = std::numeric_limits<double>::quiet_NaN();
  };
  std::unordered_map<std::string, McsAngleStats> gK0JointMomentumMCSAngleStatsByGraphName;
  std::vector<TH2F*> gK0JointObjectiveTH2Sum;
  std::vector<TH2F*> gK0JointObjectiveTH2Penalty;
  std::vector<TH2F*> gK0JointObjectiveTH2TrackLogLSum;
  std::vector<double> gK0JointObjectiveOpeningAngleDeg;
  std::vector<double> gK0JointObjectiveSigmaMGeV;
  std::vector<TGraph*> gK0JointObjectiveBestMarker;
  /// One-point marker TGraphs per joint 2D plot: TLE, true, joint (flat, 3 per TH2; nullptr if unavailable).
  std::vector<TGraph*> gK0JointObjective2DMarkerGraphs;

  void DeleteMultiGraphAndGraphs(TMultiGraph* mg) {
    if (!mg) return;
    TList* gl = mg->GetListOfGraphs();
    if (gl) gl->Delete();
    delete mg;
  }

  /// Y limits for log-likelihood pads: cap vertical span (deep tails from low-p extrapolation) and add margin.
  bool LogLYAxisRangeFromValues(const std::vector<double>& ys, double maxSpan, double marginFrac, double& yLow,
                                double& yHigh) {
    if (ys.empty()) return false;
    const double yHi = *std::max_element(ys.begin(), ys.end());
    const double yLoRaw = *std::min_element(ys.begin(), ys.end());
    if (!std::isfinite(yHi) || !std::isfinite(yLoRaw)) return false;
    double yLo = yLoRaw;
    if (!(yHi > yLo)) {
      yLo = yHi - 1.0;
    } else {
      yLo = std::max(yLo, yHi - maxSpan);
    }
    const double pad = marginFrac * std::max(1e-9, yHi - yLo);
    yLow = yLo - pad;
    yHigh = yHi + pad;
    return true;
  }

  void ClipLogLMultigraphYAxisAfterDraw(TMultiGraph* mg, double maxSpan = 450., double marginFrac = 0.08) {
    if (!mg || !gPad) return;
    TList* gl = mg->GetListOfGraphs();
    if (!gl) return;
    // Use NLL curve points only (skip 2-point vertical markers). Marker Y ranges are built per-pad from the
    // same curve samples so they stay inside this window.
    std::vector<double> ys;
    TIter iter(gl);
    TObject* obj = nullptr;
    while ((obj = iter.Next()) != nullptr) {
      auto* gr = dynamic_cast<TGraph*>(obj);
      if (!gr || gr->GetN() <= 2) continue;
      for (int i = 0; i < gr->GetN(); ++i) {
        const double py = gr->GetPointY(i);
        if (std::isfinite(py)) ys.push_back(py);
      }
    }
    double yLo = 0.;
    double yHi = 1.;
    if (!LogLYAxisRangeFromValues(ys, maxSpan, marginFrac, yLo, yHi)) return;
    gPad->Update();
    TH1* h = mg->GetHistogram();
    if (h) h->GetYaxis()->SetRangeUser(yLo, yHi);
    gPad->Modified();
    gPad->Update();
  }

  void AddLogLPadLegend(TMultiGraph* mg, bool isMcsPad) {
    if (!mg || !gPad) return;
    TList* gl = mg->GetListOfGraphs();
    if (!gl) return;
    TGraph* gCurve = nullptr;
    TGraph* gBest = nullptr;
    TGraph* gJoint = nullptr;
    TGraph* gTrue = nullptr;
    TGraph* gRms = nullptr;
    TIter iter(gl);
    TObject* obj = nullptr;
    while ((obj = iter.Next()) != nullptr) {
      auto* gr = dynamic_cast<TGraph*>(obj);
      if (!gr || !gr->GetName()) continue;
      const std::string name = gr->GetName();
      if (name.find("g_k0_joint_nll_") != std::string::npos) {
        gCurve = gr;
      } else if (name.find("g_k0_joint_ptle_") != std::string::npos || name.find("g_k0_joint_pmcs_") != std::string::npos) {
        gBest = gr;
      } else if (name.find("g_k0_joint_pjoint_") != std::string::npos) {
        gJoint = gr;
      } else if (name.find("g_k0_joint_truep_") != std::string::npos) {
        gTrue = gr;
      } else if (name.find("g_k0_joint_prms_") != std::string::npos) {
        gRms = gr;
      }
    }
    TLegend* leg = new TLegend(0.50, 0.60, 0.88, 0.89);
    if (!leg) return;
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (gCurve) leg->AddEntry(gCurve, isMcsPad ? "logL_{MCS}(p)" : "logL_{TLE}(p)", "l");
    if (gBest) leg->AddEntry(gBest, isMcsPad ? "p_{best,MCS}" : "p_{best,TLE}", "l");
    if (gJoint) leg->AddEntry(gJoint, "p_{joint}", "l");
    if (gTrue) leg->AddEntry(gTrue, "p_{true}", "l");
    if (gRms) leg->AddEntry(gRms, "p_{RMS(HL)}", "l");
    leg->Draw();
  }

  bool IsSignalLikeCode(int signalCode) {
    return (signalCode == 1 || signalCode == 5 || signalCode == 6);
  }

  /// GeV/c from theta0 ≈ RMS(Δθ), mean chord segment (cm), same Highland/Lynch-Dahl factors as MCSLikelihood.
  double PionMomentumGeVFromHighlandRmsAndMeanSeg(double rmsThetaRad, double meanSegLenCm, double x0Cm) {
    constexpr double kPionMassMeV = 139.57;
    constexpr double kHighlandMeV = 13.6;
    if (!std::isfinite(rmsThetaRad) || rmsThetaRad <= 0.) return -1.;
    if (!std::isfinite(meanSegLenCm) || meanSegLenCm <= 0.) return -1.;
    if (!std::isfinite(x0Cm) || x0Cm <= 0.) return -1.;
    const double xOverX0 = meanSegLenCm / x0Cm;
    if (!std::isfinite(xOverX0) || xOverX0 <= 0.) return -1.;
    double corr = 1.0 + 0.038 * std::log(std::max(xOverX0, 1e-12));
    if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;
    const double K = kHighlandMeV * std::sqrt(xOverX0) * corr / rmsThetaRad;
    if (!std::isfinite(K) || K <= 0.) return -1.;
    const double K2 = K * K;
    const double m2 = kPionMassMeV * kPionMassMeV;
    const double disc = K2 * K2 + 4.0 * K2 * m2;
    if (!std::isfinite(disc) || disc < 0.) return -1.;
    const double u = 0.5 * (K2 + std::sqrt(disc));
    if (!std::isfinite(u) || u <= 0.) return -1.;
    const double pMeV = std::sqrt(u);
    const double pGeV = pMeV / 1000.0;
    return (std::isfinite(pGeV) && pGeV > 0.) ? pGeV : -1.;
  }

  bool BuildMcsSegments(const AnaParticlePD& track, const pdMomReconstruction::MCSLikelihoodConfig& cfg,
                        std::vector<double>& deltaTheta, std::vector<double>& xOverX0, std::vector<double>& rrMidCm) {
    return pdMomReconstruction::BuildPionMcsScatteringSamples(track, cfg, deltaTheta, xOverX0, &rrMidCm);
  }

  /// All collection-plane (view 2) dE/dx vs residual range hits — same convention as dEdx likelihood utilities.
  bool CollectCollectionPlaneDedxVsRR(const AnaParticlePD& track, std::vector<double>& rrOut,
                                      std::vector<double>& dedxOut) {
    rrOut.clear();
    dedxOut.clear();
    if (track.Hits[2].empty()) return false;
    for (const auto& h : track.Hits[2]) {
      if (!std::isfinite(static_cast<double>(h.ResidualRange)) || h.ResidualRange < -900.f) continue;
      const double rr = static_cast<double>(h.ResidualRange);
      double dx = static_cast<double>(h.dEdx);
      if (!std::isfinite(dx) || dx <= 0.) dx = static_cast<double>(h.dEdx_NoSCE);
      if (!std::isfinite(rr) || !std::isfinite(dx) || dx <= 0.) continue;
      rrOut.push_back(rr);
      dedxOut.push_back(dx);
    }
    return !rrOut.empty();
  }

  int GetSignalCodeFromCategory(size_t candidateIndex) {
    if (!anaUtils::_categ || !anaUtils::_categ->HasCategory("signal")) return 2;
    return anaUtils::_categ->GetCategory("signal").GetObjectCode(1, static_cast<Int_t>(candidateIndex));
  }

  std::string MakeDiagDauKey(const AnaEventB& event, size_t candidateIndex, int dauIndex) {
    Int_t runId = -1;
    Int_t subRunId = -1;
    Int_t evtId = -1;
    if (event.EventInfo) {
      runId = event.EventInfo->Run;
      subRunId = event.EventInfo->SubRun;
      evtId = event.EventInfo->Event;
    }
    return Form("%d_%d_%d_%zu_%d", static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId),
                candidateIndex, dauIndex);
  }

  void DiagAnnihilationDaughterFitDirs(AnaAnnihilationVertexPD* vertex, double trackFitLength,
                                      double trackFitDistanceFromStart, TVector3& dirFit1, TVector3& dirFit2) {
    dirFit1.SetXYZ(0., 0., 0.);
    dirFit2.SetXYZ(0., 0., 0.);
    if (!vertex || vertex->Particles.size() < 2) return;
    AnaParticlePD* daughter1 = vertex->Particles[0];
    AnaParticlePD* daughter2 = vertex->Particles[1];
    if (!daughter1 || !daughter2) return;
    dirFit1.SetXYZ(daughter1->DirectionStart[0], daughter1->DirectionStart[1], daughter1->DirectionStart[2]);
    dirFit2.SetXYZ(daughter2->DirectionStart[0], daughter2->DirectionStart[1], daughter2->DirectionStart[2]);
    std::vector<double> fit1;
    std::vector<double> fit2;
    pdAnaUtils::ExtrapolateTrack(daughter1, fit1, trackFitLength, true, trackFitDistanceFromStart);
    pdAnaUtils::ExtrapolateTrack(daughter2, fit2, trackFitLength, true, trackFitDistanceFromStart);
    const bool fit1Valid = (fit1.size() >= 6 && fit1[3] > -900.0 && fit1[4] > -900.0 && fit1[5] > -900.0 &&
                            std::isfinite(fit1[3]) && std::isfinite(fit1[4]) && std::isfinite(fit1[5]));
    const bool fit2Valid = (fit2.size() >= 6 && fit2[3] > -900.0 && fit2[4] > -900.0 && fit2[5] > -900.0 &&
                            std::isfinite(fit2[3]) && std::isfinite(fit2[4]) && std::isfinite(fit2[5]));
    if (fit1Valid) dirFit1.SetXYZ(fit1[3], fit1[4], fit1[5]);
    if (fit2Valid) dirFit2.SetXYZ(fit2[3], fit2[4], fit2[5]);
  }
}

namespace neutralKaonTreeDiagnostics {

bool RegisterSignalTrueK0Id(Int_t trueK0Id) {
  return sSignalDiagTrueK0Ids.insert(trueK0Id).second;
}

void MaybeAccumulateSignalPionDedxMultiGraphs(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                              size_t candidateIndex) {
  if (!candidate) return;
  const int signalCode = GetSignalCodeFromCategory(candidateIndex);
  const bool ensureSignalOnly = ND::params().HasParameter("neutralKaonAnalysis.EnsureMomentumSignalOnly") &&
                                ND::params().GetParameterI("neutralKaonAnalysis.EnsureMomentumSignalOnly") == 1;
  if (ensureSignalOnly && !IsSignalLikeCode(signalCode)) return;

  AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
  if (!vertex || vertex->Particles.size() < 2) return;

  const int pairIdx = sK0SigDedxSerialByCode[signalCode];
  Int_t runId = -1;
  Int_t subRunId = -1;
  Int_t evtId = -1;
  if (event.EventInfo) {
    runId = event.EventInfo->Run;
    subRunId = event.EventInfo->SubRun;
    evtId = event.EventInfo->Event;
  }

  bool anyGraph = false;
  for (int idau = 0; idau < 2; ++idau) {
    AnaParticlePD* reco = static_cast<AnaParticlePD*>(vertex->Particles[idau]);
    if (!reco) continue;

    const double Lmax = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                            ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanLmaxCm")
                            : 450.;
    const double step = ND::params().HasParameter("neutralKaonAnalysis.FreeRangeScanStepCm")
                            ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeScanStepCm")
                            : 1.;
    const int minInterior =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
            ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxMinInteriorHits")
            : 15;
    const int skipFirst =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
            ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsFirst")
            : 3;
    const int skipLast =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
            ? ND::params().GetParameterI("neutralKaonAnalysis.FreeRangeDedxSkipHitsLast")
            : 3;
    const double dedxMin =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
            ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMinMeVcm")
            : 0.5;
    const double dedxMax =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
            ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxMaxMeVcm")
            : 5.0;
    const double pdfPath =
        ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
            ? ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxPdfPathCm")
            : 0.65;

    char dedxXAxis[200];
    std::snprintf(dedxXAxis, sizeof(dedxXAxis),
                  "(run %d, subrun %d, evt %d) Residual range [cm]", static_cast<int>(runId),
                  static_cast<int>(subRunId), static_cast<int>(evtId));
    TMultiGraph* mg = pdAnaUtils::MakePionFreeRangeDedxVsRRMultiGraph(
        reco, Lmax, step, skipFirst, skipLast, dedxMin, dedxMax, minInterior, pdfPath, dedxXAxis);
    if (!mg) continue;

    char biasHistTitle[240];
    std::snprintf(biasHistTitle, sizeof(biasHistTitle),
                  "(run %d, subrun %d, evt %d) #Delta(dE/dx)=measured-expected(PDF mode);#Delta(dE/dx) [MeV/cm];Entries",
                  static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId));
    TH1F* hBias = pdAnaUtils::MakePionFreeRangeDedxBiasHistogram(reco, Lmax, step, skipFirst, skipLast, dedxMin,
                                                                 dedxMax, minInterior, pdfPath, biasHistTitle);
    if (!hBias) {
      DeleteMultiGraphAndGraphs(mg);
      continue;
    }

    mg->SetName(Form("mg_sigdedx_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    hBias->SetName(Form("h_sigdedx_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    gK0SignalDedxMultiGraphs.push_back(mg);
    gK0SignalDedxBiasHistograms.push_back(hBias);
    const std::string dauKey = MakeDiagDauKey(event, candidateIndex, idau);
    sSigDedxAcceptedDauKeys.insert(dauKey);
    sSigDedxPairIndexByDauKey[dauKey] = pairIdx;
    anyGraph = true;
  }
  if (anyGraph) ++sK0SigDedxSerialByCode[signalCode];
}

void MaybeAccumulateJointMomentumLogLikelihoodGraphs(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                                     size_t candidateIndex) {
  if (!candidate) return;
  AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
  if (!vertex || vertex->Particles.size() < 2) return;
  const int signalCode = GetSignalCodeFromCategory(candidateIndex);

  const bool ensureSignalOnly = ND::params().HasParameter("neutralKaonAnalysis.EnsureMomentumSignalOnly") &&
                                ND::params().GetParameterI("neutralKaonAnalysis.EnsureMomentumSignalOnly") == 1;
  if (ensureSignalOnly) {
    if (!IsSignalLikeCode(signalCode)) return;
  }

  pdMomReconstruction::JointK0sTwoPionGridFitConfig jcfg;
  pdMomReconstruction::FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams(jcfg);
  const pdMomReconstruction::MCSLikelihoodConfig& mcsCfg = jcfg.mcs.likelihood;
  const double diagPMin =
      ND::params().HasParameter("neutralKaonAnalysis.JointK0sDiagMomentumPMinGeV")
          ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sDiagMomentumPMinGeV")
          : 0.01;

  for (int idau = 0; idau < 2; ++idau) {
    const std::string dauKey = MakeDiagDauKey(event, candidateIndex, idau);
    if (sSigDedxAcceptedDauKeys.find(dauKey) == sSigDedxAcceptedDauKeys.end()) {
      continue;
    }
    auto itPairIdx = sSigDedxPairIndexByDauKey.find(dauKey);
    if (itPairIdx == sSigDedxPairIndexByDauKey.end()) continue;
    const int pairIdx = itPairIdx->second;
    AnaParticlePD* reco = static_cast<AnaParticlePD*>(vertex->Particles[idau]);
    if (!reco) continue;

    std::vector<double> pAxisDiag;
    std::vector<double> rawLogLTle;
    std::vector<double> rawLogLMcs;
    double pBestTle = -1.;
    double pBestMcs = -1.;
    if (!pdMomReconstruction::BuildNeutralKaonJointDiagnosticsCurvesForDaughter(
            reco, jcfg.tle, mcsCfg, diagPMin, pAxisDiag, rawLogLTle, rawLogLMcs, pBestTle, pBestMcs))
      continue;
    if (pAxisDiag.empty() || rawLogLTle.size() != pAxisDiag.size()) continue;
    if (!std::isfinite(pBestTle) || pBestTle <= 0.) continue;

    const double pJoint =
        (idau == 0) ? static_cast<double>(vertex->Daughter1MomentumJoint) : static_cast<double>(vertex->Daughter2MomentumJoint);

    std::vector<double> yTleOnly;
    yTleOnly.reserve(rawLogLTle.size());
    for (double y : rawLogLTle) {
      if (std::isfinite(y)) yTleOnly.push_back(y);
    }
    constexpr double kMaxLogLSpan = 450.;
    double yMinTle = 0.0;
    double yMaxTle = 1.0;
    if (!LogLYAxisRangeFromValues(yTleOnly, kMaxLogLSpan, 0.08, yMinTle, yMaxTle)) {
      yMinTle = 0.0;
      yMaxTle = 1.0;
    }
    const double yTleLine[2] = {yMinTle, yMaxTle};

    TGraph* gNllTle = new TGraph(static_cast<int>(pAxisDiag.size()), pAxisDiag.data(), rawLogLTle.data());
    if (!gNllTle) continue;
    const double xBestTle[2] = {pBestTle, pBestTle};
    TGraph* gBestTleLine = new TGraph(2, xBestTle, yTleLine);
    if (!gBestTleLine) {
      delete gNllTle;
      continue;
    }
    TMultiGraph* mg = new TMultiGraph();
    if (!mg) {
      delete gNllTle;
      delete gBestTleLine;
      continue;
    }
    TMultiGraph* mgMcs = nullptr;

    mg->SetName(Form("mg_siglogl_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    mg->SetTitle(Form("TLEFit pion log-likelihood diagnostics (pair %d, dau %d, signalCode %d);Momentum [GeV/c];logL_{TLE}",
                      pairIdx, idau + 1, signalCode));
    gNllTle->SetName(Form("g_k0_joint_nll_tle_pair%d_dau%d", pairIdx, idau + 1));
    gNllTle->SetLineColor(kBlack);
    gNllTle->SetLineWidth(2);
    gBestTleLine->SetName(Form("g_k0_joint_ptle_pair%d_dau%d", pairIdx, idau + 1));
    gBestTleLine->SetLineColor(kBlue + 1);
    gBestTleLine->SetLineStyle(2);
    gBestTleLine->SetLineWidth(2);
    mg->Add(gNllTle, "L");
    mg->Add(gBestTleLine, "L");

    if (std::isfinite(pJoint) && pJoint > 0.) {
      const double xJoint[2] = {pJoint, pJoint};
      const double yJointLineTle[2] = {yMinTle, yMaxTle};
      TGraph* gJointLine = new TGraph(2, xJoint, yJointLineTle);
      if (gJointLine) {
        gJointLine->SetName(Form("g_k0_joint_pjoint_pair%d_dau%d", pairIdx, idau + 1));
        gJointLine->SetLineColor(kMagenta + 1);
        gJointLine->SetLineStyle(3);
        gJointLine->SetLineWidth(2);
        mg->Add(gJointLine, "L");
      }
    }

    const AnaTrueParticlePD* trueReco =
        reco->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco->TrueObject) : nullptr;
    if (trueReco && std::isfinite(static_cast<double>(trueReco->Momentum)) && trueReco->Momentum > 0.f) {
      const double xTrue[2] = {static_cast<double>(trueReco->Momentum), static_cast<double>(trueReco->Momentum)};
      const double yTrueLineTle[2] = {yMinTle, yMaxTle};
      TGraph* gTrueLine = new TGraph(2, xTrue, yTrueLineTle);
      if (gTrueLine) {
        gTrueLine->SetName(Form("g_k0_joint_truep_pair%d_dau%d", pairIdx, idau + 1));
        gTrueLine->SetLineColor(kGreen + 2);
        gTrueLine->SetLineStyle(2);
        gTrueLine->SetLineWidth(2);
        mg->Add(gTrueLine, "L");
      }
    }

    gK0JointMomentumLogLMultiGraphs.push_back(mg);

    if (!rawLogLMcs.empty()) {
      std::vector<double> dThetaObs;
      std::vector<double> xOverX0Obs;
      std::vector<double> rrObs;
      const bool haveMcsSeg = BuildMcsSegments(*reco, mcsCfg, dThetaObs, xOverX0Obs, rrObs);
      const double x0CmMcs =
          (std::isfinite(mcsCfg.radiationLengthCm) && mcsCfg.radiationLengthCm > 1e-9) ? mcsCfg.radiationLengthCm : 14.0;
      double pFromRmsGeV = -1.;
      if (haveMcsSeg && dThetaObs.size() >= 2u && !xOverX0Obs.empty()) {
        double s1 = 0.0;
        double s2 = 0.0;
        for (double dth : dThetaObs) {
          s1 += dth;
          s2 += dth * dth;
        }
        const double n = static_cast<double>(dThetaObs.size());
        const double meanDt = s1 / n;
        const double var = (s2 / n) - meanDt * meanDt;
        const double rmsDt = (var > 0.0) ? std::sqrt(var) : 0.0;
        double meanSegCm = 0.0;
        for (double xov : xOverX0Obs) meanSegCm += xov * x0CmMcs;
        meanSegCm /= static_cast<double>(xOverX0Obs.size());
        if (rmsDt > 0.0 && std::isfinite(meanSegCm) && meanSegCm > 0.0) {
          pFromRmsGeV = PionMomentumGeVFromHighlandRmsAndMeanSeg(rmsDt, meanSegCm, x0CmMcs);
        }
      }

      mgMcs = new TMultiGraph();
      if (mgMcs) {
        mgMcs->SetName(Form("mg_sigmcslogl_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
        mgMcs->SetTitle(
            Form("MCS pion log-likelihood diagnostics (pair %d, dau %d, signalCode %d);Momentum [GeV/c];logL_{MCS}",
                 pairIdx, idau + 1, signalCode));
        TGraph* gNllMcs = new TGraph(static_cast<int>(pAxisDiag.size()), pAxisDiag.data(), rawLogLMcs.data());
        if (gNllMcs) {
          gNllMcs->SetName(Form("g_k0_joint_nll_mcs_pair%d_dau%d", pairIdx, idau + 1));
          gNllMcs->SetLineColor(kOrange + 1);
          gNllMcs->SetLineWidth(2);
          mgMcs->Add(gNllMcs, "L");
        }
        std::vector<double> yMcsOnly;
        yMcsOnly.reserve(rawLogLMcs.size());
        for (double y : rawLogLMcs) {
          if (std::isfinite(y)) yMcsOnly.push_back(y);
        }
        double yMinMcs = yMinTle;
        double yMaxMcs = yMaxTle;
        if (!LogLYAxisRangeFromValues(yMcsOnly, kMaxLogLSpan, 0.08, yMinMcs, yMaxMcs)) {
          yMinMcs = yMinTle;
          yMaxMcs = yMaxTle;
        }
        const double yMcsLine[2] = {yMinMcs, yMaxMcs};
        if (std::isfinite(pBestMcs) && pBestMcs > 0.) {
          const double xMcs[2] = {pBestMcs, pBestMcs};
          TGraph* gBestMcsLine = new TGraph(2, xMcs, yMcsLine);
          if (gBestMcsLine) {
            gBestMcsLine->SetName(Form("g_k0_joint_pmcs_pair%d_dau%d", pairIdx, idau + 1));
            gBestMcsLine->SetLineColor(kOrange + 7);
            gBestMcsLine->SetLineStyle(2);
            gBestMcsLine->SetLineWidth(2);
            mgMcs->Add(gBestMcsLine, "L");
          }
        }
        if (std::isfinite(pJoint) && pJoint > 0.) {
          const double xJoint[2] = {pJoint, pJoint};
          const double yJointLineMcs[2] = {yMinMcs, yMaxMcs};
          TGraph* gJointLine = new TGraph(2, xJoint, yJointLineMcs);
          if (gJointLine) {
            gJointLine->SetName(Form("g_k0_joint_pjoint_mcs_pair%d_dau%d", pairIdx, idau + 1));
            gJointLine->SetLineColor(kMagenta + 1);
            gJointLine->SetLineStyle(3);
            gJointLine->SetLineWidth(2);
            mgMcs->Add(gJointLine, "L");
          }
        }
        const AnaTrueParticlePD* trueRecoMcs =
            reco->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco->TrueObject) : nullptr;
        if (trueRecoMcs && std::isfinite(static_cast<double>(trueRecoMcs->Momentum)) && trueRecoMcs->Momentum > 0.f) {
          const double xTrue[2] = {static_cast<double>(trueRecoMcs->Momentum), static_cast<double>(trueRecoMcs->Momentum)};
          const double yTrueLineMcs[2] = {yMinMcs, yMaxMcs};
          TGraph* gTrueLine = new TGraph(2, xTrue, yTrueLineMcs);
          if (gTrueLine) {
            gTrueLine->SetName(Form("g_k0_joint_truep_mcs_pair%d_dau%d", pairIdx, idau + 1));
            gTrueLine->SetLineColor(kGreen + 2);
            gTrueLine->SetLineStyle(2);
            gTrueLine->SetLineWidth(2);
            mgMcs->Add(gTrueLine, "L");
          }
        }
        if (std::isfinite(pFromRmsGeV) && pFromRmsGeV > 0.) {
          const double xRms[2] = {pFromRmsGeV, pFromRmsGeV};
          TGraph* gPrmsLine = new TGraph(2, xRms, yMcsLine);
          if (gPrmsLine) {
            gPrmsLine->SetName(Form("g_k0_joint_prms_mcs_pair%d_dau%d", pairIdx, idau + 1));
            gPrmsLine->SetLineColor(kCyan + 2);
            gPrmsLine->SetLineStyle(7);
            gPrmsLine->SetLineWidth(2);
            mgMcs->Add(gPrmsLine, "L");
          }
        }
        gK0JointMomentumMCSLogLMultiGraphs.push_back(mgMcs);
      }

      if (haveMcsSeg) {
        TGraph* gMcsAng = new TGraph(static_cast<int>(rrObs.size()), rrObs.data(), dThetaObs.data());
        if (gMcsAng) {
          gMcsAng->SetName(Form("g_sigmcsang_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
          gMcsAng->SetTitle(
              Form("Observed MCS #Delta#theta vs RR (pair %d, dau %d, signalCode %d);Residual range [cm];#Delta#theta [rad]",
                   pairIdx, idau + 1, signalCode));
          gMcsAng->SetMarkerStyle(kFullCircle);
          gMcsAng->SetMarkerColor(kBlack);
          gMcsAng->SetLineColor(kBlack);
          gMcsAng->SetMarkerSize(0.65);
          gK0JointMomentumMCSAngleGraphs.push_back(gMcsAng);

          std::vector<double> rrDedx, dedxHits;
          if (CollectCollectionPlaneDedxVsRR(*reco, rrDedx, dedxHits)) {
            TGraph* gMcsDedx =
                new TGraph(static_cast<int>(rrDedx.size()), rrDedx.data(), dedxHits.data());
            if (gMcsDedx) {
              gMcsDedx->SetName(Form("g_sigmcsdedx_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
              gMcsDedx->SetMarkerStyle(kFullCircle);
              gMcsDedx->SetMarkerColor(kRed);
              gMcsDedx->SetLineColor(kRed);
              gMcsDedx->SetMarkerSize(0.45);
              gK0JointMomentumMCSDedxVsRRGraphs.push_back(gMcsDedx);
            }
          }

          McsAngleStats stats;
          stats.entries = static_cast<int>(dThetaObs.size());
          stats.detectorSigmaEnabled = mcsCfg.useDetectorSigma;
          if (stats.entries > 0) {
            double s1 = 0.0;
            double s2 = 0.0;
            for (double dth : dThetaObs) {
              s1 += dth;
              s2 += dth * dth;
            }
            stats.meanObs = s1 / static_cast<double>(stats.entries);
            const double var = (s2 / static_cast<double>(stats.entries)) - stats.meanObs * stats.meanObs;
            stats.rmsObs = (var > 0.0) ? std::sqrt(var) : 0.0;
          }

          const AnaTrueParticlePD* trueRecoMcs =
              reco->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco->TrueObject) : nullptr;
          const double pTrueMcs = (trueRecoMcs && trueRecoMcs->Momentum > 0.f &&
                                   std::isfinite(static_cast<double>(trueRecoMcs->Momentum)))
                                      ? static_cast<double>(trueRecoMcs->Momentum)
                                      : -1.0;
          if (pTrueMcs > 0.0 && std::isfinite(pTrueMcs) && !xOverX0Obs.empty()) {
            constexpr double kPionMassMeV = 139.57;
            constexpr double kHighlandMeV = 13.6;
            const double pMeV = pTrueMcs * 1000.0;
            const double eMeV = std::sqrt(pMeV * pMeV + kPionMassMeV * kPionMassMeV);
            const double beta = (eMeV > 0.0 && std::isfinite(eMeV)) ? (pMeV / eMeV) : -1.0;
            if (beta > 0.0 && std::isfinite(beta)) {
              double sumTheta02 = 0.0;
              size_t nTheta0 = 0u;
              for (double xov : xOverX0Obs) {
                if (!(xov > 0.0) || !std::isfinite(xov)) continue;
                double corr = 1.0 + 0.038 * std::log(std::max(xov, 1e-12));
                if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;
                double theta0 = (kHighlandMeV / (beta * pMeV)) * std::sqrt(xov) * corr;
                if (!std::isfinite(theta0) || theta0 <= 0.0) theta0 = mcsCfg.theta0FloorRad;
                theta0 = std::max(theta0, mcsCfg.theta0FloorRad);
                sumTheta02 += theta0 * theta0;
                ++nTheta0;
              }
              if (nTheta0 > 0u && sumTheta02 >= 0.0 && std::isfinite(sumTheta02)) {
                stats.rmsTrue = std::sqrt(sumTheta02 / static_cast<double>(nTheta0));
              }
            }
          }
          if (mcsCfg.useDetectorSigma && !xOverX0Obs.empty()) {
            double sumDet2 = 0.0;
            size_t nDet = 0u;
            for (double xov : xOverX0Obs) {
              if (!(xov > 0.0) || !std::isfinite(xov)) continue;
              const double segCm = xov * x0CmMcs;
              double sigmaDet2 = 0.0;
              if (std::isfinite(segCm) && segCm > 0.0) sigmaDet2 = mcsCfg.detectorSigmaA / (segCm * segCm) + mcsCfg.detectorSigmaC;
              if (!std::isfinite(sigmaDet2) || sigmaDet2 <= 0.0) sigmaDet2 = mcsCfg.detectorSigmaFloorRad * mcsCfg.detectorSigmaFloorRad;
              const double sigmaDet = std::max(std::sqrt(sigmaDet2), mcsCfg.detectorSigmaFloorRad);
              sumDet2 += sigmaDet * sigmaDet;
              ++nDet;
            }
            if (nDet > 0u && std::isfinite(sumDet2)) {
              stats.detectorSigmaApplied = true;
              stats.rmsDetUsed = std::sqrt(sumDet2 / static_cast<double>(nDet));
              if (std::isfinite(stats.rmsObs) && stats.rmsObs >= 0.0) {
                const double mcs2 = std::max(0.0, stats.rmsObs * stats.rmsObs - stats.rmsDetUsed * stats.rmsDetUsed);
                stats.rmsMcsFromObsSubDet = std::sqrt(mcs2);
              }
            }
          }
          gK0JointMomentumMCSAngleStatsByGraphName[gMcsAng->GetName()] = stats;
        }
      }
    }
  }
}

void MaybeAccumulateJointObjective2DHeatmaps(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                           size_t candidateIndex) {
  if (!candidate) return;
  AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
  if (!vertex || vertex->Particles.size() < 2) return;
  if (vertex->JointK0sMomentumUsed != 1) return;

  const int signalCode = GetSignalCodeFromCategory(candidateIndex);
  const bool ensureSignalOnly = ND::params().HasParameter("neutralKaonAnalysis.EnsureMomentumSignalOnly") &&
                                ND::params().GetParameterI("neutralKaonAnalysis.EnsureMomentumSignalOnly") == 1;
  if (ensureSignalOnly && !IsSignalLikeCode(signalCode)) return;

  const std::string key0 = MakeDiagDauKey(event, candidateIndex, 0);
  const std::string key1 = MakeDiagDauKey(event, candidateIndex, 1);
  if (sSigDedxAcceptedDauKeys.find(key0) == sSigDedxAcceptedDauKeys.end() ||
      sSigDedxAcceptedDauKeys.find(key1) == sSigDedxAcceptedDauKeys.end()) {
    return;
  }
  auto it0 = sSigDedxPairIndexByDauKey.find(key0);
  auto it1 = sSigDedxPairIndexByDauKey.find(key1);
  if (it0 == sSigDedxPairIndexByDauKey.end() || it1 == sSigDedxPairIndexByDauKey.end() ||
      it0->second != it1->second) {
    return;
  }
  const int pairIdx = it0->second;

  AnaParticlePD* reco1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* reco2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
  if (!reco1 || !reco2) return;

  pdMomReconstruction::JointK0sTwoPionGridFitConfig jcfgHeat;
  pdMomReconstruction::FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams(jcfgHeat);
  const pdMomReconstruction::PionTLEFitConfig& tle = jcfgHeat.tle;

  std::vector<double> p1v, logL1, p2v, logL2;
  if (!pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
          reco1, tle.scanLmaxCm, tle.scanStepCm, tle.minInteriorHits, tle.skipHitsFirst, tle.skipHitsLast,
          tle.dedxMinMeVcm, tle.dedxMaxMeVcm, tle.dedxPdfPathCm, p1v, logL1, tle.scanStepFineCm,
          tle.lowPMomentumRefineGeV) ||
      !pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
          reco2, tle.scanLmaxCm, tle.scanStepCm, tle.minInteriorHits, tle.skipHitsFirst, tle.skipHitsLast,
          tle.dedxMinMeVcm, tle.dedxMaxMeVcm, tle.dedxPdfPathCm, p2v, logL2, tle.scanStepFineCm,
          tle.lowPMomentumRefineGeV)) {
    return;
  }
  std::vector<double> logL1Joint = logL1;
  std::vector<double> logL2Joint = logL2;
  if (jcfgHeat.useMCS && std::isfinite(jcfgHeat.mcsWeight) && jcfgHeat.mcsWeight > 0.) {
    const pdMomReconstruction::MCSLikelihoodConfig& mcsCfg = jcfgHeat.mcs.likelihood;
    std::vector<double> logLMcs1;
    std::vector<double> logLMcs2;
    if (pdMomReconstruction::BuildPionMCSLogLikelihoodVsMomentumCurve(*reco1, p1v, mcsCfg, logLMcs1) &&
        pdMomReconstruction::BuildPionMCSLogLikelihoodVsMomentumCurve(*reco2, p2v, mcsCfg, logLMcs2) &&
        logLMcs1.size() == logL1Joint.size() && logLMcs2.size() == logL2Joint.size()) {
      for (size_t i = 0; i < logL1Joint.size(); ++i) logL1Joint[i] += jcfgHeat.mcsWeight * logLMcs1[i];
      for (size_t i = 0; i < logL2Joint.size(); ++i) logL2Joint[i] += jcfgHeat.mcsWeight * logLMcs2[i];
    }
  }

  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart =
      ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
  TVector3 dirFit1, dirFit2;
  DiagAnnihilationDaughterFitDirs(vertex, trackFitLength, trackFitDistanceFromStart, dirFit1, dirFit2);

  const double pMinGeV = jcfgHeat.pMinGeV;
  const double pMaxGeV = jcfgHeat.pMaxGeV;
  const double pStepGeV = jcfgHeat.pStepGeV;
  const double kK0sMassGeV = jcfgHeat.mK0sMassGeV;
  double sigma_m_gev = jcfgHeat.sigmaMassGeV;
  if (jcfgHeat.useEventSigmaM && std::isfinite(static_cast<double>(vertex->JointK0sSigmaMEventGeV)) &&
      vertex->JointK0sSigmaMEventGeV > 0.f) {
    sigma_m_gev = static_cast<double>(vertex->JointK0sSigmaMEventGeV);
  }
  const double penaltyScale = jcfgHeat.massPenaltyScale;

  Int_t runId = -1;
  Int_t subRunId = -1;
  Int_t evtId = -1;
  if (event.EventInfo) {
    runId = event.EventInfo->Run;
    subRunId = event.EventInfo->SubRun;
    evtId = event.EventInfo->Event;
  }

  char titleS[512];
  char titleP[512];
  std::snprintf(titleS, sizeof(titleS),
                "(run %d subrun %d evt %d) Joint objective S=logL1+logL2-pen "
                "(coarse grid, signalCode %d pair %d);p1 [GeV/c];p2 [GeV/c];S",
                static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId), signalCode, pairIdx);
  std::snprintf(titleP, sizeof(titleP),
                "(run %d subrun %d evt %d) Mass penalty 0.5*scale*((M-m_{K0s})/#sigma)^{2} "
                "(coarse grid, signalCode %d pair %d);p1 [GeV/c];p2 [GeV/c];penalty",
                static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId), signalCode, pairIdx);
  char titleLL[512];
  std::snprintf(titleLL, sizeof(titleLL),
                "(run %d subrun %d evt %d) track logL sum used in joint fit (TLE + optional MCS, no mass term) "
                "(coarse grid, signalCode %d pair %d);p1 [GeV/c];p2 [GeV/c];logL1+logL2",
                static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId), signalCode, pairIdx);

  char nameBufS[96];
  char nameBufP[96];
  char nameBufLL[96];
  std::snprintf(nameBufS, sizeof(nameBufS), "h_sigjoint2d_%d_%d", signalCode, pairIdx);
  std::snprintf(nameBufP, sizeof(nameBufP), "h_sigjointpen_%d_%d", signalCode, pairIdx);
  std::snprintf(nameBufLL, sizeof(nameBufLL), "h_sigjointlogl_%d_%d", signalCode, pairIdx);

  TH2F* hObj = nullptr;
  TH2F* hPen = nullptr;
  TH2F* hLogL = nullptr;
  if (!pdJointK0sPionMomentum::MakeJointK0sObjectiveTH2CoarsePass(
          p1v, logL1Joint, p2v, logL2Joint, dirFit1, dirFit2, pMinGeV, pMaxGeV, pStepGeV, kK0sMassGeV, sigma_m_gev,
          penaltyScale,
          nameBufS, titleS, nameBufP, titleP, nameBufLL, titleLL, hObj, hPen, hLogL)) {
    return;
  }

  gK0JointObjectiveTH2Sum.push_back(hObj);
  gK0JointObjectiveTH2Penalty.push_back(hPen);
  gK0JointObjectiveTH2TrackLogLSum.push_back(hLogL);
  double openingAngleDeg = -999.;
  if (dirFit1.Mag2() > 0. && dirFit2.Mag2() > 0.) {
    openingAngleDeg = dirFit1.Angle(dirFit2) * 57.29577951308232;
  }
  gK0JointObjectiveOpeningAngleDeg.push_back(openingAngleDeg);
  gK0JointObjectiveSigmaMGeV.push_back(sigma_m_gev);

  size_t im1 = 0;
  for (size_t i = 1; i < logL1.size(); ++i) {
    if (logL1[i] > logL1[im1]) im1 = i;
  }
  size_t im2 = 0;
  for (size_t i = 1; i < logL2.size(); ++i) {
    if (logL2[i] > logL2[im2]) im2 = i;
  }
  const double p1TleMarg = (im1 < p1v.size()) ? p1v[im1] : -1.;
  const double p2TleMarg = (im2 < p2v.size()) ? p2v[im2] : -1.;

  const AnaTrueParticlePD* true1 =
      reco1->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco1->TrueObject) : nullptr;
  const AnaTrueParticlePD* true2 =
      reco2->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco2->TrueObject) : nullptr;
  const double p1True =
      (true1 && std::isfinite(static_cast<double>(true1->Momentum)) && true1->Momentum > 0.f)
          ? static_cast<double>(true1->Momentum)
          : -1.;
  const double p2True =
      (true2 && std::isfinite(static_cast<double>(true2->Momentum)) && true2->Momentum > 0.f)
          ? static_cast<double>(true2->Momentum)
          : -1.;

  const double p1Joint =
      (std::isfinite(static_cast<double>(reco1->Momentum)) && reco1->Momentum > 0.f)
          ? static_cast<double>(reco1->Momentum)
          : -1.;
  const double p2Joint =
      (std::isfinite(static_cast<double>(reco2->Momentum)) && reco2->Momentum > 0.f)
          ? static_cast<double>(reco2->Momentum)
          : -1.;

  auto makeIntersectionMark = [](double p1, double p2, const char* name, const char* title, Color_t col, Style_t mst,
                                 Size_t msz) -> TGraph* {
    if (!std::isfinite(p1) || !std::isfinite(p2) || p1 <= 0. || p2 <= 0.) return nullptr;
    TGraph* g = new TGraph(1);
    g->SetPoint(0, p1, p2);
    g->SetName(name);
    g->SetTitle(title);
    g->SetMarkerColor(col);
    g->SetMarkerStyle(mst);
    g->SetMarkerSize(msz);
    g->SetLineWidth(0);
    return g;
  };

  char nTle[96], nTr[96], nJn[96];
  std::snprintf(nTle, sizeof(nTle), "g_sigjoint2d_mrk_tle_%d_%d", signalCode, pairIdx);
  std::snprintf(nTr, sizeof(nTr), "g_sigjoint2d_mrk_true_%d_%d", signalCode, pairIdx);
  std::snprintf(nJn, sizeof(nJn), "g_sigjoint2d_mrk_joint_%d_%d", signalCode, pairIdx);

  gK0JointObjective2DMarkerGraphs.push_back(makeIntersectionMark(
      p1TleMarg, p2TleMarg, nTle, "TLE marginal (p1 max, p2 max);p1 [GeV/c];p2 [GeV/c]", kRed + 1, kFullCircle, 1.4));
  gK0JointObjective2DMarkerGraphs.push_back(makeIntersectionMark(
      p1True, p2True, nTr, "True (p1,p2);p1 [GeV/c];p2 [GeV/c]", kGreen + 2, kFullSquare, 1.4));
  gK0JointObjective2DMarkerGraphs.push_back(makeIntersectionMark(
      p1Joint, p2Joint, nJn, "Joint fit (p1,p2);p1 [GeV/c];p2 [GeV/c]", kMagenta + 1, kStar, 1.6));

  gK0JointObjectiveBestMarker.push_back(nullptr);
}

void ResetSignalK0DiagnosticsIfNewEvent(const AnaEventB& event) {
  Int_t r = -1;
  Int_t s = -1;
  Int_t e = -1;
  if (event.EventInfo) {
    r = event.EventInfo->Run;
    s = event.EventInfo->SubRun;
    e = event.EventInfo->Event;
  }
  if (r != sSignalDiagRun || s != sSignalDiagSub || e != sSignalDiagEvt) {
    sSignalDiagTrueK0Ids.clear();
    sSigDedxAcceptedDauKeys.clear();
    sSigDedxPairIndexByDauKey.clear();
    sSignalDiagRun = r;
    sSignalDiagSub = s;
    sSignalDiagEvt = e;
  }
}

void WriteSignalPionDedxDiagnostics(OutputManager& output) {
  TTree* tree = output.GetTree();
  if (!tree) return;
  TFile* file = tree->GetCurrentFile();
  if (!file) return;

  file->cd();
  for (TH1F* h : gK0SignalDedxBiasHistograms) {
    if (h) {
      delete h;
    }
  }
  // Build one canvas per daughter with 4 panels:
  // 1) sigdedx, 2) observed MCS angle distribution, 3) siglogl(TLEFit), 4) siglogl(MCS).
  // This keeps each top-row observable aligned with its corresponding likelihood below.
  // Only canvases are persisted; source objects are kept transient.
  struct SigMomCanvasBundle {
    TMultiGraph* dedx = nullptr;
    TMultiGraph* tle = nullptr;
    TMultiGraph* mcs = nullptr;
    TGraph* mcsAngleGraph = nullptr;
    TGraph* mcsDedxVsRRGraph = nullptr;
  };
  std::map<std::tuple<int, int, int>, SigMomCanvasBundle> mgBySignalPairDau;
  for (TMultiGraph* mg : gK0SignalDedxMultiGraphs) {
    if (!mg || !mg->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(mg->GetName(), "mg_sigdedx_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].dedx = mg;
    }
  }
  for (TMultiGraph* mg : gK0JointMomentumLogLMultiGraphs) {
    if (!mg || !mg->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(mg->GetName(), "mg_siglogl_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].tle = mg;
    }
  }
  for (TMultiGraph* mg : gK0JointMomentumMCSLogLMultiGraphs) {
    if (!mg || !mg->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(mg->GetName(), "mg_sigmcslogl_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].mcs = mg;
    }
  }
  for (TGraph* g : gK0JointMomentumMCSAngleGraphs) {
    if (!g || !g->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(g->GetName(), "g_sigmcsang_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].mcsAngleGraph = g;
    }
  }
  for (TGraph* g : gK0JointMomentumMCSDedxVsRRGraphs) {
    if (!g || !g->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(g->GetName(), "g_sigmcsdedx_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].mcsDedxVsRRGraph = g;
    }
  }
  for (const auto& it : mgBySignalPairDau) {
    const int sc = std::get<0>(it.first);
    const int pj = std::get<1>(it.first);
    const int dau = std::get<2>(it.first);
    TMultiGraph* mgDedx = it.second.dedx;
    TMultiGraph* mgLogL = it.second.tle;
    TMultiGraph* mgMcs = it.second.mcs;
    TGraph* gMcsAngles = it.second.mcsAngleGraph;
    TGraph* gMcsDedxVsRR = it.second.mcsDedxVsRRGraph;
    if (!mgDedx && !mgLogL && !mgMcs && !gMcsAngles) continue;

    TCanvas* cv = new TCanvas(Form("c_sigmomdiag_%d_dau%d_%d", sc, dau, pj),
                              Form("Signal momentum diagnostics (signal %d, pair %d, dau %d)", sc, pj, dau),
                              1300, 950);
    cv->Divide(2, 2);
    cv->cd(1);
    if (mgDedx) mgDedx->Draw("A");
    cv->cd(2);
    TPad* pad2 = static_cast<TPad*>(gPad);
    if (pad2) pad2->SetRightMargin(0.26);
    if (gMcsAngles) {
      double xRRlo = 0.;
      double xRRhi = 0.;
      bool haveRR = false;
      auto expandXRr = [&haveRR, &xRRlo, &xRRhi](const TGraph* g) {
        if (!g || g->GetN() <= 0) return;
        const double* xx = g->GetX();
        for (int i = 0; i < g->GetN(); ++i) {
          if (!std::isfinite(xx[i])) continue;
          if (!haveRR) {
            xRRlo = xRRhi = xx[i];
            haveRR = true;
          } else {
            if (xx[i] < xRRlo) xRRlo = xx[i];
            if (xx[i] > xRRhi) xRRhi = xx[i];
          }
        }
      };
      expandXRr(gMcsAngles);
      expandXRr(gMcsDedxVsRR);
      if (haveRR && xRRhi > xRRlo) {
        const double mx = 0.03 * (xRRhi - xRRlo);
        gMcsAngles->GetXaxis()->SetLimits(xRRlo - mx, xRRhi + mx);
        gMcsAngles->GetXaxis()->SetRangeUser(xRRlo - mx, xRRhi + mx);
      }
      gMcsAngles->Draw("AP");
      gPad->Update();
      if (gMcsDedxVsRR && gMcsDedxVsRR->GetN() > 0) {
        constexpr double kDedxAxisMinMeVcm = 0.;
        constexpr double kDedxAxisMaxMeVcm = 10.;
        const int n2 = gMcsDedxVsRR->GetN();
        double* x2 = gMcsDedxVsRR->GetX();
        double* y2 = gMcsDedxVsRR->GetY();
        const double y1min = gPad->GetUymin();
        const double y1max = gPad->GetUymax();
        const double xrmax = gPad->GetUxmax();
        const double dSpan = kDedxAxisMaxMeVcm - kDedxAxisMinMeVcm;
        std::vector<double> yScaled(static_cast<size_t>(n2));
        for (int i = 0; i < n2; ++i) {
          const double yc =
              std::max(kDedxAxisMinMeVcm, std::min(kDedxAxisMaxMeVcm, static_cast<double>(y2[i])));
          yScaled[static_cast<size_t>(i)] =
              y1min + (yc - kDedxAxisMinMeVcm) / dSpan * (y1max - y1min);
        }
        auto* gDedxScaled = new TGraph(n2, x2, yScaled.data());
        gDedxScaled->SetMarkerColor(kRed);
        gDedxScaled->SetLineColor(kRed);
        gDedxScaled->SetMarkerStyle(kFullCircle);
        gDedxScaled->SetMarkerSize(0.45);
        gDedxScaled->Draw("P SAME");
        // Heap-allocated: a stack TGaxis is destroyed before the canvas is written, so the axis never persists.
        auto* axisDedx =
            new TGaxis(xrmax, y1min, xrmax, y1max, kDedxAxisMinMeVcm, kDedxAxisMaxMeVcm, 510, "+L");
        axisDedx->SetLineColor(kRed);
        axisDedx->SetLabelColor(kRed);
        axisDedx->SetTitleColor(kRed);
        axisDedx->SetLabelSize(0.045);
        axisDedx->SetTitleSize(0.048);
        axisDedx->SetTitleOffset(1.15);
        axisDedx->SetTitle("dE/dx [MeV/cm]");
        axisDedx->Draw();
        gPad->Modified();
        gPad->Update();
      }
      TLegend* leg = new TLegend(0.50, 0.58, 0.88, 0.89);
      if (leg) {
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gMcsAngles, "Observed #Delta#theta vs RR", "p");
        if (gMcsDedxVsRR && gMcsDedxVsRR->GetN() > 0) {
          leg->AddEntry(gMcsDedxVsRR, "dE/dx vs RR (collection)", "p");
        }
        McsAngleStats stats;
        auto itStats = gK0JointMomentumMCSAngleStatsByGraphName.find(gMcsAngles->GetName());
        if (itStats != gK0JointMomentumMCSAngleStatsByGraphName.end()) stats = itStats->second;
        leg->AddEntry((TObject*)nullptr, Form("N=%d", stats.entries), "");
        if (std::isfinite(stats.meanObs)) {
          leg->AddEntry((TObject*)nullptr, Form("Mean(obs)=%.4f rad", stats.meanObs), "");
        } else {
          leg->AddEntry((TObject*)nullptr, "Mean(obs)=n/a", "");
        }
        if (std::isfinite(stats.rmsObs)) {
          leg->AddEntry((TObject*)nullptr, Form("RMS(obs)=%.4f rad", stats.rmsObs), "");
        } else {
          leg->AddEntry((TObject*)nullptr, "RMS(obs)=n/a", "");
        }
        if (std::isfinite(stats.rmsTrue) && stats.rmsTrue > 0.0) {
          leg->AddEntry((TObject*)nullptr, Form("RMS(exp @ p_{true})=%.4f rad", stats.rmsTrue), "");
        } else {
          leg->AddEntry((TObject*)nullptr, "RMS(exp @ p_{true})=n/a", "");
        }
        if (stats.detectorSigmaEnabled) {
          if (stats.detectorSigmaApplied && std::isfinite(stats.rmsDetUsed) && stats.rmsDetUsed >= 0.0) {
            leg->AddEntry((TObject*)nullptr, Form("RMS(det used)=%.4f rad", stats.rmsDetUsed), "");
          } else {
            leg->AddEntry((TObject*)nullptr, "RMS(det used)=n/a", "");
          }
          if (std::isfinite(stats.rmsMcsFromObsSubDet) && stats.rmsMcsFromObsSubDet >= 0.0) {
            leg->AddEntry((TObject*)nullptr, Form("RMS(MCS from obs-det)=%.4f rad", stats.rmsMcsFromObsSubDet), "");
          } else {
            leg->AddEntry((TObject*)nullptr, "RMS(MCS from obs-det)=n/a", "");
          }
        } else {
          leg->AddEntry((TObject*)nullptr, "RMS(det used)=off", "");
        }
        leg->Draw();
      }
    }
    cv->cd(3);
    if (mgLogL) {
      mgLogL->Draw("A");
      ClipLogLMultigraphYAxisAfterDraw(mgLogL);
      AddLogLPadLegend(mgLogL, false);
    }
    cv->cd(4);
    if (mgMcs) {
      mgMcs->Draw("A");
      ClipLogLMultigraphYAxisAfterDraw(mgMcs);
      AddLogLPadLegend(mgMcs, true);
    }
    file->cd();
    cv->Write(cv->GetName(), TObject::kOverwrite);
    delete cv;
  }
  for (TMultiGraph* mg : gK0SignalDedxMultiGraphs) {
    DeleteMultiGraphAndGraphs(mg);
  }
  for (TMultiGraph* mg : gK0JointMomentumLogLMultiGraphs) {
    DeleteMultiGraphAndGraphs(mg);
  }
  for (TMultiGraph* mg : gK0JointMomentumMCSLogLMultiGraphs) {
    DeleteMultiGraphAndGraphs(mg);
  }
  for (TGraph* g : gK0JointMomentumMCSAngleGraphs) {
    if (g) delete g;
  }
  for (TGraph* g : gK0JointMomentumMCSDedxVsRRGraphs) {
    if (g) delete g;
  }

  const size_t nJoint2d = gK0JointObjectiveTH2Sum.size();
  for (size_t ij = 0; ij < nJoint2d; ++ij) {
    TH2F* hS = gK0JointObjectiveTH2Sum[ij];
    TH2F* hP = (ij < gK0JointObjectiveTH2Penalty.size()) ? gK0JointObjectiveTH2Penalty[ij] : nullptr;
    TH2F* hLL = (ij < gK0JointObjectiveTH2TrackLogLSum.size()) ? gK0JointObjectiveTH2TrackLogLSum[ij] : nullptr;
    TGraph* bestM = (ij < gK0JointObjectiveBestMarker.size()) ? gK0JointObjectiveBestMarker[ij] : nullptr;
    TGraph* m0 = (3 * ij + 0 < gK0JointObjective2DMarkerGraphs.size()) ? gK0JointObjective2DMarkerGraphs[3 * ij + 0] : nullptr;
    TGraph* m1 = (3 * ij + 1 < gK0JointObjective2DMarkerGraphs.size()) ? gK0JointObjective2DMarkerGraphs[3 * ij + 1] : nullptr;
    TGraph* m2 = (3 * ij + 2 < gK0JointObjective2DMarkerGraphs.size()) ? gK0JointObjective2DMarkerGraphs[3 * ij + 2] : nullptr;

    int sc = -1;
    int pj = -1;
    if (hS && std::sscanf(hS->GetName(), "h_sigjoint2d_%d_%d", &sc, &pj) == 2 && sc >= 0 && pj >= 0) {
      const double openingAngleDeg =
          (ij < gK0JointObjectiveOpeningAngleDeg.size()) ? gK0JointObjectiveOpeningAngleDeg[ij] : -999.;
      const double sigmaMGeV =
          (ij < gK0JointObjectiveSigmaMGeV.size()) ? gK0JointObjectiveSigmaMGeV[ij] : -999.;
      const double sigmaMMeV = (std::isfinite(sigmaMGeV) && sigmaMGeV > 0.) ? sigmaMGeV * 1e3 : -999.;
      TCanvas* cv2d = new TCanvas(Form("c_sigjoint_bundle_%d_%d", sc, pj),
                                  Form("Joint 2D diagnostics (signal %d pair %d, theta12=%.2f deg, sigma_m=%.2f MeV)",
                                       sc, pj, openingAngleDeg, sigmaMMeV),
                                  1800, 550);
      cv2d->Divide(3, 1);
      TPad* p1 = static_cast<TPad*>(cv2d->cd(1));
      if (p1) p1->SetRightMargin(0.18);
      if (hS) hS->Draw("COLZ");
      if (m0) m0->Draw("P SAME");
      if (m1) m1->Draw("P SAME");
      if (m2) m2->Draw("P SAME");
      TPad* p2 = static_cast<TPad*>(cv2d->cd(2));
      if (p2) p2->SetRightMargin(0.18);
      if (hP) hP->Draw("COLZ");
      if (m0) m0->Draw("P SAME");
      if (m1) m1->Draw("P SAME");
      if (m2) m2->Draw("P SAME");
      TPad* p3 = static_cast<TPad*>(cv2d->cd(3));
      if (p3) p3->SetRightMargin(0.18);
      if (hLL) hLL->Draw("COLZ");
      if (m0) m0->Draw("P SAME");
      if (m1) m1->Draw("P SAME");
      if (m2) m2->Draw("P SAME");
      file->cd();
      cv2d->Write(cv2d->GetName(), TObject::kOverwrite);
      delete cv2d;
    }

    if (hS) delete hS;
    if (hP) delete hP;
    if (hLL) delete hLL;
    if (bestM) delete bestM;
    for (TGraph* g : {m0, m1, m2}) {
      if (g) delete g;
    }
  }
  gK0SignalDedxMultiGraphs.clear();
  gK0SignalDedxBiasHistograms.clear();
  gK0JointMomentumLogLMultiGraphs.clear();
  gK0JointMomentumMCSLogLMultiGraphs.clear();
  gK0JointMomentumMCSAngleGraphs.clear();
  gK0JointMomentumMCSDedxVsRRGraphs.clear();
  gK0JointMomentumMCSAngleStatsByGraphName.clear();
  gK0JointObjectiveTH2Sum.clear();
  gK0JointObjectiveTH2Penalty.clear();
  gK0JointObjectiveTH2TrackLogLSum.clear();
  gK0JointObjectiveOpeningAngleDeg.clear();
  gK0JointObjectiveSigmaMGeV.clear();
  gK0JointObjectiveBestMarker.clear();
  gK0JointObjective2DMarkerGraphs.clear();
  sSigDedxAcceptedDauKeys.clear();
  sSigDedxPairIndexByDauKey.clear();
  sK0SigDedxSerialByCode.clear();
}

} // namespace neutralKaonTreeDiagnostics
