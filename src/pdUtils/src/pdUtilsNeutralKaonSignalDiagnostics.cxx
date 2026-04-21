#include "pdUtilsNeutralKaonSignalDiagnostics.hxx"

#include "CategoryManager.hxx"
#include "Parameters.hxx"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "neutralKaonAnalysisUtils.hxx"
#include "pdJointK0sPionMomentum.hxx"
#include "pdUtilsDEdx.hxx"
#include "pdUtilsLineFit.hxx"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <string>
#include <tuple>
#include <unordered_set>
#include <TH2F.h>
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

  bool IsSignalLikeCode(int signalCode) {
    return (signalCode == 1 || signalCode == 5 || signalCode == 6);
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

    std::vector<double> pAxis;
    std::vector<double> logL;
    if (!pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(
            reco, Lmax, step, minInterior, skipFirst, skipLast, dedxMin, dedxMax, pdfPath, pAxis, logL)) {
      continue;
    }
    if (pAxis.empty() || pAxis.size() != logL.size()) continue;

    auto itBest = std::max_element(logL.begin(), logL.end());
    if (itBest == logL.end()) continue;
    const size_t bestIdx = static_cast<size_t>(std::distance(logL.begin(), itBest));
    if (bestIdx >= pAxis.size()) continue;
    const double pBest = pAxis[bestIdx];
    const double yBest = logL[bestIdx];

    double yMin = yBest - 1.0;
    double yMax = yBest + 1.0;
    bool firstFinite = true;
    for (double y : logL) {
      if (!std::isfinite(y)) continue;
      if (firstFinite) {
        yMin = y;
        yMax = y;
        firstFinite = false;
      } else {
        yMin = std::min(yMin, y);
        yMax = std::max(yMax, y);
      }
    }
    if (!std::isfinite(yMin) || !std::isfinite(yMax) || yMin == yMax) {
      yMin = yBest - 1.0;
      yMax = yBest + 1.0;
    }

    TGraph* gLogL = new TGraph(static_cast<int>(pAxis.size()), pAxis.data(), logL.data());
    if (!gLogL) continue;
    const double xBest[2] = {pBest, pBest};
    const double yBestLine[2] = {yMin, yMax};
    TGraph* gBestLine = new TGraph(2, xBest, yBestLine);
    if (!gBestLine) {
      delete gLogL;
      continue;
    }
    TMultiGraph* mg = new TMultiGraph();
    if (!mg) {
      delete gLogL;
      delete gBestLine;
      continue;
    }

    mg->SetName(Form("mg_siglogl_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    mg->SetTitle(Form("Joint pion logL (pair %d, dau %d, signalCode %d);Momentum [GeV/c];log L",
                      pairIdx, idau + 1, signalCode));
    gLogL->SetName(Form("g_k0_joint_logl_pair%d_dau%d", pairIdx, idau + 1));
    gLogL->SetLineColor(kBlue + 1);
    gLogL->SetLineWidth(2);
    gBestLine->SetName(Form("g_k0_joint_bestp_pair%d_dau%d", pairIdx, idau + 1));
    gBestLine->SetLineColor(kRed + 1);
    gBestLine->SetLineWidth(2);
    mg->Add(gLogL, "L");
    mg->Add(gBestLine, "L");

    if (std::isfinite(static_cast<double>(reco->Momentum)) && reco->Momentum > 0.f) {
      const double xReco[2] = {static_cast<double>(reco->Momentum), static_cast<double>(reco->Momentum)};
      const double yRecoLine[2] = {yMin, yMax};
      TGraph* gRecoLine = new TGraph(2, xReco, yRecoLine);
      if (gRecoLine) {
        gRecoLine->SetName(Form("g_k0_joint_recop_pair%d_dau%d", pairIdx, idau + 1));
        gRecoLine->SetLineColor(kMagenta + 1);
        gRecoLine->SetLineStyle(3);
        gRecoLine->SetLineWidth(2);
        mg->Add(gRecoLine, "L");
      }
    }

    const AnaTrueParticlePD* trueReco =
        reco->TrueObject ? static_cast<const AnaTrueParticlePD*>(reco->TrueObject) : nullptr;
    if (trueReco && std::isfinite(static_cast<double>(trueReco->Momentum)) && trueReco->Momentum > 0.f) {
      const double xTrue[2] = {static_cast<double>(trueReco->Momentum), static_cast<double>(trueReco->Momentum)};
      const double yTrueLine[2] = {yMin, yMax};
      TGraph* gTrueLine = new TGraph(2, xTrue, yTrueLine);
      if (gTrueLine) {
        gTrueLine->SetName(Form("g_k0_joint_truep_pair%d_dau%d", pairIdx, idau + 1));
        gTrueLine->SetLineColor(kGreen + 2);
        gTrueLine->SetLineStyle(2);
        gTrueLine->SetLineWidth(2);
        mg->Add(gTrueLine, "L");
      }
    }

    gK0JointMomentumLogLMultiGraphs.push_back(mg);
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

  std::vector<double> p1v, logL1, p2v, logL2;
  if (!pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(reco1, Lmax, step, minInterior, skipFirst, skipLast,
                                                                  dedxMin, dedxMax, pdfPath, p1v, logL1) ||
      !pdAnaUtils::BuildPionFreeRangeLogLikelihoodVsMomentumCurve(reco2, Lmax, step, minInterior, skipFirst, skipLast,
                                                                  dedxMin, dedxMax, pdfPath, p2v, logL2)) {
    return;
  }

  const double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
  const double trackFitDistanceFromStart =
      ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");
  TVector3 dirFit1, dirFit2;
  DiagAnnihilationDaughterFitDirs(vertex, trackFitLength, trackFitDistanceFromStart, dirFit1, dirFit2);

  const double pMinGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                             ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMinGeV")
                             : 0.05;
  const double pMaxGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                             ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPMaxGeV")
                             : 2.0;
  const double pStepGeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                              ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMomentumPStepGeV")
                              : 0.05;
  constexpr double kK0sMassGeV = 0.497611;
  const double sigmaMassMeV = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                  ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassSigmaMeV")
                                  : 10.0;
  double sigma_m_gev = sigmaMassMeV * 1e-3;
  const bool useEventSigmaM =
      !ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") ||
      ND::params().GetParameterI("neutralKaonAnalysis.JointK0sMassSigmaEventPropagation") != 0;
  if (useEventSigmaM && std::isfinite(static_cast<double>(vertex->JointK0sSigmaMEventGeV)) &&
      vertex->JointK0sSigmaMEventGeV > 0.f) {
    sigma_m_gev = static_cast<double>(vertex->JointK0sSigmaMEventGeV);
  }
  const double penaltyScale = ND::params().HasParameter("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                                  ? ND::params().GetParameterD("neutralKaonAnalysis.JointK0sMassPenaltyScale")
                                  : 1.0;

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
                "(run %d subrun %d evt %d) logL1+logL2 free-range TLE (no mass term) "
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
          p1v, logL1, p2v, logL2, dirFit1, dirFit2, pMinGeV, pMaxGeV, pStepGeV, kK0sMassGeV, sigma_m_gev, penaltyScale,
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
  // Build one canvas per daughter with sigdedx + siglogl stacked.
  // Only canvases are persisted; source objects are kept transient.
  std::map<std::tuple<int, int, int>, std::pair<TMultiGraph*, TMultiGraph*> > mgBySignalPairDau;
  for (TMultiGraph* mg : gK0SignalDedxMultiGraphs) {
    if (!mg || !mg->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(mg->GetName(), "mg_sigdedx_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].first = mg;
    }
  }
  for (TMultiGraph* mg : gK0JointMomentumLogLMultiGraphs) {
    if (!mg || !mg->GetName()) continue;
    int sc = -1;
    int dau = -1;
    int pj = -1;
    if (std::sscanf(mg->GetName(), "mg_siglogl_%d_dau%d_%d", &sc, &dau, &pj) == 3) {
      mgBySignalPairDau[std::make_tuple(sc, pj, dau)].second = mg;
    }
  }
  for (const auto& it : mgBySignalPairDau) {
    const int sc = std::get<0>(it.first);
    const int pj = std::get<1>(it.first);
    const int dau = std::get<2>(it.first);
    TMultiGraph* mgDedx = it.second.first;
    TMultiGraph* mgLogL = it.second.second;
    if (!mgDedx && !mgLogL) continue;

    TCanvas* cv = new TCanvas(Form("c_sigmomdiag_%d_dau%d_%d", sc, dau, pj),
                              Form("Signal momentum diagnostics (signal %d, pair %d, dau %d)", sc, pj, dau),
                              900, 900);
    cv->Divide(1, 2);
    cv->cd(1);
    if (mgDedx) mgDedx->Draw("A");
    cv->cd(2);
    if (mgLogL) mgLogL->Draw("A");
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
