#include "pionMomentumAnalysisUtils.hxx"

#include "CategoryManager.hxx"
#include "pdUtilsDEdx.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TTree.h"
#include "TGaxis.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace pionMomentumAnaUtils {

namespace {
  struct DiagBundle {
    int signalCode = -999;
    int serialInCode = -1;
    TMultiGraph* dedx = nullptr;
    TGraph* mcsAngles = nullptr;
    TGraph* mcsAnglesKept = nullptr;
    TGraph* mcsDedxVsRR = nullptr;
    TMultiGraph* tleLogL = nullptr;
    TMultiGraph* mcsLogL = nullptr;
    int mcsN = 0;
    int mcsNKept = 0;
    double mcsMeanObs = -999.;
    double mcsRmsObs = -999.;
    double mcsMeanObsAll = -999.;
    double mcsRmsObsAll = -999.;
    double mcsRmsTrue = -999.;
  };

  std::map<std::tuple<int, int, int>, DiagBundle> sDiagByEvent;
  std::map<int, int> sDiagSerialByCode;

  void DeleteMultiGraphAndGraphs(TMultiGraph* mg) {
    if (!mg) return;
    TList* gl = mg->GetListOfGraphs();
    if (gl) gl->Delete();
    delete mg;
  }

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
}

//********************************************************************
void AddCustomCategories() {
//********************************************************************
  if (!anaUtils::_categ || anaUtils::_categ->HasCategory("maintruestopping")) return;

  std::string part_types[] = {"stopping", "non_stopping", "not_true_pion", NAMEOTHER};
  int part_codes[] = {1, 2, 3, CATOTHER};
  int part_colors[] = {2, 46, 7, COLOTHER};
  const int npart = sizeof(part_types) / sizeof(part_types[0]);

  std::reverse(part_types, part_types + npart);
  std::reverse(part_codes, part_codes + npart);
  std::reverse(part_colors, part_colors + npart);

  anaUtils::_categ->AddCategory("maintruestopping", npart, part_types, part_codes, part_colors);
}

//********************************************************************
int FillMainTrueStoppingCategory(AnaParticlePD* mainTrack, double stoppingMaxTrueEndMomentumGeV) {
//********************************************************************
  int stopCode = CATNOTRUTH;
  if (mainTrack && mainTrack->TrueObject) {
    const AnaTrueParticlePD* tpart = static_cast<const AnaTrueParticlePD*>(mainTrack->TrueObject);
    if (tpart) {
      if (std::abs(tpart->PDG) != 211) {
        stopCode = 3;
      } else {
        const double pend = static_cast<double>(tpart->MomentumEnd);
        if (std::isfinite(pend) && pend >= 0.) {
          stopCode = (pend <= stoppingMaxTrueEndMomentumGeV) ? 1 : 2;
        }
      }
    }
  }

  if (anaUtils::_categ && anaUtils::_categ->HasCategory("maintruestopping")) {
    anaUtils::_categ->SetCode("maintruestopping", stopCode, CATOTHER);
  }
  return stopCode;
}

//********************************************************************
void MaybeAccumulateMainTrackMomentumDiagnostics(AnaParticlePD* mainTrack, const AnaEventB& event, int stoppingCode,
                                                 const MomentumDiagConfig& cfg, double bestTleMomentumGeV,
                                                 double bestMcsMomentumGeV, const std::vector<double>& pAxisTleGeV,
                                                 const std::vector<double>& logLTle,
                                                 const std::vector<double>& pAxisMcsGeV,
                                                 const std::vector<double>& logLMcs,
                                                 const std::vector<double>& mcsDeltaTheta,
                                                 const std::vector<double>& mcsSegmentLengthCm,
                                                 const std::vector<double>& mcsDeltaThetaKept,
                                                 const std::vector<double>& mcsResidualRangeKeptCm) {
//********************************************************************
  if (!cfg.enableMomentumDiagnosticMultigraphs || !mainTrack) return;
  if (cfg.ensureMomentumSignalOnly && !(stoppingCode == 1 || stoppingCode == 2)) return;
  if (pAxisTleGeV.empty() || pAxisTleGeV.size() != logLTle.size()) return;

  Int_t runId = -1, subRunId = -1, evtId = -1;
  if (event.EventInfo) {
    runId = event.EventInfo->Run;
    subRunId = event.EventInfo->SubRun;
    evtId = event.EventInfo->Event;
  }
  const auto key = std::make_tuple(static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId));
  if (sDiagByEvent.find(key) != sDiagByEvent.end()) return;

  DiagBundle bundle;
  bundle.signalCode = stoppingCode;
  bundle.serialInCode = sDiagSerialByCode[stoppingCode]++;

  char dedxXAxis[200];
  std::snprintf(dedxXAxis, sizeof(dedxXAxis), "(run %d, subrun %d, evt %d) Residual range [cm]", static_cast<int>(runId),
                static_cast<int>(subRunId), static_cast<int>(evtId));
  bundle.dedx = pdAnaUtils::MakePionFreeRangeDedxVsRRMultiGraph(
      mainTrack, cfg.freeRangeScanLmaxCm, cfg.freeRangeScanStepCm, cfg.freeRangeDedxSkipHitsFirst,
      cfg.freeRangeDedxSkipHitsLast, cfg.freeRangeDedxMinMeVcm, cfg.freeRangeDedxMaxMeVcm, cfg.freeRangeDedxMinInteriorHits,
      cfg.freeRangeDedxPdfPathCm, dedxXAxis);

  bundle.tleLogL = new TMultiGraph();
  if (bundle.tleLogL) {
    bundle.tleLogL->SetName(Form("mg_main_siglogl_%d_%d", stoppingCode, bundle.serialInCode));
    bundle.tleLogL->SetTitle(Form("Main-track TLE log-likelihood (category code %d);Momentum [GeV/c];logL_{TLE}",
                                  stoppingCode));
    TGraph* gTle = new TGraph(static_cast<int>(pAxisTleGeV.size()), pAxisTleGeV.data(), logLTle.data());
    gTle->SetLineColor(kBlack);
    gTle->SetLineWidth(2);
    bundle.tleLogL->Add(gTle, "L");
    if (std::isfinite(bestTleMomentumGeV) && bestTleMomentumGeV > 0.) {
      const double x[2] = {bestTleMomentumGeV, bestTleMomentumGeV};
      const double y[2] = {*std::min_element(logLTle.begin(), logLTle.end()), *std::max_element(logLTle.begin(), logLTle.end())};
      TGraph* gBest = new TGraph(2, x, y);
      gBest->SetLineColor(kBlue + 1);
      gBest->SetLineStyle(2);
      gBest->SetLineWidth(2);
      bundle.tleLogL->Add(gBest, "L");
    }
    const AnaTrueParticlePD* truePart = mainTrack->TrueObject ? static_cast<const AnaTrueParticlePD*>(mainTrack->TrueObject) : nullptr;
    if (truePart && std::isfinite(static_cast<double>(truePart->Momentum)) && truePart->Momentum > 0.f) {
      const double pTrue = static_cast<double>(truePart->Momentum);
      const double x[2] = {pTrue, pTrue};
      const double y[2] = {*std::min_element(logLTle.begin(), logLTle.end()), *std::max_element(logLTle.begin(), logLTle.end())};
      TGraph* gTrue = new TGraph(2, x, y);
      gTrue->SetLineColor(kGreen + 2);
      gTrue->SetLineStyle(7);
      gTrue->SetLineWidth(2);
      bundle.tleLogL->Add(gTrue, "L");
    }
  }

  if (!pAxisMcsGeV.empty() && pAxisMcsGeV.size() == logLMcs.size()) {
    bundle.mcsLogL = new TMultiGraph();
    if (bundle.mcsLogL) {
      bundle.mcsLogL->SetName(Form("mg_main_sigmcslogl_%d_%d", stoppingCode, bundle.serialInCode));
      bundle.mcsLogL->SetTitle(Form("Main-track MCS log-likelihood (category code %d);Momentum [GeV/c];logL_{MCS}",
                                    stoppingCode));
      TGraph* gMcs = new TGraph(static_cast<int>(pAxisMcsGeV.size()), pAxisMcsGeV.data(), logLMcs.data());
      gMcs->SetLineColor(kOrange + 1);
      gMcs->SetLineWidth(2);
      bundle.mcsLogL->Add(gMcs, "L");
      if (std::isfinite(bestMcsMomentumGeV) && bestMcsMomentumGeV > 0.) {
        const double x[2] = {bestMcsMomentumGeV, bestMcsMomentumGeV};
        const double y[2] = {*std::min_element(logLMcs.begin(), logLMcs.end()), *std::max_element(logLMcs.begin(), logLMcs.end())};
        TGraph* gBest = new TGraph(2, x, y);
        gBest->SetLineColor(kOrange + 7);
        gBest->SetLineStyle(2);
        gBest->SetLineWidth(2);
        bundle.mcsLogL->Add(gBest, "L");
      }
      const AnaTrueParticlePD* truePart = mainTrack->TrueObject ? static_cast<const AnaTrueParticlePD*>(mainTrack->TrueObject) : nullptr;
      if (truePart && std::isfinite(static_cast<double>(truePart->Momentum)) && truePart->Momentum > 0.f) {
        const double pTrue = static_cast<double>(truePart->Momentum);
        const double x[2] = {pTrue, pTrue};
        const double y[2] = {*std::min_element(logLMcs.begin(), logLMcs.end()), *std::max_element(logLMcs.begin(), logLMcs.end())};
        TGraph* gTrue = new TGraph(2, x, y);
        gTrue->SetLineColor(kGreen + 2);
        gTrue->SetLineStyle(7);
        gTrue->SetLineWidth(2);
        bundle.mcsLogL->Add(gTrue, "L");
      }
    }
  }

  if (!mcsDeltaTheta.empty() && mcsDeltaTheta.size() == mcsSegmentLengthCm.size()) {
    std::vector<double> rrMcs;
    rrMcs.reserve(mcsSegmentLengthCm.size());
    double totalLen = 0.0;
    for (double seg : mcsSegmentLengthCm) totalLen += seg;
    double accum = 0.0;
    for (double seg : mcsSegmentLengthCm) {
      const double segCenter = accum + 0.5 * seg;
      rrMcs.push_back(std::max(0.0, totalLen - segCenter));
      accum += seg;
    }
    bundle.mcsAngles = new TGraph(static_cast<int>(rrMcs.size()), rrMcs.data(), mcsDeltaTheta.data());
    if (bundle.mcsAngles) {
      bundle.mcsAngles->SetName(Form("g_main_sigmcsang_%d_%d", stoppingCode, bundle.serialInCode));
      bundle.mcsAngles->SetTitle(
          Form("Observed MCS #Delta#theta vs residual range (category code %d);Residual range [cm];#Delta#theta [rad]",
               stoppingCode));
      bundle.mcsAngles->SetMarkerStyle(kFullCircle);
      bundle.mcsAngles->SetMarkerSize(0.65);
      bundle.mcsAngles->SetLineColor(kBlack);
      bundle.mcsAngles->SetMarkerColor(kBlack);
    }
    bundle.mcsN = static_cast<int>(mcsDeltaTheta.size());
    if (bundle.mcsN > 0) {
      double s1 = 0.0;
      double s2 = 0.0;
      for (double dth : mcsDeltaTheta) {
        s1 += dth;
        s2 += dth * dth;
      }
      bundle.mcsMeanObsAll = s1 / static_cast<double>(bundle.mcsN);
      const double varAll =
          (s2 / static_cast<double>(bundle.mcsN)) - bundle.mcsMeanObsAll * bundle.mcsMeanObsAll;
      bundle.mcsRmsObsAll = (varAll > 0.) ? std::sqrt(varAll) : 0.0;
    }

    if (!mcsDeltaThetaKept.empty() && mcsDeltaThetaKept.size() == mcsResidualRangeKeptCm.size()) {
      bundle.mcsNKept = static_cast<int>(mcsDeltaThetaKept.size());
      if (bundle.mcsNKept > 0) {
        double s1k = 0.0;
        double s2k = 0.0;
        for (double dth : mcsDeltaThetaKept) {
          s1k += dth;
          s2k += dth * dth;
        }
        bundle.mcsMeanObs = s1k / static_cast<double>(bundle.mcsNKept);
        const double varK =
            (s2k / static_cast<double>(bundle.mcsNKept)) - bundle.mcsMeanObs * bundle.mcsMeanObs;
        bundle.mcsRmsObs = (varK > 0.) ? std::sqrt(varK) : 0.0;
      }
      bundle.mcsAnglesKept =
          new TGraph(static_cast<int>(mcsResidualRangeKeptCm.size()), mcsResidualRangeKeptCm.data(),
                     mcsDeltaThetaKept.data());
      if (bundle.mcsAnglesKept) {
        bundle.mcsAnglesKept->SetName(Form("g_main_sigmcsangkept_%d_%d", stoppingCode, bundle.serialInCode));
        bundle.mcsAnglesKept->SetMarkerStyle(kFullCircle);
        bundle.mcsAnglesKept->SetMarkerSize(0.75);
        bundle.mcsAnglesKept->SetLineColor(kGreen + 2);
        bundle.mcsAnglesKept->SetMarkerColor(kGreen + 2);
      }
    } else {
      bundle.mcsNKept = bundle.mcsN;
      bundle.mcsMeanObs = bundle.mcsMeanObsAll;
      bundle.mcsRmsObs = bundle.mcsRmsObsAll;
    }

    const AnaTrueParticlePD* truePart =
        mainTrack->TrueObject ? static_cast<const AnaTrueParticlePD*>(mainTrack->TrueObject) : nullptr;
    if (truePart && std::isfinite(static_cast<double>(truePart->Momentum)) && truePart->Momentum > 0.f &&
        std::isfinite(cfg.mcsRadiationLengthCm) && cfg.mcsRadiationLengthCm > 0.) {
      constexpr double kPionMassMeV = 139.57;
      constexpr double kHighlandMeV = 13.6;
      const double pMeV = static_cast<double>(truePart->Momentum) * 1000.0;
      const double eMeV = std::sqrt(pMeV * pMeV + kPionMassMeV * kPionMassMeV);
      const double beta = (eMeV > 0.0 && std::isfinite(eMeV)) ? (pMeV / eMeV) : -1.0;
      if (beta > 0.0 && std::isfinite(beta)) {
        double sumTheta02 = 0.0;
        int nTheta0 = 0;
        for (double segCm : mcsSegmentLengthCm) {
          if (!(segCm > 0.0) || !std::isfinite(segCm)) continue;
          const double xOverX0 = segCm / cfg.mcsRadiationLengthCm;
          if (!(xOverX0 > 0.0) || !std::isfinite(xOverX0)) continue;
          double corr = 1.0 + 0.038 * std::log(std::max(xOverX0, 1e-12));
          if (!std::isfinite(corr) || corr < 0.1) corr = 0.1;
          const double theta0 = (kHighlandMeV / (beta * pMeV)) * std::sqrt(xOverX0) * corr;
          if (!std::isfinite(theta0) || theta0 <= 0.0) continue;
          sumTheta02 += theta0 * theta0;
          ++nTheta0;
        }
        if (nTheta0 > 0 && std::isfinite(sumTheta02)) {
          bundle.mcsRmsTrue = std::sqrt(sumTheta02 / static_cast<double>(nTheta0));
        }
      }
    }
  }

  std::vector<double> rrDedx, dedx;
  if (CollectCollectionPlaneDedxVsRR(*mainTrack, rrDedx, dedx)) {
    bundle.mcsDedxVsRR = new TGraph(static_cast<int>(rrDedx.size()), rrDedx.data(), dedx.data());
    if (bundle.mcsDedxVsRR) {
      bundle.mcsDedxVsRR->SetName(Form("g_main_sigmcsdedx_%d_%d", stoppingCode, bundle.serialInCode));
      bundle.mcsDedxVsRR->SetMarkerStyle(kFullCircle);
      bundle.mcsDedxVsRR->SetMarkerSize(0.45);
      bundle.mcsDedxVsRR->SetMarkerColor(kRed);
      bundle.mcsDedxVsRR->SetLineColor(kRed);
    }
  }

  sDiagByEvent[key] = bundle;
}

//********************************************************************
void WriteMomentumDiagnostics(OutputManager& output) {
//********************************************************************
  TTree* tree = output.GetTree();
  if (!tree) return;
  TFile* file = tree->GetCurrentFile();
  if (!file) return;
  file->cd();

  for (auto& it : sDiagByEvent) {
    DiagBundle& b = it.second;
    if (!b.dedx && !b.mcsAngles && !b.tleLogL && !b.mcsLogL) continue;

    TCanvas* cv = new TCanvas(Form("c_mainmomdiag_%d_%d", b.signalCode, b.serialInCode),
                              Form("Main momentum diagnostics (category code %d, index %d)", b.signalCode,
                                   b.serialInCode),
                              1300, 950);
    cv->Divide(2, 2);
    cv->cd(1);
    if (b.dedx) b.dedx->Draw("A");
    cv->cd(2);
    if (b.mcsAngles) {
      TPad* pad2 = static_cast<TPad*>(gPad);
      if (pad2) pad2->SetRightMargin(0.26);
      b.mcsAngles->Draw("AP");
      if (b.mcsAnglesKept && b.mcsAnglesKept->GetN() > 0) {
        b.mcsAnglesKept->Draw("P SAME");
      }
      gPad->Update();
      if (b.mcsDedxVsRR && b.mcsDedxVsRR->GetN() > 0) {
        const int n2 = b.mcsDedxVsRR->GetN();
        double* x2 = b.mcsDedxVsRR->GetX();
        double* y2 = b.mcsDedxVsRR->GetY();
        const double y1min = gPad->GetUymin(), y1max = gPad->GetUymax(), xrmax = gPad->GetUxmax();
        constexpr double kDedxMin = 0., kDedxMax = 10.;
        std::vector<double> ys(static_cast<size_t>(n2));
        for (int i = 0; i < n2; ++i) {
          const double yc = std::max(kDedxMin, std::min(kDedxMax, y2[i]));
          ys[static_cast<size_t>(i)] = y1min + (yc - kDedxMin) / (kDedxMax - kDedxMin) * (y1max - y1min);
        }
        TGraph* gScaled = new TGraph(n2, x2, ys.data());
        gScaled->SetMarkerColor(kRed);
        gScaled->SetLineColor(kRed);
        gScaled->SetMarkerStyle(kFullCircle);
        gScaled->SetMarkerSize(0.45);
        gScaled->Draw("P SAME");
        TGaxis* ax = new TGaxis(xrmax, y1min, xrmax, y1max, kDedxMin, kDedxMax, 510, "+L");
        ax->SetLineColor(kRed);
        ax->SetLabelColor(kRed);
        ax->SetTitleColor(kRed);
        ax->SetLabelSize(0.045);
        ax->SetTitleSize(0.048);
        ax->SetTitleOffset(1.15);
        ax->SetTitle("dE/dx [MeV/cm]");
        ax->Draw();
        TLegend* leg = new TLegend(0.50, 0.58, 0.88, 0.89);
        if (leg) {
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->AddEntry(b.mcsAngles, "Observed #Delta#theta vs RR (all MCS samples)", "p");
          if (b.mcsAnglesKept && b.mcsAnglesKept->GetN() > 0) {
            leg->AddEntry(b.mcsAnglesKept, "MCS samples used in fit", "p");
          }
          leg->AddEntry(b.mcsDedxVsRR, "dE/dx vs RR (collection)", "p");
          if (b.mcsNKept > 0 && b.mcsN > 0 && b.mcsNKept != b.mcsN) {
            leg->AddEntry((TObject*)nullptr, Form("N(all/kept)=%d/%d", b.mcsN, b.mcsNKept), "");
          } else {
            leg->AddEntry((TObject*)nullptr, Form("N=%d", b.mcsN), "");
          }
          if (std::isfinite(b.mcsMeanObs) && b.mcsMeanObs > -998.) {
            leg->AddEntry((TObject*)nullptr, Form("Mean(obs,kept)=%.4f rad", b.mcsMeanObs), "");
          } else {
            leg->AddEntry((TObject*)nullptr, "Mean(obs,kept)=n/a", "");
          }
          if (std::isfinite(b.mcsRmsObs) && b.mcsRmsObs > -998.) {
            leg->AddEntry((TObject*)nullptr, Form("RMS(obs,kept)=%.4f rad", b.mcsRmsObs), "");
          } else {
            leg->AddEntry((TObject*)nullptr, "RMS(obs,kept)=n/a", "");
          }
          if (b.mcsNKept > 0 && b.mcsN > 0 && b.mcsNKept != b.mcsN) {
            if (std::isfinite(b.mcsRmsObsAll) && b.mcsRmsObsAll > -998.) {
              leg->AddEntry((TObject*)nullptr, Form("RMS(obs,all)=%.4f rad", b.mcsRmsObsAll), "");
            } else {
              leg->AddEntry((TObject*)nullptr, "RMS(obs,all)=n/a", "");
            }
          }
          if (std::isfinite(b.mcsRmsTrue) && b.mcsRmsTrue > -998.) {
            leg->AddEntry((TObject*)nullptr, Form("RMS(true @ p_{true})=%.4f rad", b.mcsRmsTrue), "");
          } else {
            leg->AddEntry((TObject*)nullptr, "RMS(true @ p_{true})=n/a", "");
          }
          leg->Draw();
        }
        gPad->Modified();
        gPad->Update();
      }
    }
    cv->cd(3);
    if (b.tleLogL) b.tleLogL->Draw("A");
    cv->cd(4);
    if (b.mcsLogL) b.mcsLogL->Draw("A");
    file->cd();
    cv->Write(cv->GetName(), TObject::kOverwrite);
    delete cv;

    DeleteMultiGraphAndGraphs(b.dedx);
    DeleteMultiGraphAndGraphs(b.tleLogL);
    DeleteMultiGraphAndGraphs(b.mcsLogL);
    if (b.mcsAngles) delete b.mcsAngles;
    if (b.mcsAnglesKept) delete b.mcsAnglesKept;
    if (b.mcsDedxVsRR) delete b.mcsDedxVsRR;
  }
  sDiagByEvent.clear();
  sDiagSerialByCode.clear();
}

}  // namespace pionMomentumAnaUtils
