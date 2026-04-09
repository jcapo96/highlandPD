#include "pdUtilsNeutralKaonSignalDiagnostics.hxx"

#include "CategoryManager.hxx"
#include "Parameters.hxx"
#include "TFile.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "neutralKaonAnalysisUtils.hxx"
#include "pdUtilsDEdx.hxx"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_set>

namespace {
  Int_t sSignalDiagRun = -999999999;
  Int_t sSignalDiagSub = -999999999;
  Int_t sSignalDiagEvt = -999999999;
  std::unordered_set<Int_t> sSignalDiagTrueK0Ids;

  std::vector<TMultiGraph*> gK0SignalDedxMultiGraphs;
  std::vector<TH1F*> gK0SignalDedxBiasHistograms;
  /// Running index per signal subtype (two_stopping=1, one_stopping=5, interacting=6).
  int sK0SigDedxSerialByCode[3] = {0, 0, 0};

  int SignalCategoryToSerialSlot(int signalCode) {
    if (signalCode == 1) return 0;
    if (signalCode == 5) return 1;
    if (signalCode == 6) return 2;
    return -1;
  }

  void DeleteMultiGraphAndGraphs(TMultiGraph* mg) {
    if (!mg) return;
    TList* gl = mg->GetListOfGraphs();
    if (gl) gl->Delete();
    delete mg;
  }
}

namespace neutralKaonTreeDiagnostics {

bool RegisterSignalTrueK0Id(Int_t trueK0Id) {
  return sSignalDiagTrueK0Ids.insert(trueK0Id).second;
}

void MaybeAccumulateSignalPionDedxMultiGraphs(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                              size_t candidateIndex) {
  if (!candidate || !anaUtils::_categ || !anaUtils::_categ->HasCategory("signal")) return;
  const int signalCode = anaUtils::_categ->GetCategory("signal").GetObjectCode(1, static_cast<Int_t>(candidateIndex));
  if (signalCode != 1 && signalCode != 5 && signalCode != 6) return;

  const int codeSlot = SignalCategoryToSerialSlot(signalCode);
  if (codeSlot < 0) return;

  AnaAnnihilationVertexPD* vertex = candidate->AnnihilationVertex;
  if (!vertex || vertex->Particles.size() < 2) return;

  const int pairIdx = sK0SigDedxSerialByCode[codeSlot];
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
    const AnaTrueParticlePD* tpart = reco->TrueObject ? static_cast<AnaTrueParticlePD*>(reco->TrueObject) : nullptr;
    if (!tpart || std::abs(tpart->PDG) != 211) continue;

    double truncMinRR = 0.;
    double truncFrac = 0.;
    if (ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxLandauTruncMinRRCm")) {
      truncMinRR = ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxLandauTruncMinRRCm");
    }
    if (ND::params().HasParameter("neutralKaonAnalysis.FreeRangeDedxLandauTruncFraction")) {
      truncFrac = ND::params().GetParameterD("neutralKaonAnalysis.FreeRangeDedxLandauTruncFraction");
    }

    char dedxXAxis[200];
    std::snprintf(dedxXAxis, sizeof(dedxXAxis),
                  "(run %d, subrun %d, evt %d) Residual range [cm]", static_cast<int>(runId),
                  static_cast<int>(subRunId), static_cast<int>(evtId));
    TMultiGraph* mg = pdAnaUtils::MakePionFreeRangeDedxVsRRMultiGraph(reco, 500., 0.5, truncMinRR, truncFrac,
                                                                      dedxXAxis);
    if (!mg) continue;

    char biasHistTitle[240];
    std::snprintf(biasHistTitle, sizeof(biasHistTitle),
                  "(run %d, subrun %d, evt %d) #Delta(dE/dx)=measured-expected(PDF mode);#Delta(dE/dx) [MeV/cm];Entries",
                  static_cast<int>(runId), static_cast<int>(subRunId), static_cast<int>(evtId));
    TH1F* hBias = pdAnaUtils::MakePionFreeRangeDedxBiasHistogram(reco, 500., 0.5, truncMinRR, truncFrac,
                                                                 biasHistTitle);
    if (!hBias) {
      DeleteMultiGraphAndGraphs(mg);
      continue;
    }

    mg->SetName(Form("mg_sigdedx_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    hBias->SetName(Form("h_sigdedx_%d_dau%d_%d", signalCode, idau + 1, pairIdx));
    gK0SignalDedxMultiGraphs.push_back(mg);
    gK0SignalDedxBiasHistograms.push_back(hBias);
    anyGraph = true;
  }
  if (anyGraph) ++sK0SigDedxSerialByCode[codeSlot];
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
  for (TMultiGraph* mg : gK0SignalDedxMultiGraphs) {
    if (mg) {
      file->cd();
      mg->Write(mg->GetName(), TObject::kOverwrite);
      DeleteMultiGraphAndGraphs(mg);
    }
  }
  for (TH1F* h : gK0SignalDedxBiasHistograms) {
    if (h) {
      file->cd();
      h->Write(h->GetName(), TObject::kOverwrite);
      delete h;
    }
  }
  gK0SignalDedxMultiGraphs.clear();
  gK0SignalDedxBiasHistograms.clear();
  sK0SigDedxSerialByCode[0] = sK0SigDedxSerialByCode[1] = sK0SigDedxSerialByCode[2] = 0;
}

} // namespace neutralKaonTreeDiagnostics
