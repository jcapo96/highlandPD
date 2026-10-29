#include "pionMomentumDump.hxx"
#include "OutputManager.hxx"
#include "pdUtilsDEdx.hxx"
#include "pdUtilsTrack.hxx"
#include <TFile.h>
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace pionMomentumDump {
namespace {

bool ValidXYZ(const TVector3& p) {
  return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900. && p.Y() > -900. &&
         p.Z() > -900.;
}

// True momentum [GeV/c] at the true trajectory point closest to a reco position (-999 if no trajectory).
float TrueMomentumNear(const AnaTrueParticlePD* t, const TVector3& recoPos) {
  if (!t || t->TrjPoints.empty() || !ValidXYZ(recoPos)) return -999.f;
  double best = 1e30;
  float p = -999.f;
  for (const auto& tp : t->TrjPoints) {
    const TVector3& q = ValidXYZ(tp.Position_NoSCE) ? tp.Position_NoSCE : tp.Position;
    if (!ValidXYZ(q)) continue;
    const double d = (q - recoPos).Mag2();
    if (d < best) {
      best = d;
      p = tp.Momentum;
    }
  }
  return p;
}

}  // namespace

void TrackDump::Book(OutputManager& output) {
  TTree* anaTree = output.GetTree();
  TFile* file = anaTree ? anaTree->GetCurrentFile() : nullptr;
  if (file) file->cd();
  std::cout << "[piondump] booking tree in file " << (file ? file->GetName() : "(none)") << std::endl;
  _tree = new TTree("piondump", "Per-track dump for pion momentum estimator studies");
  _tree->Branch("anaentry", &_anaentry);
  _tree->Branch("run", &_run);
  _tree->Branch("subrun", &_subrun);
  _tree->Branch("evt", &_evt);
  _tree->Branch("ismc", &_ismc);
  _tree->Branch("kind", &_kind);
  _tree->Branch("beammom", &_beammom);
  _tree->Branch("beammomtpc", &_beammomtpc);
  _tree->Branch("beampdgpi", &_beampdgpi);
  _tree->Branch("truepdg", &_truepdg);
  _tree->Branch("trueid", &_trueid);
  _tree->Branch("trueparentpdg", &_trueparentpdg);
  _tree->Branch("trueparentid", &_trueparentid);
  _tree->Branch("trueprocend", &_trueprocend);
  _tree->Branch("truep0", &_truep0);
  _tree->Branch("truepend", &_truepend);
  _tree->Branch("truepintpc", &_truepintpc);
  _tree->Branch("truepatstart", &_truepatstart);
  _tree->Branch("truelength", &_truelength);
  _tree->Branch("truestart", _truestart, "truestart[3]/F");
  _tree->Branch("trueend", _trueend, "trueend[3]/F");
  _tree->Branch("truedir", _truedir, "truedir[3]/F");
  _tree->Branch("length", &_length);
  _tree->Branch("start", _start, "start[3]/F");
  _tree->Branch("end", _end, "end[3]/F");
  _tree->Branch("dirstart", _dirstart, "dirstart[3]/F");
  _tree->Branch("ndau", &_ndau);
  _tree->Branch("ndesc", &_ndesc);
  _tree->Branch("ndeschits", &_ndeschits);
  _tree->Branch("descevis", &_descevis);
  _tree->Branch("isbeam", &_isbeam);
  _tree->Branch("chi2proton", &_chi2proton);
  _tree->Branch("chi2ndf", &_chi2ndf);
  _tree->Branch("hrr", &_hrr);
  _tree->Branch("hdedx", &_hdedx);
  _tree->Branch("hpitch", &_hpitch);
  _tree->Branch("hx", &_hx);
  _tree->Branch("hy", &_hy);
  _tree->Branch("hz", &_hz);
  _tree->Branch("hwire", &_hwire);
  _tree->Branch("htick", &_htick);
  _tree->Branch("htpc", &_htpc);
  _tree->Branch("tx", &_tx);
  _tree->Branch("ty", &_ty);
  _tree->Branch("tz", &_tz);
  _tree->Branch("ptle", &_ptle);
  _tree->Branch("pmcs", &_pmcs);
  _tree->Branch("pcsda", &_pcsda);
  _tree->Branch("prec", &_prec);
  _tree->Branch("precmeasured", &_precmeasured);
  _tree->Branch("tlep", &_tlep);
  _tree->Branch("tlell", &_tlell);
}

void TrackDump::FillTrack(AnaParticlePD* track, int kind, AnaBeamPD* beam, const DumpConfig& cfg) {
  _kind = kind;
  // PDSPAnalyzerTreeConverter stores the beam-instrumentation momentum (beam_inst_P) in BeamParticle->Momentum.
  _beammom = (beam && beam->BeamParticle) ? static_cast<Float_t>(beam->BeamParticle->Momentum) : -999.f;
  _beammomtpc = beam ? static_cast<Float_t>(beam->BeamMomentumInTPC) : -999.f;
  _beampdgpi = 0;
  if (beam)
    for (int pdg : beam->PDGs)
      if (std::abs(pdg) == 211) _beampdgpi = 1;

  const auto* t = track->TrueObject ? static_cast<const AnaTrueParticlePD*>(track->TrueObject) : nullptr;
  _truepdg = t ? t->PDG : -999;
  _trueid = t ? t->ID : -999;
  _trueparentpdg = t ? t->ParentPDG : -999;
  _trueparentid = t ? t->ParentID : -999;
  _trueprocend = t ? static_cast<Int_t>(t->ProcessEnd) : -999;
  _truep0 = t ? t->Momentum : -999.f;
  _truepend = t ? t->MomentumEnd : -999.f;
  _truepintpc = t ? t->MomentumInTPC : -999.f;
  _truelength = t ? t->Length : -999.f;
  for (int i = 0; i < 3; ++i) {
    _truestart[i] = t ? t->Position[i] : -999.f;
    _trueend[i] = t ? t->PositionEnd[i] : -999.f;
    _truedir[i] = t ? t->Direction[i] : -999.f;
    _start[i] = track->PositionStart[i];
    _end[i] = track->PositionEnd[i];
    _dirstart[i] = track->DirectionStart[i];
  }
  _truepatstart = TrueMomentumNear(t, TVector3(_start[0], _start[1], _start[2]));
  _length = track->Length;
  _ndau = static_cast<Int_t>(track->Daughters.size());
  _isbeam = track->isPandora ? 1 : 0;
  _chi2proton = track->Chi2Proton;
  _chi2ndf = track->Chi2ndf;

  std::vector<AnaParticlePD*> desc;
  pdMomReconstruction::CollectAllDescendants(track, desc);
  _ndesc = 0;
  _ndeschits = 0;
  _descevis = 0.f;
  for (AnaParticlePD* d : desc) {
    if (!d || d == track) continue;
    ++_ndesc;
    for (const AnaHitPD& h : d->Hits[2]) {
      if (!(h.dEdx > 0.f) || h.dEdx > 1000.f || !(h.Pitch > 0.f) || h.Pitch > 100.f) continue;
      _descevis += h.dEdx * h.Pitch;
      ++_ndeschits;
    }
  }

  _hrr.clear();
  _hdedx.clear();
  _hpitch.clear();
  _hx.clear();
  _hy.clear();
  _hz.clear();
  _hwire.clear();
  _htick.clear();
  _htpc.clear();
  if (track->Hits[2].size() >= 2 && track->Hits[2][0].ResidualRange <= 0.f) pdAnaUtils::ComputeResidualRange(track);
  for (const AnaHitPD& h : track->Hits[2]) {
    _hrr.push_back(h.ResidualRange);
    _hdedx.push_back(h.dEdx);
    _hpitch.push_back(h.Pitch);
    _hx.push_back(h.Position.X());
    _hy.push_back(h.Position.Y());
    _hz.push_back(h.Position.Z());
    _hwire.push_back(h.WireID.Wire);
    _htick.push_back(h.PeakTime);
    _htpc.push_back(h.WireID.TPC);
  }
  // Same point selection as BuildPionMcsScatteringSamples.
  _tx.clear();
  _ty.clear();
  _tz.clear();
  for (size_t i = 1; i < track->TrjPoints.size() && static_cast<int>(_tx.size()) < cfg.maxTrajPoints; ++i) {
    const AnaTrajectoryPointPD& tp = track->TrjPoints[i];
    if (!tp.IsInTPC || !ValidXYZ(tp.Position)) continue;
    _tx.push_back(tp.Position.X());
    _ty.push_back(tp.Position.Y());
    _tz.push_back(tp.Position.Z());
  }

  _ptle = -999.f;
  _prec = -999.f;
  _precmeasured = 0;
  _tlep.clear();
  _tlell.clear();
  _pcsda = pdMomReconstruction::EstimatePionMomentumFromCSDA(track);
  if (kind != 0 || cfg.tleOnBeamTrack) {
    pdMomReconstruction::PionMomentumRecoResult rec;
    if (pdMomReconstruction::EstimatePionMomentumRecommended(track, cfg.reco, rec)) {
      _prec = static_cast<Float_t>(rec.momentumGeV);
      _precmeasured = rec.measured ? 1 : 0;
      if (rec.tleValid) _ptle = static_cast<Float_t>(rec.tleMomentumGeV);
    }
    if (cfg.storeTleCurve) {
      pdMomReconstruction::PionTLEFitResult tleOut;
      std::vector<double> pAxis, logL;
      if (pdMomReconstruction::EstimatePionMomentumTLEFit(track, cfg.tle, tleOut, &pAxis, &logL) && tleOut.valid) {
        for (size_t i = 0; i < pAxis.size(); ++i) {
          if (pAxis[i] > 2.5) break;
          _tlep.push_back(static_cast<float>(pAxis[i]));
          _tlell.push_back(static_cast<float>(logL[i]));
        }
      }
    }
  }
  _pmcs = -999.f;
  pdMomReconstruction::PionMCSResult mcsOut;
  pdMomReconstruction::PionMCSRecoBuffers mcsBuf;
  if (pdMomReconstruction::EstimatePionMomentumMCSReco(track, cfg.mcs, mcsOut, mcsBuf) && mcsOut.valid)
    _pmcs = static_cast<Float_t>(mcsOut.bestMomentumGeV);

  _tree->Fill();
}

void TrackDump::FillEvent(OutputManager& output, AnaEventB& event, AnaParticlePD* mainTrack, AnaBeamPD* beam,
                          Long64_t anaEntry, const DumpConfig& cfg) {
  if (!_tree) Book(output);
  _anaentry = anaEntry;
  _run = event.EventInfo ? event.EventInfo->Run : -999;
  _subrun = event.EventInfo ? event.EventInfo->SubRun : -999;
  _evt = event.EventInfo ? event.EventInfo->Event : -999;
  _ismc = (event.EventInfo && event.GetIsMC()) ? 1 : 0;

  if (mainTrack) {
    FillTrack(mainTrack, 0, beam, cfg);
    if (_ismc) {
      for (AnaRecObjectC* obj : mainTrack->Daughters) {
        auto* dau = dynamic_cast<AnaParticlePD*>(obj);
        if (!dau || dau->Type != AnaParticlePD::kTrack || !dau->TrueObject) continue;
        const auto* t = static_cast<const AnaTrueParticlePD*>(dau->TrueObject);
        if (std::abs(t->PDG) != 211 || t->ParentPDG == 310) continue;  // K0S daughters are kind 2
        FillTrack(dau, 1, beam, cfg);
      }
    }
  }
  if (_ismc) {
    for (Int_t ip = 0; ip < event.nParticles; ++ip) {
      auto* part = static_cast<AnaParticlePD*>(event.Particles[ip]);
      if (!part || part == mainTrack || part->Type != AnaParticlePD::kTrack || !part->TrueObject) continue;
      const auto* t = static_cast<const AnaTrueParticlePD*>(part->TrueObject);
      if (std::abs(t->PDG) != 211 || t->ParentPDG != 310) continue;
      FillTrack(part, 2, beam, cfg);
    }
  }
}

void TrackDump::Write(OutputManager& output) {
  if (!_tree) return;
  TFile* file = _tree->GetCurrentFile();
  if (!file) {
    TTree* anaTree = output.GetTree();
    file = anaTree ? anaTree->GetCurrentFile() : nullptr;
  }
  std::cout << "[piondump] writing " << _tree->GetEntries() << " entries to " << (file ? file->GetName() : "(none)")
            << std::endl;
  if (!file) return;
  file->cd();
  _tree->Write("", TObject::kOverwrite);
  _tree->ResetBranchAddresses();
}

}  // namespace pionMomentumDump
