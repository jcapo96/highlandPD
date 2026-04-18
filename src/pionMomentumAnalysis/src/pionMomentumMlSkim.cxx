#include "pionMomentumMlSkim.hxx"
#include "Parameters.hxx"
#include "pdMomReconstruction.hxx"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

inline bool TrueCoordsDefined(const float p[3]) {
  return p[0] > -900.f && p[1] > -900.f && p[2] > -900.f;
}

inline bool InBox(const float p[3], double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
  const double x = p[0];
  const double y = p[1];
  const double z = p[2];
  return x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin && z <= zmax;
}

}  // namespace

//********************************************************************
pionMomentumMlSkim::pionMomentumMlSkim()
    : _file(nullptr),
      _tree(nullptr),
      _truthSignalPdg(211),
      _truthMaxStartMomentumGeV(1.0),
      _applyTruthContainment(0),
      _cxMin(0),
      _cxMax(0),
      _cyMin(0),
      _cyMax(0),
      _czMin(0),
      _czMax(0) {
//********************************************************************
}

//********************************************************************
pionMomentumMlSkim::~pionMomentumMlSkim() {
//********************************************************************
  Finalize();
}

//********************************************************************
bool pionMomentumMlSkim::Initialize(const std::string& paramPrefix) {
//********************************************************************
  _prefix = paramPrefix;
  const std::string p = _prefix + ".";

  _truthSignalPdg = ND::params().GetParameterI(p + "TruthSignalPdg");
  _truthMaxStartMomentumGeV = ND::params().GetParameterD(p + "TruthMaxStartMomentumGeV");
  _applyTruthContainment = ND::params().GetParameterI(p + "ApplyTruthContainmentCut");
  _cxMin = ND::params().GetParameterD(p + "ContainmentXMin");
  _cxMax = ND::params().GetParameterD(p + "ContainmentXMax");
  _cyMin = ND::params().GetParameterD(p + "ContainmentYMin");
  _cyMax = ND::params().GetParameterD(p + "ContainmentYMax");
  _czMin = ND::params().GetParameterD(p + "ContainmentZMin");
  _czMax = ND::params().GetParameterD(p + "ContainmentZMax");

  const std::string outPath = ND::params().GetParameterS(p + "MlSkimOutputFile");
  if (outPath.empty()) {
    std::cerr << "pionMomentumMlSkim: " << p << "MlSkimOutputFile is empty" << std::endl;
    return false;
  }

  _file = TFile::Open(outPath.c_str(), "RECREATE");
  if (!_file || _file->IsZombie()) {
    std::cerr << "pionMomentumMlSkim: cannot create " << outPath << std::endl;
    return false;
  }
  SetupTree();
  return true;
}

//********************************************************************
void pionMomentumMlSkim::Finalize() {
//********************************************************************
  if (_file && _file->IsOpen()) {
    _file->cd();
    if (_tree) _tree->Write();
    _file->Write();
    _file->Close();
  }
  delete _file;
  _file = nullptr;
  _tree = nullptr;
}

//********************************************************************
void pionMomentumMlSkim::SetupTree() {
//********************************************************************
  _file->cd();
  _tree = new TTree("mlskim", "Truth pion rows + reco subgraph (one entry per truth seed)");
  _tree->Branch("run", &_br_run, "run/I");
  _tree->Branch("subrun", &_br_subrun, "subrun/I");
  _tree->Branch("event", &_br_evt, "event/I");
  _tree->Branch("seed_true_id", &_br_seed_true_id, "seed_true_id/I");
  _tree->Branch("true_pdg", &_br_true_pdg, "true_pdg/I");
  _tree->Branch("true_start_momentum", &_br_true_start_momentum, "true_start_momentum/F");
  _tree->Branch("true_end_momentum", &_br_true_end_momentum, "true_end_momentum/F");
  _tree->Branch("true_start_x", &_br_true_start_x, "true_start_x/F");
  _tree->Branch("true_start_y", &_br_true_start_y, "true_start_y/F");
  _tree->Branch("true_start_z", &_br_true_start_z, "true_start_z/F");
  _tree->Branch("true_end_x", &_br_true_end_x, "true_end_x/F");
  _tree->Branch("true_end_y", &_br_true_end_y, "true_end_y/F");
  _tree->Branch("true_end_z", &_br_true_end_z, "true_end_z/F");
  _tree->Branch("true_parent_pdg", &_br_true_parent_pdg, "true_parent_pdg/I");
  _tree->Branch("true_gparent_pdg", &_br_true_gparent_pdg, "true_gparent_pdg/I");
  _tree->Branch("true_parent_id", &_br_true_parent_id, "true_parent_id/I");
  _tree->Branch("reco_has_seed", &_br_reco_has_seed, "reco_has_seed/I");
  _tree->Branch("reco_ambiguous_match", &_br_reco_ambiguous_match, "reco_ambiguous_match/I");
  _tree->Branch("reco_seed_unique_id", &_br_reco_seed_unique_id, "reco_seed_unique_id/I");
  _tree->Branch("n_reco_nodes", &_br_n_reco_nodes, "n_reco_nodes/I");
  _tree->Branch("reco_subtree_visible_e_gev", &_br_reco_subtree_visible_e_gev, "reco_subtree_visible_e_gev/F");
  _tree->Branch("reco_unique_id", _br_reco_unique_id, "reco_unique_id[n_reco_nodes]/I");
  _tree->Branch("reco_length", _br_reco_length, "reco_length[n_reco_nodes]/F");
  _tree->Branch("reco_momentum", _br_reco_momentum, "reco_momentum[n_reco_nodes]/F");
  _tree->Branch("reco_nhits_coll", _br_reco_nhits_coll, "reco_nhits_coll[n_reco_nodes]/I");
  _tree->Branch("reco_visible_e_gev", _br_reco_visible_e_gev, "reco_visible_e_gev[n_reco_nodes]/F");
}

//********************************************************************
void pionMomentumMlSkim::ResetBranches() {
//********************************************************************
  _br_run = _br_subrun = _br_evt = -999;
  _br_seed_true_id = -999;
  _br_true_pdg = -999;
  _br_true_start_momentum = -999.f;
  _br_true_end_momentum = -999.f;
  _br_true_start_x = _br_true_start_y = _br_true_start_z = -999.f;
  _br_true_end_x = _br_true_end_y = _br_true_end_z = -999.f;
  _br_true_parent_pdg = _br_true_gparent_pdg = _br_true_parent_id = -999;
  _br_reco_has_seed = 0;
  _br_reco_ambiguous_match = 0;
  _br_reco_seed_unique_id = -999;
  _br_n_reco_nodes = 0;
  _br_reco_subtree_visible_e_gev = -999.f;
  const int nmax = pionMomentumMlSkimConstants::kMaxRecoNodes;
  for (int i = 0; i < nmax; ++i) {
    _br_reco_unique_id[i] = -999;
    _br_reco_length[i] = -999.f;
    _br_reco_momentum[i] = -999.f;
    _br_reco_nhits_coll[i] = -999;
    _br_reco_visible_e_gev[i] = -999.f;
  }
}

//********************************************************************
bool pionMomentumMlSkim::TruePositionValid(const float pos[3]) const {
//********************************************************************
  (void)this;
  return TrueCoordsDefined(pos);
}

//********************************************************************
bool pionMomentumMlSkim::PassesContainment(const AnaTrueParticlePD& t) const {
//********************************************************************
  if (!_applyTruthContainment) return true;
  if (!TruePositionValid(t.Position) || !TruePositionValid(t.PositionEnd)) return false;
  if (!InBox(t.Position, _cxMin, _cxMax, _cyMin, _cyMax, _czMin, _czMax)) return false;
  if (!InBox(t.PositionEnd, _cxMin, _cxMax, _cyMin, _cyMax, _czMin, _czMax)) return false;
  return true;
}

//********************************************************************
AnaParticlePD* pionMomentumMlSkim::FindRecoSeed(const AnaEventPD& event, const AnaTrueParticlePD& trueSeed,
                                                bool& ambiguous) {
//********************************************************************
  ambiguous = false;
  AnaParticlePD* found = nullptr;
  AnaParticleB** parts = event.Particles;
  const Int_t n = event.nParticles;
  for (Int_t i = 0; i < n; ++i) {
    auto* p = static_cast<AnaParticlePD*>(parts[i]);
    if (!p || !p->TrueObject) continue;
    auto* tp = static_cast<AnaTrueParticlePD*>(p->TrueObject);
    if (!tp || tp->ID != trueSeed.ID) continue;
    if (found) {
      ambiguous = true;
      return found;
    }
    found = p;
  }
  return found;
}

//********************************************************************
double pionMomentumMlSkim::SumSubtreeCollectionVisibleEnergyGeV(const AnaParticlePD* part) {
//********************************************************************
  if (!part) return 0.;
  std::vector<AnaParticlePD*> nodes;
  pdMomReconstruction::CollectAllDescendants(const_cast<AnaParticlePD*>(part), nodes);
  double sumMeV = 0.;
  for (AnaParticlePD* n : nodes) {
    const double e = pdMomReconstruction::CalculateDepositedEnergy(n, 2);
    if (e > 0. && e != -999.) sumMeV += e;
  }
  return sumMeV / 1000.;
}

//********************************************************************
void pionMomentumMlSkim::FillOneRow(AnaEventPD& event, AnaSpillPD& spill, const AnaTrueParticlePD& trueSeed) {
//********************************************************************
  ResetBranches();
  if (spill.EventInfo) {
    _br_run = spill.EventInfo->Run;
    _br_subrun = spill.EventInfo->SubRun;
    _br_evt = spill.EventInfo->Event;
  }
  _br_seed_true_id = trueSeed.ID;
  _br_true_pdg = trueSeed.PDG;
  _br_true_start_momentum = trueSeed.Momentum;
  _br_true_end_momentum = trueSeed.MomentumEnd;
  _br_true_start_x = trueSeed.Position[0];
  _br_true_start_y = trueSeed.Position[1];
  _br_true_start_z = trueSeed.Position[2];
  _br_true_end_x = trueSeed.PositionEnd[0];
  _br_true_end_y = trueSeed.PositionEnd[1];
  _br_true_end_z = trueSeed.PositionEnd[2];
  _br_true_parent_pdg = trueSeed.ParentPDG;
  _br_true_gparent_pdg = trueSeed.GParentPDG;
  _br_true_parent_id = trueSeed.ParentID;

  bool ambiguous = false;
  AnaParticlePD* seed = FindRecoSeed(event, trueSeed, ambiguous);
  _br_reco_ambiguous_match = ambiguous ? 1 : 0;
  if (seed) {
    _br_reco_has_seed = 1;
    _br_reco_seed_unique_id = seed->UniqueID;
    _br_reco_subtree_visible_e_gev = static_cast<Float_t>(SumSubtreeCollectionVisibleEnergyGeV(seed));

    std::vector<AnaParticlePD*> nodes;
    pdMomReconstruction::CollectAllDescendants(seed, nodes);
    const int nAll = static_cast<int>(nodes.size());
    _br_n_reco_nodes = std::min(nAll, pionMomentumMlSkimConstants::kMaxRecoNodes);
    for (int i = 0; i < _br_n_reco_nodes; ++i) {
      AnaParticlePD* d = nodes[static_cast<size_t>(i)];
      _br_reco_unique_id[i] = d->UniqueID;
      _br_reco_length[i] = d->Length;
      _br_reco_momentum[i] = d->Momentum;
      _br_reco_nhits_coll[i] = d->NHitsPerPlane[2];
      const double eOne = pdMomReconstruction::CalculateDepositedEnergy(d, 2);
      _br_reco_visible_e_gev[i] =
          (eOne > 0. && eOne != -999.) ? static_cast<Float_t>(eOne / 1000.) : -999.f;
    }
    if (nAll > pionMomentumMlSkimConstants::kMaxRecoNodes) {
      std::cerr << "pionMomentumMlSkim: truncated reco subgraph from " << nAll << " nodes to "
                << pionMomentumMlSkimConstants::kMaxRecoNodes << " (event " << _br_evt << " true id " << trueSeed.ID
                << ")" << std::endl;
    }
  }

  if (_tree) _tree->Fill();
}

//********************************************************************
void pionMomentumMlSkim::FillForBunch(AnaEventPD& event, AnaSpillPD& spill, const AnaBunchPD& bunch) {
//********************************************************************
  if (!_tree) return;

  const Int_t bunchIndex = bunch.Bunch;
  std::vector<AnaTrueParticleB*>& trueParts = spill.TrueParticles;
  for (AnaTrueParticleB* tp : trueParts) {
    if (!tp) continue;
    auto& t = *static_cast<AnaTrueParticlePD*>(tp);
    if (t.Bunch != bunchIndex) continue;
    if (t.PDG != _truthSignalPdg) continue;
    if (!std::isfinite(t.Momentum) || static_cast<double>(t.Momentum) > _truthMaxStartMomentumGeV) continue;
    if (!PassesContainment(t)) continue;

    FillOneRow(event, spill, t);
  }
}
