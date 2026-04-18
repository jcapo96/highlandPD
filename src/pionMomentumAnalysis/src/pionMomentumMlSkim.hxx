#ifndef pionMomentumMlSkim_h
#define pionMomentumMlSkim_h

#include "pdDataClasses.hxx"

#include <string>

class TFile;
class TTree;

namespace pionMomentumMlSkimConstants {
  const Int_t kMaxRecoNodes = 128;
}

/// Writes flat TTree "mlskim" (truth seed rows + reco subgraph). Parameter keys: `<prefix>.*`.
class pionMomentumMlSkim {
 public:
  pionMomentumMlSkim();
  ~pionMomentumMlSkim();

  bool Initialize(const std::string& paramPrefix);
  void Finalize();

  void FillForBunch(AnaEventPD& event, AnaSpillPD& spill, const AnaBunchPD& bunch);

 private:
  void SetupTree();
  void ResetBranches();
  void FillOneRow(AnaEventPD& event, AnaSpillPD& spill, const AnaTrueParticlePD& trueSeed);
  bool TruePositionValid(const float pos[3]) const;
  bool PassesContainment(const AnaTrueParticlePD& t) const;
  static AnaParticlePD* FindRecoSeed(const AnaEventPD& event, const AnaTrueParticlePD& trueSeed, bool& ambiguous);
  static double SumSubtreeCollectionVisibleEnergyGeV(const AnaParticlePD* seed);

  std::string _prefix;
  TFile* _file;
  TTree* _tree;

  Int_t _truthSignalPdg;
  Double_t _truthMaxStartMomentumGeV;
  Int_t _applyTruthContainment;
  Double_t _cxMin, _cxMax, _cyMin, _cyMax, _czMin, _czMax;

  Int_t _br_run;
  Int_t _br_subrun;
  Int_t _br_evt;
  Int_t _br_seed_true_id;
  Int_t _br_true_pdg;
  Float_t _br_true_start_momentum;
  Float_t _br_true_end_momentum;
  Float_t _br_true_start_x;
  Float_t _br_true_start_y;
  Float_t _br_true_start_z;
  Float_t _br_true_end_x;
  Float_t _br_true_end_y;
  Float_t _br_true_end_z;
  Int_t _br_true_parent_pdg;
  Int_t _br_true_gparent_pdg;
  Int_t _br_true_parent_id;
  Int_t _br_reco_has_seed;
  Int_t _br_reco_ambiguous_match;
  Int_t _br_reco_seed_unique_id;
  Int_t _br_n_reco_nodes;
  Float_t _br_reco_subtree_visible_e_gev;
  Int_t _br_reco_unique_id[pionMomentumMlSkimConstants::kMaxRecoNodes];
  Float_t _br_reco_length[pionMomentumMlSkimConstants::kMaxRecoNodes];
  Float_t _br_reco_momentum[pionMomentumMlSkimConstants::kMaxRecoNodes];
  Int_t _br_reco_nhits_coll[pionMomentumMlSkimConstants::kMaxRecoNodes];
  Float_t _br_reco_visible_e_gev[pionMomentumMlSkimConstants::kMaxRecoNodes];
};

#endif
