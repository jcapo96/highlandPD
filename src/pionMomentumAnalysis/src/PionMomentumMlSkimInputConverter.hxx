#ifndef PionMomentumMlSkimInputConverter_h
#define PionMomentumMlSkimInputConverter_h

#include "InputConverter.hxx"
#include "pdDataClasses.hxx"

/// Reads flat `mlskim` TTree (pionMomentumSkim.exe output) into a minimal AnaSpillPD with one AnaTrueParticlePD
/// (the skim seed) per entry for downstream highland loops.
class PionMomentumMlSkimInputConverter : public InputConverter {
 public:
  PionMomentumMlSkimInputConverter();

  bool Initialize() override;
  Int_t GetSpill(Long64_t& entry, AnaSpillC*& spill) override;
  Int_t GetEvent(Long64_t& entry, AnaEventC*& event) override;
  bool AddFileToTChain(const std::string& inputString) override;
  void IncrementPOTBySpill() override;
  Int_t GetNEvents(Int_t entries = -1) override;
  void Reset() override;

 private:
  void SetBranchAddresses();

  Int_t _r_run;
  Int_t _r_subrun;
  Int_t _r_evt;
  Int_t _r_seed_true_id;
  Int_t _r_true_pdg;
  Float_t _r_true_start_momentum;
  Float_t _r_true_end_momentum;
  Float_t _r_true_start_x;
  Float_t _r_true_start_y;
  Float_t _r_true_start_z;
  Float_t _r_true_end_x;
  Float_t _r_true_end_y;
  Float_t _r_true_end_z;
  Int_t _r_true_parent_pdg;
  Int_t _r_true_gparent_pdg;
  Int_t _r_true_parent_id;
};

#endif
