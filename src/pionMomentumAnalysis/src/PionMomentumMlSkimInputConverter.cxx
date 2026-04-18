#include "PionMomentumMlSkimInputConverter.hxx"

#include <iostream>

//********************************************************************
PionMomentumMlSkimInputConverter::PionMomentumMlSkimInputConverter() : InputConverter("mlskim") {
//********************************************************************
  AddChain("mlskim");
  fChain = GetChain("mlskim");
  _nentries = 0;
}

//********************************************************************
void PionMomentumMlSkimInputConverter::Reset() {
//********************************************************************
  InputConverter::Reset();
}

//********************************************************************
bool PionMomentumMlSkimInputConverter::Initialize() {
//********************************************************************
  if (!fChain) return false;
  fChain->SetMakeClass(1);
  fCurrent = -1;
  SetBranchAddresses();
  return true;
}

//********************************************************************
void PionMomentumMlSkimInputConverter::SetBranchAddresses() {
//********************************************************************
  fChain->SetBranchAddress("run", &_r_run);
  fChain->SetBranchAddress("subrun", &_r_subrun);
  fChain->SetBranchAddress("event", &_r_evt);
  fChain->SetBranchAddress("seed_true_id", &_r_seed_true_id);
  fChain->SetBranchAddress("true_pdg", &_r_true_pdg);
  fChain->SetBranchAddress("true_start_momentum", &_r_true_start_momentum);
  fChain->SetBranchAddress("true_end_momentum", &_r_true_end_momentum);
  fChain->SetBranchAddress("true_start_x", &_r_true_start_x);
  fChain->SetBranchAddress("true_start_y", &_r_true_start_y);
  fChain->SetBranchAddress("true_start_z", &_r_true_start_z);
  fChain->SetBranchAddress("true_end_x", &_r_true_end_x);
  fChain->SetBranchAddress("true_end_y", &_r_true_end_y);
  fChain->SetBranchAddress("true_end_z", &_r_true_end_z);
  fChain->SetBranchAddress("true_parent_pdg", &_r_true_parent_pdg);
  fChain->SetBranchAddress("true_gparent_pdg", &_r_true_gparent_pdg);
  fChain->SetBranchAddress("true_parent_id", &_r_true_parent_id);
}

//********************************************************************
bool PionMomentumMlSkimInputConverter::AddFileToTChain(const std::string& inputString) {
//********************************************************************
  std::cout << "PionMomentumMlSkimInputConverter::AddFileToTChain " << inputString << std::endl;
  TChain probe("mlskim");
  probe.AddFile(inputString.c_str());
  if (probe.GetEntries() == 0) {
    std::cout << "      ----> No mlskim entries. IGNORED." << std::endl;
    return true;
  }
  return fChain->Add(inputString.c_str()) >= 0;
}

//********************************************************************
Int_t PionMomentumMlSkimInputConverter::GetNEvents(Int_t entries) {
//********************************************************************
  if (entries == -1) return static_cast<Int_t>(GetEntries());
  return entries;
}

//********************************************************************
void PionMomentumMlSkimInputConverter::IncrementPOTBySpill() {
//********************************************************************
}

//********************************************************************
Int_t PionMomentumMlSkimInputConverter::GetEvent(Long64_t& /*entry*/, AnaEventC*& /*event*/) {
//********************************************************************
  return 0;
}

//********************************************************************
Int_t PionMomentumMlSkimInputConverter::GetSpill(Long64_t& entry, AnaSpillC*& spill) {
//********************************************************************
  if (!fChain) return 0;

  const Int_t nbytes = fChain->GetEntry(entry);
  if (nbytes <= 0) return 0;

  auto* sp = new AnaSpillPD();
  spill = sp;

  auto* info = new AnaEventInfoPD();
  info->Run = _r_run;
  info->SubRun = _r_subrun;
  info->Event = _r_evt;
  info->IsMC = true;
  sp->EventInfo = info;

  auto* tp = new AnaTrueParticlePD();
  tp->ID = _r_seed_true_id;
  tp->PDG = _r_true_pdg;
  tp->ParentPDG = _r_true_parent_pdg;
  tp->GParentPDG = _r_true_gparent_pdg;
  tp->ParentID = _r_true_parent_id;
  tp->Momentum = _r_true_start_momentum;
  tp->MomentumEnd = _r_true_end_momentum;
  tp->Position[0] = _r_true_start_x;
  tp->Position[1] = _r_true_start_y;
  tp->Position[2] = _r_true_start_z;
  tp->Position[3] = 0.f;
  tp->PositionEnd[0] = _r_true_end_x;
  tp->PositionEnd[1] = _r_true_end_y;
  tp->PositionEnd[2] = _r_true_end_z;
  tp->PositionEnd[3] = 0.f;
  tp->Bunch = 0;

  sp->TrueParticles.push_back(tp);

  auto* bunch = new AnaBunchPD();
  bunch->Bunch = 0;
  sp->Bunches.push_back(bunch);

  sp->Beam = nullptr;

  sp->RedoLinks();

  entry++;
  return nbytes;
}
