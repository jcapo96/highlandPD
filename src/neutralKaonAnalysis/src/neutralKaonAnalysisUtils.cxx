#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

namespace {

bool GetRecoDaughterTrueObjects(AnaNeutralParticlePD* neutralParticle,
                                AnaTrueParticlePD*& trueDaughter1,
                                AnaTrueParticlePD*& trueDaughter2) {
  trueDaughter1 = nullptr;
  trueDaughter2 = nullptr;

  if(!neutralParticle || !neutralParticle->AnnihilationVertex ||
     neutralParticle->AnnihilationVertex->Particles.size() < 2) {
    return false;
  }

  AnaParticlePD* daughter1 = neutralParticle->AnnihilationVertex->Particles[0];
  AnaParticlePD* daughter2 = neutralParticle->AnnihilationVertex->Particles[1];
  if(!daughter1 || !daughter2) {
    return false;
  }

  trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
  return trueDaughter1 && trueDaughter2;
}

}

//********************************************************************
void neutralKaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddSignalCandidateCategory();
  AddLooseSignalCandidateCategory();
}



//********************************************************************
void neutralKaonAnaUtils::AddSignalCandidateCategory(){

  std::string part_types[] = {
    "signal",                    // 1 - True K0 -> pi+ pi- with parent match
    "k0-no-decay",              // 2 - True K0 that do not decay
    "k0-decay-not-2pi",         // 3 - True K0 that decays but not into 2 pions
    "k0-decay-2pi-no-parent",   // 4 - True K0 that decays into two pions but parent doesn't match
    "background",               // 5 - Background (not a true K0 or other cases)
    NAMEOTHER};
  int part_codes[]         = {1, 2, 3, 4, 5, CATOTHER};
  int part_colors[]        = {2, 3, 4, 5, 6, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("signal", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::AddLooseSignalCandidateCategory(){

  std::string part_types[] = {
    "loose-signal",              // 1 - reco daughters share a true K0 parent decaying to pi+ pi-
    "true-k0-other-topology",    // 2 - common true parent is a K0 but topology differs
    "common-parent-not-k0",      // 3 - reco daughters share a true parent, but it is not a K0
    "no-common-true-parent",     // 4 - true daughters exist but do not share a parent
    "background",                // 5 - missing reco or true daughter information
    NAMEOTHER};
  int part_codes[]         = {1, 2, 3, 4, 5, CATOTHER};
  int part_colors[]        = {2, 3, 4, 5, 6, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("signal_loose", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  // Default to background
  int signalCode = 5; // background

  int run = -1;
  int subrun = -1;
  int evt = -1;
  if (event.EventInfo) {
    run = event.EventInfo->Run;
    subrun = event.EventInfo->SubRun;
    evt = event.EventInfo->Event;
  }

  auto printCategoryDecision = [&](const char* reason,
                                   int code,
                                   const AnaTrueParticlePD* trueNeutral,
                                   const AnaTrueParticlePD* trueDau1,
                                   const AnaTrueParticlePD* trueDau2,
                                   const AnaTrueParticlePD* parentTrue) {
    if (code != 1) return;
    std::cout << "[signal-category] run=" << run
              << " subrun=" << subrun
              << " event=" << evt
              << " neutralID=" << (neutralParticle ? neutralParticle->UniqueID : -1)
              << " code=" << code
              << " reason=" << reason
              << " truePDG=" << (trueNeutral ? trueNeutral->PDG : -999)
              << " parentRecoID=" << ((neutralParticle && neutralParticle->Parent) ? neutralParticle->Parent->UniqueID : -1)
              << " parentTrueID=" << (parentTrue ? parentTrue->ID : -1)
              << " dauPDG=(" << (trueDau1 ? trueDau1->PDG : -999)
              << "," << (trueDau2 ? trueDau2->PDG : -999) << ")"
              << std::endl;
  };

  // Check if neutral particle exists
  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("null-neutral-particle", signalCode, nullptr, nullptr, nullptr, nullptr);
    return;
  }

  // Get the true neutral particle
  AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  if(!trueNeutralParticle) {
    // No true object - background
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("no-true-object", signalCode, nullptr, nullptr, nullptr, nullptr);
    return;
  }

  // Check if the neutral particle is a K0 (PDG == 310)
  if(trueNeutralParticle->PDG != 310) {
    // Not a K0 - background
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("true-object-not-k0", signalCode, trueNeutralParticle, nullptr, nullptr, nullptr);
    return;
  }

  // We have a true K0 - now categorize based on decay properties
  // Check 1: Does the K0 decay? (ProcessEnd == 2 means decay)
  if(trueNeutralParticle->ProcessEnd != 2) {
    // True K0 that does not decay -> code 2
    signalCode = 2; // k0-no-decay
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("k0-no-decay", signalCode, trueNeutralParticle, nullptr, nullptr, nullptr);
    return;
  }

  // K0 decays - check decay products
  // Check 2: Does it decay into exactly 2 particles?
  if(trueNeutralParticle->Daughters.size() != 2) {
    // True K0 that decays but not into 2 particles -> code 3
    signalCode = 3; // k0-decay-not-2pi
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("k0-not-two-daughters", signalCode, trueNeutralParticle, nullptr, nullptr, nullptr);
    return;
  }

  // K0 decays into 2 particles - check if they are pions
  // Need to check annihilation vertex and daughter particles
  if(!neutralParticle->AnnihilationVertex ||
     neutralParticle->AnnihilationVertex->Particles.size() < 2) {
    // Can't check daughters - treat as decay but not 2 pions
    signalCode = 3; // k0-decay-not-2pi
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-annihilation-daughters", signalCode, trueNeutralParticle, nullptr, nullptr, nullptr);
    return;
  }

  // Get the two daughter particles from the annihilation vertex
  AnaParticlePD* daughter1 = neutralParticle->AnnihilationVertex->Particles[0];
  AnaParticlePD* daughter2 = neutralParticle->AnnihilationVertex->Particles[1];

  if(!daughter1 || !daughter2) {
    signalCode = 3; // k0-decay-not-2pi
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("null-reco-daughters", signalCode, trueNeutralParticle, nullptr, nullptr, nullptr);
    return;
  }

  // Get true objects of the daughters
  AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

  if(!trueDaughter1 || !trueDaughter2) {
    signalCode = 3; // k0-decay-not-2pi
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-true-daughters", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  // Check if the two daughters are pi+ and pi- (or vice versa)
  // PDG codes: 211 = pi+, -211 = pi-
  bool arePions = ((trueDaughter1->PDG == 211 && trueDaughter2->PDG == -211) ||
                    (trueDaughter1->PDG == -211 && trueDaughter2->PDG == 211));
  if(!arePions) {
    // True K0 that decays but not into 2 pions -> code 3
    signalCode = 3; // k0-decay-not-2pi
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("daughters-not-pi+-pi-", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  // K0 decays into 2 pions - now check parent match
  // Check if both daughters have the same true parent (same ParentID)
  if(trueDaughter1->ParentID <= 0 || trueDaughter2->ParentID <= 0) {
    // Can't check parent - treat as no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("daughter-parentid-invalid", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  if(trueDaughter1->ParentID != trueDaughter2->ParentID) {
    // Daughters have different parents - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("daughter-parents-differ", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  // Check if the reconstructed parent's true object matches the daughters' true parent
  if(!neutralParticle->Parent) {
    // No reconstructed parent - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-reco-parent", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  AnaTrueParticlePD* parentTrueObject = static_cast<AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
  if(!parentTrueObject) {
    // No true object for reconstructed parent - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-parent-trueobject", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, nullptr);
    return;
  }

  // Get the TrueParent of the true daughter and check if its ID matches parentTrueObject->ID
  AnaTrueParticlePD* trueDaughter1TrueParent = pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter1->ParentID);
  if(!trueDaughter1TrueParent) {
    // Could not get the true parent of the daughter - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-daughter1-true-parent", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, parentTrueObject);
    return;
  }

  AnaTrueParticlePD* trueDaughter2TrueParent = pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter2->ParentID);
  if(!trueDaughter2TrueParent) {
    // Could not get the true parent of the daughter - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-daughter2-true-parent", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, parentTrueObject);
    return;
  }

  AnaTrueParticlePD* trueGrandParent = pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter1TrueParent->ParentID);
  if(!trueGrandParent) {
    // Could not get the true grand parent - no parent match
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("missing-true-grandparent", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, parentTrueObject);
    return;
  }

  // Check if the TrueGrandParent's ID matches the parent's true object ID
  if(trueGrandParent->ID != parentTrueObject->ID) {
    // TrueGrandParent's ID doesn't match the parent's true object ID
    // True K0 that decays into two pions but parent doesn't match -> code 4
    signalCode = 4; // k0-decay-2pi-no-parent
    anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
    printCategoryDecision("grandparent-parent-mismatch", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, parentTrueObject);
    return;
  }

  signalCode = 1; // signal
  anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
  printCategoryDecision("signal", signalCode, trueNeutralParticle, trueDaughter1, trueDaughter2, parentTrueObject);
}

//********************************************************************
void neutralKaonAnaUtils::FillLooseSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  int signalCode = 5;

  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if(!GetRecoDaughterTrueObjects(neutralParticle, trueDaughter1, trueDaughter2)) {
    anaUtils::_categ->SetObjectCode("signal_loose", signalCode, CATOTHER, -1);
    return;
  }

  if(trueDaughter1->ParentID <= 0 || trueDaughter2->ParentID <= 0 ||
     trueDaughter1->ParentID != trueDaughter2->ParentID) {
    signalCode = 4;
    anaUtils::_categ->SetObjectCode("signal_loose", signalCode, CATOTHER, -1);
    return;
  }

  AnaTrueParticlePD* trueParent = pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter1->ParentID);
  if(!trueParent) {
    signalCode = 4;
    anaUtils::_categ->SetObjectCode("signal_loose", signalCode, CATOTHER, -1);
    return;
  }

  const bool arePions =
      ((trueDaughter1->PDG == 211 && trueDaughter2->PDG == -211) ||
       (trueDaughter1->PDG == -211 && trueDaughter2->PDG == 211));

  if(trueParent->PDG == 310) {
    signalCode = arePions ? 1 : 2;
  } else {
    signalCode = 3;
  }

  anaUtils::_categ->SetObjectCode("signal_loose", signalCode, CATOTHER, -1);
}

