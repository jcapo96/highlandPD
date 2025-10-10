#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void neutralKaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddNeutralParticleSignalBackgroundCategory();
  AddNeutralParticlePDGCategory();
  AddNeutralParticleChargeCategory();
  AddNeutralParticleNoTruthCategory();
  AddSignalCandidateCategory();
}



//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticleSignalBackgroundCategory(){
  //********************************************************************

  std::string part_types[] = {
    // With truth information
    "k0-decay-pi+pi-",          // 1 - K0 decay -> pi+ pi-
    "k0-decay-gamma+gamma",     // 2 - K0 decay -> gamma gamma
    "k0-decay-other",           // 3 - K0 decay -> other particles
    "k0-decay-3+body",          // 4 - K0 decay -> 3+ particles
    "k0-no-decay-2pi",          // 5 - K0 no decay, vertex has 2 pions
    "k0-no-decay-1pi",          // 6 - K0 no decay, vertex has 1 pion
    "k0-no-decay-proton",       // 7 - K0 no decay, vertex has proton
    "k0-no-decay-other",        // 8 - K0 no decay, other vertex content
    "pi0-decay-gamma+gamma",    // 9 - pi0 decay -> gamma gamma
    "pi0-decay-other",          // 10 - pi0 decay -> other
    "other-parent",             // 11 - Other parent type
    // Without truth information - detailed parent analysis
    "same-parent-charged",     // 12 - No truth, same charged parent
    "same-parent-neutral",     // 13 - No truth, same neutral parent
    "diff-parent-charged",     // 14 - No truth, different charged parents
    "diff-parent-neutral",     // 15 - No truth, different neutral parents
    "mixed-parents",           // 16 - No truth, mixed charged/neutral parents
    "no-parent",               // 17 - No truth, no parent info
    NAMEOTHER};                         // 18 - Other cases
  int part_codes[]         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, CATOTHER};
  int part_colors[]        = {2, 3, 4, 5, 6, 7, 8, 9, 15, 25, 35, 45, 55, 65, 75, 85, 95, COLOTHER};
    const int NPART = sizeof(part_types)/sizeof(part_types[0]);

    std::reverse(part_types,  part_types  + NPART);
    std::reverse(part_codes,  part_codes  + NPART);
    std::reverse(part_colors, part_colors + NPART);

    // name to be called in terminal, counter to follow, variable name
    anaUtils::_categ->AddObjectCategory("isk0", neutralKaonTree::nk0, "nk0",
                NPART, part_types, part_codes, part_colors,
                1, -1000);
  }

//********************************************************************
void neutralKaonAnaUtils::AddSignalCandidateCategory(){

  std::string part_types[] = {"signal", "notruth", "k0-nodecay-pi+pi-", "nok0-decay-pi+pi-", "nok0-nodecay-pi+pi-", "1pion-1other", "2other", "nbody-decay", "nbody-other", NAMEOTHER};
  int part_codes[]         = {1, 2, 3, 4, 5, 6, 7, 8, 9, CATOTHER};
  int part_colors[]        = {2, 3, 4, 5, 6, 7, 8, 9, 31, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("signal", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1); // notruth
    return;
  }

  AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  if(!trueNeutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1); // notruth
    return;
  }
  else{
    AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(neutralParticle->Vertex->Particles[0]->TrueObject);
    AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(neutralParticle->Vertex->Particles[1]->TrueObject);
    AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);

    if(!trueDaughter1 || !trueDaughter2 || !trueParent) {
      anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1); // notruth
      return;
    }
    else{
      // If the candidate has 2 daughters
      if(trueNeutralParticle->Daughters.size() == 2) {
      // If the daughters are pions
      if((trueDaughter1->PDG==211 && trueDaughter2->PDG==-211) || (trueDaughter1->PDG==-211 && trueDaughter2->PDG==211)) {
        // If the candidate has 2 daughters
          if(trueNeutralParticle->PDG == 310) {
            // If the candidate decays
            if(trueNeutralParticle->ProcessEnd == 2) { // Decay
              std::cout << "SIGNAL CANDIDATE FOUND: K0 -> pi+ pi- decay" << std::endl;
              anaUtils::_categ->SetObjectCode("signal", 1, CATOTHER, -1); // k0
              return;
            }
            else{
              anaUtils::_categ->SetObjectCode("signal", 3, CATOTHER, -1); // k0 not decays but generates exactly 2 pions
              return;
            }
          }
          else{
            if(trueNeutralParticle->ProcessEnd == 2) {
              anaUtils::_categ->SetObjectCode("signal", 4, CATOTHER, -1); // no k0 but decays into exactly 2 pions
              return;
            }
            else{
              anaUtils::_categ->SetObjectCode("signal", 5, CATOTHER, -1); // no k0 but and does not decay, but generates exactly 2 pions
              return;
            }
          }
        }
        else if ((abs(trueDaughter1->PDG) == 211 && abs(trueDaughter2->PDG) != 211) || (abs(trueDaughter2->PDG) == 211 && abs(trueDaughter1->PDG) != 211)){
          anaUtils::_categ->SetObjectCode("signal", 6, CATOTHER, -1); // 1pion-1other
          return;
        }
        else{
          anaUtils::_categ->SetObjectCode("signal", 7, CATOTHER, -1); // 2other
          return;
        }
      }
      else{
        if(trueNeutralParticle->ProcessEnd == 2) {
          anaUtils::_categ->SetObjectCode("signal", 8, CATOTHER, -1); // nbody-decay
          return;
        }
        else{
          anaUtils::_categ->SetObjectCode("signal", 9, CATOTHER, -1); // nbody-other
          return;
        }
      }
    }
  }



}

//********************************************************************
void neutralKaonAnaUtils::FillNeutralParticleSignalBackgroundCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
  //********************************************************************

    if(!neutralParticle) {
      anaUtils::_categ->SetObjectCode("isk0", 17, CATOTHER, -1); // no-truth-no-parent
      return;
    }

    AnaTrueParticlePD* trueneutralparticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);

    // FIRST LEVEL: Check if we have truth information
    if(!trueneutralparticle) {
      // No truth information - check vertex particles
      if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
        anaUtils::_categ->SetObjectCode("isk0", 17, CATOTHER, -1); // no-truth-no-parent
        return;
      }

      AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
      AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

      if(!particle1 || !particle2) {
        anaUtils::_categ->SetObjectCode("isk0", 17, CATOTHER, -1); // no-truth-no-parent
        return;
      }

      // Check if particles have truth and same parent
      AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
      AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);

      if(!trueParticle1 || !trueParticle2) {
        anaUtils::_categ->SetObjectCode("isk0", 17, CATOTHER, -1); // no-truth-no-parent
        return;
      }

      // Helper function to determine if a particle is charged or neutral
      auto isChargedParticle = [](Int_t pdg) -> bool {
        Int_t abs_pdg = abs(pdg);
        // Charged particles: pions (211), kaons (321), protons (2212), muons (13), electrons (11), etc.
        // Neutral particles: neutrons (2112), photons (22), neutrinos, etc.
        // Special case: K0 (310) and pi0 (111) are neutral
        if (abs_pdg == 310 || abs_pdg == 111 || abs_pdg == 22 || abs_pdg == 2112) {
          return false; // neutral
        }
        // Most other particles with non-zero PDG are charged
        return (pdg != 0);
      };

      bool parent1_charged = isChargedParticle(trueParticle1->ParentPDG);
      bool parent2_charged = isChargedParticle(trueParticle2->ParentPDG);

      if(trueParticle1->ParentPDG == trueParticle2->ParentPDG) {
        // Same parent
        if (parent1_charged) {
          if(trueParticle1->ID == trueParticle2->ID) {
          anaUtils::_categ->SetObjectCode("isk0", 12, CATOTHER, -1); // same-parent-charged
          }
          else{
            anaUtils::_categ->SetObjectCode("isk0", 14, CATOTHER, -1); // no-truth-diff-parent-charged
          }
        } else {
          if(trueParticle1->ID == trueParticle2->ID) {
            anaUtils::_categ->SetObjectCode("isk0", 13, CATOTHER, -1); // same-parent-neutral
          }
          else{
            anaUtils::_categ->SetObjectCode("isk0", 15, CATOTHER, -1); // no-truth-diff-parent-neutral
          }
        }
      } else {
        // Different parents
        if (parent1_charged && parent2_charged) {
          anaUtils::_categ->SetObjectCode("isk0", 14, CATOTHER, -1); // no-truth-diff-parent-charged
        } else if (!parent1_charged && !parent2_charged) {
          anaUtils::_categ->SetObjectCode("isk0", 15, CATOTHER, -1); // no-truth-diff-parent-neutral
        } else {
          anaUtils::_categ->SetObjectCode("isk0", 16, CATOTHER, -1); // no-truth-mixed-parents
        }
      }
      return;
    }

    // SECOND LEVEL: We have truth information - check parent type
    if(trueneutralparticle->PDG == 310) { // K0
      // THIRD LEVEL: Check if K0 decays
      if(trueneutralparticle->ProcessEnd == AnaTrueParticleB::Decay) {
        // FOURTH LEVEL: Check decay products
        if(trueneutralparticle->Daughters.size() == 2) {
          // Find the true particles associated with the daughter IDs
          AnaTrueParticlePD* dau1 = nullptr;
          AnaTrueParticlePD* dau2 = nullptr;

          for(int i = 0; i < event.nTrueParticles; i++) {
            AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(event.TrueParticles[i]);
            if(truePart && truePart->ID == trueneutralparticle->Daughters[0]) {
              dau1 = truePart;
            }
            if(truePart && truePart->ID == trueneutralparticle->Daughters[1]) {
              dau2 = truePart;
            }
          }

          if(dau1 && dau2) {
            // Check for signal: pi+ and pi-
            if((dau1->PDG == 211 && dau2->PDG == -211) || (dau1->PDG == -211 && dau2->PDG == 211)) {
              anaUtils::_categ->SetObjectCode("isk0", 1, CATOTHER, -1); // k0-decay-pi+pi-
              return;
            }
            // Check for gamma gamma
            else if((dau1->PDG == 22 && dau2->PDG == 22)) {
              anaUtils::_categ->SetObjectCode("isk0", 2, CATOTHER, -1); // k0-decay-gamma+gamma
              return;
            }
            // Other particles
            else {
              anaUtils::_categ->SetObjectCode("isk0", 3, CATOTHER, -1); // k0-decay-other
              return;
            }
          } else {
            anaUtils::_categ->SetObjectCode("isk0", 3, CATOTHER, -1); // k0-decay-other
            return;
          }
        } else {
          anaUtils::_categ->SetObjectCode("isk0", 4, CATOTHER, -1); // k0-decay-3+body
          return;
        }
      } else {
        // K0 no decay - check vertex content
        if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
          anaUtils::_categ->SetObjectCode("isk0", 8, CATOTHER, -1); // k0-no-decay-other
          return;
        }

        AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
        AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

        if(!particle1 || !particle2) {
          anaUtils::_categ->SetObjectCode("isk0", 8, CATOTHER, -1); // k0-no-decay-other
          return;
        }

        AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
        AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);

        if(!trueParticle1 || !trueParticle2) {
          anaUtils::_categ->SetObjectCode("isk0", 8, CATOTHER, -1); // k0-no-decay-other
          return;
        }

        // Check vertex content
        bool hasPion1 = (trueParticle1->PDG == 211 || trueParticle1->PDG == -211);
        bool hasPion2 = (trueParticle2->PDG == 211 || trueParticle2->PDG == -211);
        bool hasProton1 = (trueParticle1->PDG == 2212);
        bool hasProton2 = (trueParticle2->PDG == 2212);

        if(hasPion1 && hasPion2) {
          anaUtils::_categ->SetObjectCode("isk0", 5, CATOTHER, -1); // k0-no-decay-2pi
        } else if(hasPion1 || hasPion2) {
          anaUtils::_categ->SetObjectCode("isk0", 6, CATOTHER, -1); // k0-no-decay-1pi
        } else if(hasProton1 || hasProton2) {
          anaUtils::_categ->SetObjectCode("isk0", 7, CATOTHER, -1); // k0-no-decay-proton
        } else {
          anaUtils::_categ->SetObjectCode("isk0", 8, CATOTHER, -1); // k0-no-decay-other
        }
        return;
      }
    } else if(trueneutralparticle->PDG == 111) { // pi0
      // Check vertex content for pi0
      if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
        anaUtils::_categ->SetObjectCode("isk0", 10, CATOTHER, -1); // pi0-decay-other
        return;
      }

      AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
      AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

      if(!particle1 || !particle2) {
        anaUtils::_categ->SetObjectCode("isk0", 10, CATOTHER, -1); // pi0-decay-other
        return;
      }

      AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
      AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);

      if(!trueParticle1 || !trueParticle2) {
        anaUtils::_categ->SetObjectCode("isk0", 10, CATOTHER, -1); // pi0-decay-other
        return;
      }

      // Check if vertex has 2 gammas
      if(trueParticle1->PDG == 22 && trueParticle2->PDG == 22) {
        anaUtils::_categ->SetObjectCode("isk0", 9, CATOTHER, -1); // pi0-decay-gamma+gamma
      } else {
        anaUtils::_categ->SetObjectCode("isk0", 10, CATOTHER, -1); // pi0-decay-other
      }
      return;
    } else {
      anaUtils::_categ->SetObjectCode("isk0", 11, CATOTHER, -1); // other-parent
      return;
    }
}



//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticlePDGCategory(){
//********************************************************************

  std::string part_types[] = {
    "notruth",           // 1 - No true object associated
    "k0",                // 2 - K0 (310)
    "pi0",               // 3 - pi0 (111)
    "gamma",             // 4 - gamma (22)
    "neutron",           // 5 - neutron (2112)
    "lambda",            // 6 - lambda (3122)
    "sigma0",            // 7 - sigma0 (3212)
    "k+",                // 8 - K+ (321)
    "k-",                // 9 - K- (-321)
    "pi+",               // 10 - pi+ (211)
    "pi-",               // 11 - pi- (-211)
    "e-",                // 12 - electron (11)
    "k0l",               // 13 - K0L (130)
    "proton",            // 15 - proton (2212)
    "eta",               // 16 - eta (221)
    "other",             // 17 - Other PDG codes
    NAMEOTHER};          // 18 - Other cases
  int part_codes[]         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, CATOTHER};
  int part_colors[]        = {1, 2, 3, 4, 5, 6, 7, 8, 9, 54, 21, 32, 43, 65, 76, 87, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  // name to be called in terminal, counter to follow, variable name
  anaUtils::_categ->AddObjectCategory("neutralpdg", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticleChargeCategory(){
//********************************************************************

  std::string part_types[] = {
    "notruth",           // 1 - No true object associated
    "charged",           // 2 - Charged particle
    "neutral",           // 3 - Neutral particle
    NAMEOTHER};          // 4 - Other cases
  int part_codes[]         = {1, 2, 3, CATOTHER};
  int part_colors[]        = {1, 2, 3, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  // name to be called in terminal, counter to follow, variable name
  anaUtils::_categ->AddObjectCategory("neutralcharge", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::FillNeutralParticlePDGCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 1, CATOTHER, -1); // notruth
    return;
  }

  // Check if we have true equivalent neutral particle
  AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  if(!trueNeutralParticle) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 1, CATOTHER, -1); // notruth
    return;
  }

  // Get the PDG code
  Int_t pdg = trueNeutralParticle->PDG;
  Int_t abs_pdg = abs(pdg);

  // Categorize based on PDG
  if(abs_pdg == 310) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 2, CATOTHER, -1); // k0/k0s
  } else if(abs_pdg == 111) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 3, CATOTHER, -1); // pi0
  } else if(abs_pdg == 22) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 4, CATOTHER, -1); // gamma
  } else if(abs_pdg == 2112) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 5, CATOTHER, -1); // neutron
  } else if(abs_pdg == 3122) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 6, CATOTHER, -1); // lambda
  } else if(abs_pdg == 3212) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 7, CATOTHER, -1); // sigma0
  } else if(abs_pdg == 321) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 8, CATOTHER, -1); // k+/k-
  } else if(abs_pdg == 211) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 10, CATOTHER, -1); // pi+/pi-
  } else if(abs_pdg == 11) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 12, CATOTHER, -1); // e-
  } else if(abs_pdg == 130) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 13, CATOTHER, -1); // k0l
  } else if(abs_pdg == 2212) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 15, CATOTHER, -1); // proton
  } else if(abs_pdg == 221) {
    anaUtils::_categ->SetObjectCode("neutralpdg", 16, CATOTHER, -1); // eta
  } else {
    anaUtils::_categ->SetObjectCode("neutralpdg", 17, CATOTHER, -1); // other
  }
}

//********************************************************************
void neutralKaonAnaUtils::FillNeutralParticleChargeCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("neutralcharge", 1, CATOTHER, -1); // notruth
    return;
  }

  // Check if we have true equivalent neutral particle
  AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  if(!trueNeutralParticle) {
    anaUtils::_categ->SetObjectCode("neutralcharge", 1, CATOTHER, -1); // notruth
    return;
  }

  // Get the PDG code
  Int_t pdg = trueNeutralParticle->PDG;
  Int_t abs_pdg = abs(pdg);

  // Helper function to determine if a particle is charged or neutral
  auto isChargedParticle = [](Int_t pdg) -> bool {
    Int_t abs_pdg = abs(pdg);
    // Charged particles: pions (211), kaons (321), protons (2212), muons (13), electrons (11), etc.
    // Neutral particles: neutrons (2112), photons (22), neutrinos, K0 (310), pi0 (111), etc.
    if (abs_pdg == 310 || abs_pdg == 111 || abs_pdg == 22 || abs_pdg == 2112 || abs_pdg==130 || abs_pdg==221) {
      return false; // neutral
    }
    // Most other particles with non-zero PDG are charged
    return (pdg != 0);
  };

  if(isChargedParticle(pdg)) {
    anaUtils::_categ->SetObjectCode("neutralcharge", 2, CATOTHER, -1); // charged
  } else {
    anaUtils::_categ->SetObjectCode("neutralcharge", 3, CATOTHER, -1); // neutral
  }
}

//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticleNoTruthCategory(){
//********************************************************************

  std::string part_types[] = {
    // No truth object cases - same true parent - charged parent
    "notruth-same-charged-reco-match-parent-is-daughter",      // 1
    "notruth-same-charged-reco-match-parent-not-daughter",     // 2
    "notruth-same-charged-no-match-parent-is-daughter",        // 3
    "notruth-same-charged-no-match-parent-not-daughter",       // 4
    // No truth object cases - same true parent - neutral parent
    "notruth-same-neutral-reco-match-parent-is-daughter",      // 5
    "notruth-same-neutral-reco-match-parent-not-daughter",     // 6
    "notruth-same-neutral-no-match-parent-is-daughter",        // 7
    "notruth-same-neutral-no-match-parent-not-daughter",       // 8
    // No truth object cases - different true parents
    "notruth-diff-parent-reco-match",                          // 9
    "notruth-diff-parent-no-match",                            // 10
    "notruth-no-parent-info",                                  // 11
    // Has truth object cases
    "hastruth-charged",                                        // 12
    "hastruth-neutral",                                        // 13
    NAMEOTHER};                                                // 14
  int part_codes[]         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, CATOTHER};
  int part_colors[]        = {2, 3, 4, 5, 6, 7, 8, 9, 25, 31, 42, 53, 64, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  // name to be called in terminal, counter to follow, variable name
  anaUtils::_categ->AddObjectCategory("notruth", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
void neutralKaonAnaUtils::FillNeutralParticleNoTruthCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("notruth", 11, CATOTHER, -1); // no-parent-info
    return;
  }

  AnaTrueParticlePD* trueNeutralParticle = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);

  // FIRST LEVEL: Check if neutralParticle has TrueObject
  if(!trueNeutralParticle) {
    // No truth object - analyze vertex particles
    if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
      anaUtils::_categ->SetObjectCode("notruth", 11, CATOTHER, -1); // no-parent-info
      return;
    }

    AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
    AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

    if(!particle1 || !particle2) {
      anaUtils::_categ->SetObjectCode("notruth", 11, CATOTHER, -1); // no-parent-info
      return;
    }

    // Check if particles have truth information
    AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
    AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);

    if(!trueParticle1 || !trueParticle2) {
      anaUtils::_categ->SetObjectCode("notruth", 11, CATOTHER, -1); // no-parent-info
      return;
    }

    // Helper function to determine if a particle is charged or neutral
    auto isChargedParticle = [](Int_t pdg) -> bool {
      Int_t abs_pdg = abs(pdg);
      // Neutral particles: K0 (310), pi0 (111), photon (22), neutron (2112), K0L (130), etc.
      if (abs_pdg == 310 || abs_pdg == 111 || abs_pdg == 22 || abs_pdg == 2112 ||
          abs_pdg == 130 || abs_pdg == 221 || abs_pdg == 3212) {
        return false; // neutral
      }
      // Most other particles with non-zero PDG are charged
      return (pdg != 0);
    };

    // Helper function to check if both particles are reco daughters of neutralParticle->Parent
    auto checkRecoMatch = [&]() -> bool {
      if(!neutralParticle->Parent) return false;

      AnaParticlePD* parent = neutralParticle->Parent;
      bool particle1IsDaughter = false;
      bool particle2IsDaughter = false;

      for(int i = 0; i < parent->Daughters.size(); i++) {
        if(parent->Daughters[i]) {
          if(parent->Daughters[i]->UniqueID == particle1->UniqueID) {
            particle1IsDaughter = true;
          }
          if(parent->Daughters[i]->UniqueID == particle2->UniqueID) {
            particle2IsDaughter = true;
          }
        }
      }

      return (particle1IsDaughter && particle2IsDaughter);
    };

    // Check if the two particles have the same true parent
    if(trueParticle1->ParentPDG == trueParticle2->ParentPDG && trueParticle1->ParentID == trueParticle2->ParentID) {
      // Same true parent - check if parent is charged or neutral AND check reco match
      bool parentCharged = isChargedParticle(trueParticle1->ParentPDG);
      bool recoMatch = checkRecoMatch();

      // Check if neutralParticle->Parent is a true daughter of the common true parent
      bool parentIsDaughter = false;
      if(neutralParticle->Parent) {
        AnaTrueParticlePD* trueRecoParent = static_cast<AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
        if(trueRecoParent) {
          // Check if the true object of neutralParticle->Parent has the same parent as the vertex particles
          if(trueRecoParent->ParentID == trueParticle1->ParentID &&
             trueRecoParent->ParentPDG == trueParticle1->ParentPDG) {
            parentIsDaughter = true;
          }
        }
      }

      if(parentCharged) {
        if(recoMatch) {
          if(parentIsDaughter) {
            anaUtils::_categ->SetObjectCode("notruth", 1, CATOTHER, -1); // same-charged-reco-match-parent-is-daughter
          } else {
            anaUtils::_categ->SetObjectCode("notruth", 2, CATOTHER, -1); // same-charged-reco-match-parent-not-daughter
          }
        } else {
          if(parentIsDaughter) {
            anaUtils::_categ->SetObjectCode("notruth", 3, CATOTHER, -1); // same-charged-no-match-parent-is-daughter
          } else {
            anaUtils::_categ->SetObjectCode("notruth", 4, CATOTHER, -1); // same-charged-no-match-parent-not-daughter
          }
        }
      } else {
        if(recoMatch) {
          if(parentIsDaughter) {
            anaUtils::_categ->SetObjectCode("notruth", 5, CATOTHER, -1); // same-neutral-reco-match-parent-is-daughter
          } else {
            anaUtils::_categ->SetObjectCode("notruth", 6, CATOTHER, -1); // same-neutral-reco-match-parent-not-daughter
          }
        } else {
          if(parentIsDaughter) {
            anaUtils::_categ->SetObjectCode("notruth", 7, CATOTHER, -1); // same-neutral-no-match-parent-is-daughter
          } else {
            anaUtils::_categ->SetObjectCode("notruth", 8, CATOTHER, -1); // same-neutral-no-match-parent-not-daughter
          }
        }
      }
      return;
    } else {
      // Different true parents - check if associated at reco level as daughters of neutralParticle->Parent
      bool recoMatch = checkRecoMatch();

      if(recoMatch) {
        anaUtils::_categ->SetObjectCode("notruth", 9, CATOTHER, -1); // diff-parent-reco-match
      } else {
        anaUtils::_categ->SetObjectCode("notruth", 10, CATOTHER, -1); // diff-parent-no-match
      }
      return;
    }
  } else {
    // Has truth object - check if it's charged or neutral
    Int_t pdg = trueNeutralParticle->PDG;
    Int_t abs_pdg = abs(pdg);

    // Helper function to determine if a particle is charged or neutral
    auto isChargedParticle = [](Int_t pdg) -> bool {
      Int_t abs_pdg = abs(pdg);
      // Neutral particles: K0 (310), pi0 (111), photon (22), neutron (2112), K0L (130), etc.
      if (abs_pdg == 310 || abs_pdg == 111 || abs_pdg == 22 || abs_pdg == 2112 ||
          abs_pdg == 130 || abs_pdg == 221 || abs_pdg == 3212) {
        return false; // neutral
      }
      // Most other particles with non-zero PDG are charged
      return (pdg != 0);
    };

    if(isChargedParticle(pdg)) {
      anaUtils::_categ->SetObjectCode("notruth", 12, CATOTHER, -1); // hastruth-charged
    } else {
      anaUtils::_categ->SetObjectCode("notruth", 13, CATOTHER, -1); // hastruth-neutral
    }
    return;
  }
}