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
  AddNeutralParticleDetailedBackgroundCategory();
}

//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticleDetailedBackgroundCategory(){
//********************************************************************

  std::string part_types[] = {
    "signal",           // 1 - True K0 -> pi+ pi- decay
    "2-parent vtx",        // 2 - Different parents in the vertex
    "k0s-decay",        // 3 - Parent is K0S
    "pi0 parent",              // 4 - Parent is pi0
    "gamma parent",            // 5 - Parent is gamma
    "proton parent",          // 6 - Parent is charged particle
    "k0->pi0+pi0",      // 7 - K0 -> pi0 + pi0 decay
    "3-body decay",     // 8 - Parent has more than 2 daughters in the vertex
    "no truth",         // 9 - No truth information
    "k0s-non-decay",    // 10 - Parent is K0S but not a decay
    "p in vtx",         // 12 - proton in the vertex
    NAMEOTHER           // 11 - Other cases
  };

  int part_codes[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, CATOTHER};
  int part_colors[] = {2, 3, 4, 5, 6, 7, 8, 9, 1, 11, 18, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  // name to be called in terminal, counter to follow, variable name
  anaUtils::_categ->AddObjectCategory("isk0_detailed", neutralKaonTree::nk0, "nk0",
				      NPART, part_types, part_codes, part_colors,
				      1, -1000);
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
              std::cout << "SIGNAL CANDIDATE FOUND: K0 -> pi+ pi- decay" << std::endl;
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
void neutralKaonAnaUtils::FillNeutralParticleDetailedBackgroundCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) return;

  // Check if we have a vertex and particles
  if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", CATOTHER, CATOTHER, -1);
    return;
  }

  // Get the two particles from the vertex
  AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
  AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

  if(!particle1 || !particle2) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", CATOTHER, CATOTHER, -1);
  }

  // Get true particles
  AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
  AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);
  Int_t dau1PDG = trueParticle1->PDG;
  Int_t dau2PDG = trueParticle2->PDG;
  if(!trueParticle1 || !trueParticle2) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", 9, CATOTHER, -1); // Category no truth
  }

  // 1. Check if the two particles have the same parent
  if(trueParticle1->ParentID != trueParticle2->ParentID) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", 2, CATOTHER, -1); // Category diff. par
  }

  // 2. Find the parent particle
  AnaTrueParticlePD* parentTrueParticle = nullptr;
  for(int i = 0; i < event.nTrueParticles; i++) {
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(event.TrueParticles[i]);
    if(truePart && truePart->ID == trueParticle1->ParentID) {
      parentTrueParticle = truePart;
      break;
    }
  }

  if(!parentTrueParticle) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", 9, CATOTHER, -1); // Category no truth
    return;
  }

  // 3. Check if parent has exactly 2 daughters
  if(parentTrueParticle->Daughters.size() != 2) {
    anaUtils::_categ->SetObjectCode("isk0_detailed", 8, CATOTHER, -1); // Category 3-body decay
    return;
  }

  // 4. Check parent PDG
  Int_t parentPDG = parentTrueParticle->PDG;

  // 5. Check if parent is K0S (310)
  if(parentPDG == 310) {
    // Check if daughters are pi+ and pi-
    Int_t dau1PDG = trueParticle1->PDG;
    Int_t dau2PDG = trueParticle2->PDG;

    if((dau1PDG == 211 && dau2PDG == -211) || (dau1PDG == -211 && dau2PDG == 211)) {
      if(parentTrueParticle->ProcessEnd == AnaTrueParticleB::Decay) {
        anaUtils::_categ->SetObjectCode("isk0_detailed", 1, CATOTHER, -1); // Category signal
      } else {
        anaUtils::_categ->SetObjectCode("isk0_detailed", 10, CATOTHER, -1); // Category k0s-non-decay
      }
    } else if(dau1PDG == 111 && dau2PDG == 111) {
      anaUtils::_categ->SetObjectCode("isk0_detailed", 7, CATOTHER, -1); // Category k0->pi0+pi0
    } else {
      if(parentTrueParticle->ProcessEnd == AnaTrueParticleB::Decay) {
        anaUtils::_categ->SetObjectCode("isk0_detailed", 3, CATOTHER, -1);
        if(dau1PDG == 2212 || dau2PDG == 2212) {
          anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
        }
      } else {
        anaUtils::_categ->SetObjectCode("isk0_detailed", 10, CATOTHER, -1); // Category k0s-non-decay
        if(dau1PDG == 2212 || dau2PDG == 2212) {
          anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
        }
      }
    }
  }
  // 6. Check other parent types
  else if(parentPDG == 111) { // pi0
    anaUtils::_categ->SetObjectCode("isk0_detailed", 4, CATOTHER, -1); // Category pi0
    if(dau1PDG == 2212 || dau2PDG == 2212) {
      anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
    }
  } else if(parentPDG == 22) { // gamma
    anaUtils::_categ->SetObjectCode("isk0_detailed", 5, CATOTHER, -1); // Category gamma
    if(dau1PDG == 2212 || dau2PDG == 2212) {
      anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
    }
  } else if(parentPDG == 2212) { // proton
    anaUtils::_categ->SetObjectCode("isk0_detailed", 6, CATOTHER, -1); // Category proton
    if(dau1PDG == 2212 || dau2PDG == 2212) {
      anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
    }
  } else {
    anaUtils::_categ->SetObjectCode("isk0_detailed", CATOTHER, CATOTHER, -1); // Category other
    if(dau1PDG == 2212 || dau2PDG == 2212) {
      anaUtils::_categ->SetObjectCode("isk0_detailed", 12, CATOTHER, -1); // Category p in vtx
    }
  }
}