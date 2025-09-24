#include "neutralKaonTree.hxx"
#include "neutralKaonAnalysis.hxx"
#include "pdAnalysisUtils.hxx"
#include <cmath>
#include <set>
#include <iostream>

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_VertexCandidates(OutputManager& output, UInt_t nmax){
    //********************************************************************

  // Number of vertex candidates
  AddVarI(output, n_recovtx_candidates, "number of reconstructed vertex candidates");
  AddVarI(output, n_true_vertex_candidates, "number of true vertex candidates");

  // Reconstructed vertex candidates
  AddVarMaxSizeVI(output, recovtx_recoparticles, "reconstructed vertex number of particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_trueparticles, "true vertex number of particles", n_recovtx_candidates, nmax);
  AddVarMaxSize3MF(output, recovtx_recoposition, "reconstructed vertex position", n_recovtx_candidates, nmax);
  AddVarMaxSize3MF(output, recovtx_trueposition, "true vertex position", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_recoquality, "reconstructed vertex quality", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_truequality, "true vertex quality", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_par_truepdg, "reconstructed vertex parent true PDG", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_par_recopdg, "reconstructed vertex parent reconstructed PDG", n_recovtx_candidates, nmax);
  AddVarMaxSize3MF(output, recovtx_par_recodir, "reconstructed vertex parent reconstructed direction", n_recovtx_candidates, nmax);
  AddVarMaxSize3MF(output, recovtx_par_truedir, "reconstructed vertex parent true direction", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_recoangles, "reconstructed angles between parent end and daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_trueangles, "true angles between parent end and daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodistances, "reconstructed distances between daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_truedistances, "true distances between daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodistances_ss, "true distances between daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodistances_se, "true distances between daughter start/end positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodistances_es, "true distances between daughter end/start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodistances_ee, "true distances between daughter end/end positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_recodistances, "reconstructed distances between parent end and daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_truedistances, "true distances between parent end and daughter start positions", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_recosyst_angle, "reconstructed angle between parent direction and daughter system direction", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_par_dau_truesyst_angle, "true angle between parent direction and daughter system direction", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recolengths, "reconstructed lengths of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_truelengths, "true lengths of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_truepdgs, "true PDG codes of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_recopdgs, "reconstructed PDG codes of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_trueendproc, "true end processes of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_recoendproc, "reconstructed end processes of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_par_recoangles, "reconstructed angles between daughters and parent direction", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_par_trueangles, "true angles between daughters and parent direction", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_par_truepdg, "true PDG codes of parents of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVI(output, recovtx_dau_trueproc, "true process codes of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recomom, "reconstructed momenta of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_truemom, "true momenta of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_recodedx, "reconstructed dE/dx of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_dau_truededx, "true dE/dx of daughter particles", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_recomass, "reconstructed invariant mass of vertex parent particle", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_truemass, "true invariant mass of vertex parent particle", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_daupi_recodistance, "reconstructed distance between pion daughters from same K0", n_recovtx_candidates, nmax);
  AddVarMaxSizeVD(output, recovtx_daupi_truedistance, "true distance between pion daughters from same K0", n_recovtx_candidates, nmax);

  // True vertex candidates
  AddVarMaxSizeVI(output, true_vertex_nparticles, "true vertex number of particles", n_true_vertex_candidates, nmax);
  AddVarMaxSize3MF(output, true_vertex_position, "true vertex position", n_true_vertex_candidates, nmax);
  AddVarMaxSizeVD(output, true_vertex_quality, "true vertex quality", n_true_vertex_candidates, nmax);
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_VertexCandidates(OutputManager& output, const std::vector<AnaVertexPD*>& reconCandidates, const std::vector<AnaTrueVertexPD*>& trueCandidates){
  //********************************************************************

  // Fill reconstructed vertex candidates
  for (const auto& reconVertex : reconCandidates) {
    if (!reconVertex) continue;

    // Additional safety check: ensure all particles in vertex have unique IDs
    std::set<int> particleIDs;
    bool hasUniqueParticleIDs = true;
    for (int i = 0; i < reconVertex->NParticles; i++) {
      if (reconVertex->Particles[i]) {
        if (particleIDs.find(reconVertex->Particles[i]->UniqueID) != particleIDs.end()) {
          hasUniqueParticleIDs = false;
          break;
        }
        particleIDs.insert(reconVertex->Particles[i]->UniqueID);
      }
    }

    // Skip this vertex if particles have duplicate IDs
    if (!hasUniqueParticleIDs) {
      std::cout << "Warning: Skipping vertex with duplicate particle IDs" << std::endl;
      continue;
    }

    output.FillVectorVar(recovtx_recoparticles, (Int_t)reconVertex->NParticles);
    output.FillMatrixVarFromArray(recovtx_recoposition, reconVertex->Position, 3);
    output.FillVectorVar(recovtx_recoquality, (Double_t)reconVertex->NParticles); // Simple quality metric based on number of particles

    // Fill true information for reconstructed vertex (if available)
    Int_t trueParticles = -999;
    Float_t truePosition[3] = {-999.0, -999.0, -999.0};
    Double_t trueQuality = -999.0;

    // Try to get true vertex information from the first true particle if available
    if (reconVertex->NParticles > 0 && reconVertex->Particles[0]) {
      AnaTrueParticleB* trueParticle = reconVertex->Particles[0]->GetTrueParticle();
      if (trueParticle) {
        // Count true particles in this vertex
        trueParticles = reconVertex->NParticles; // For now, assume same count
        // Use reconstructed position as approximation for true position
        truePosition[0] = reconVertex->Position[0];
        truePosition[1] = reconVertex->Position[1];
        truePosition[2] = reconVertex->Position[2];
        trueQuality = (Double_t)trueParticles;
      }
    }

    output.FillVectorVar(recovtx_trueparticles, trueParticles);
    output.FillMatrixVarFromArray(recovtx_trueposition, truePosition, 3);
    output.FillVectorVar(recovtx_truequality, trueQuality);

    // Fill parent PDG (from true particle associated to reconstructed parent)
    Int_t parentTruePDG = -999;
    Int_t parentReconPDG = -999;
    if (reconVertex->Parent) {
      // Get the true particle associated with the reconstructed parent
      AnaTrueParticleB* trueParent = reconVertex->Parent->GetTrueParticle();
      if (trueParent) {
        parentTruePDG = trueParent->PDG; // Use the true PDG
      }
      // Always try to get reconstructed PDG
      parentReconPDG = reconVertex->Parent->ReconPDG[0];
    }
    output.FillVectorVar(recovtx_par_truepdg, parentTruePDG);
    output.FillVectorVar(recovtx_par_recopdg, parentReconPDG);

    // Fill particle direction (average direction of all particles in vertex)
    Float_t avgReconDirection[3] = {0.0, 0.0, 0.0};
    Float_t avgTrueDirection[3] = {-999.0, -999.0, -999.0};
    if (reconVertex->NParticles > 0) {
      for (const auto& particle : reconVertex->Particles) {
        if (particle) {
          avgReconDirection[0] += particle->DirectionStart[0];
          avgReconDirection[1] += particle->DirectionStart[1];
          avgReconDirection[2] += particle->DirectionStart[2];

          // Try to get true direction
          AnaTrueParticleB* trueParticle = particle->GetTrueParticle();
          if (trueParticle && avgTrueDirection[0] == -999.0) { // Only fill once
            avgTrueDirection[0] = trueParticle->Direction[0];
            avgTrueDirection[1] = trueParticle->Direction[1];
            avgTrueDirection[2] = trueParticle->Direction[2];
          }
        }
      }
      // Normalize reconstructed direction
      Float_t norm = sqrt(avgReconDirection[0]*avgReconDirection[0] + avgReconDirection[1]*avgReconDirection[1] + avgReconDirection[2]*avgReconDirection[2]);
      if (norm > 0) {
        avgReconDirection[0] /= norm;
        avgReconDirection[1] /= norm;
        avgReconDirection[2] /= norm;
      }
    }
    output.FillMatrixVarFromArray(recovtx_par_recodir, avgReconDirection, 3);
    output.FillMatrixVarFromArray(recovtx_par_truedir, avgTrueDirection, 3);

    // Calculate angles between parent end position and daughter start positions
    if (reconVertex->Parent && reconVertex->NParticles > 0) {
      for (const auto& daughter : reconVertex->Particles) {
        if (daughter) {
          // Use parent end position if available, otherwise use vertex position
          Float_t parentPos[3];
          if (reconVertex->Parent->PositionEnd[0] != -999 &&
              reconVertex->Parent->PositionEnd[1] != -999 &&
              reconVertex->Parent->PositionEnd[2] != -999) {
            parentPos[0] = reconVertex->Parent->PositionEnd[0];
            parentPos[1] = reconVertex->Parent->PositionEnd[1];
            parentPos[2] = reconVertex->Parent->PositionEnd[2];
          } else {
            // Fallback to vertex position
            parentPos[0] = reconVertex->Position[0];
            parentPos[1] = reconVertex->Position[1];
            parentPos[2] = reconVertex->Position[2];
          }

          // Calculate direction from parent end to daughter start
          Float_t direction[3] = {
            daughter->PositionStart[0] - parentPos[0],
            daughter->PositionStart[1] - parentPos[1],
            daughter->PositionStart[2] - parentPos[2]
          };

          // Normalize direction
          Float_t norm = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
          if (norm > 0) {
            direction[0] /= norm;
            direction[1] /= norm;
            direction[2] /= norm;
          }

          // Calculate angle with parent direction at end (if available)
          Float_t dotProduct = 0.0;
          if (reconVertex->Parent->DirectionEnd[0] != -999 &&
              reconVertex->Parent->DirectionEnd[1] != -999 &&
              reconVertex->Parent->DirectionEnd[2] != -999) {
            dotProduct = direction[0] * reconVertex->Parent->DirectionEnd[0] +
                        direction[1] * reconVertex->Parent->DirectionEnd[1] +
                        direction[2] * reconVertex->Parent->DirectionEnd[2];
          } else {
            // Fallback: use direction from parent start to end
            Float_t parentDir[3] = {
              reconVertex->Parent->PositionEnd[0] - reconVertex->Parent->PositionStart[0],
              reconVertex->Parent->PositionEnd[1] - reconVertex->Parent->PositionStart[1],
              reconVertex->Parent->PositionEnd[2] - reconVertex->Parent->PositionStart[2]
            };
            Float_t parentNorm = sqrt(parentDir[0]*parentDir[0] + parentDir[1]*parentDir[1] + parentDir[2]*parentDir[2]);
            if (parentNorm > 0) {
              parentDir[0] /= parentNorm;
              parentDir[1] /= parentNorm;
              parentDir[2] /= parentNorm;
              dotProduct = direction[0] * parentDir[0] + direction[1] * parentDir[1] + direction[2] * parentDir[2];
            }
          }

          // Clamp dot product to avoid numerical issues
          dotProduct = std::max(-1.0f, std::min(1.0f, dotProduct));
          Double_t angle = acos(dotProduct) * 180.0 / M_PI; // Convert to degrees
          output.FillVectorVar(recovtx_par_dau_recoangles, angle);

          // Calculate true angle if true particle is available
          Double_t trueAngle = -999.0;
          AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
          if (trueDaughter && reconVertex->Parent) {
            AnaTrueParticleB* trueParent = reconVertex->Parent->GetTrueParticle();
            if (trueParent) {
              // Use true positions and directions
              Float_t trueParentPos[3];
              if (trueParent->PositionEnd[0] != -999 &&
                  trueParent->PositionEnd[1] != -999 &&
                  trueParent->PositionEnd[2] != -999) {
                trueParentPos[0] = trueParent->PositionEnd[0];
                trueParentPos[1] = trueParent->PositionEnd[1];
                trueParentPos[2] = trueParent->PositionEnd[2];
              } else {
                trueParentPos[0] = trueParent->Position[0];
                trueParentPos[1] = trueParent->Position[1];
                trueParentPos[2] = trueParent->Position[2];
              }

              // Calculate direction from true parent end to true daughter start
              Float_t trueDirection[3] = {
                trueDaughter->Position[0] - trueParentPos[0],
                trueDaughter->Position[1] - trueParentPos[1],
                trueDaughter->Position[2] - trueParentPos[2]
              };

              // Normalize true direction
              Float_t trueNorm = sqrt(trueDirection[0]*trueDirection[0] + trueDirection[1]*trueDirection[1] + trueDirection[2]*trueDirection[2]);
              if (trueNorm > 0) {
                trueDirection[0] /= trueNorm;
                trueDirection[1] /= trueNorm;
                trueDirection[2] /= trueNorm;
              }

              // Calculate angle with true parent direction
              Float_t trueDotProduct = 0.0;
              if (trueParent->DirectionEnd[0] != -999 &&
                  trueParent->DirectionEnd[1] != -999 &&
                  trueParent->DirectionEnd[2] != -999) {
                trueDotProduct = trueDirection[0] * trueParent->DirectionEnd[0] +
                                trueDirection[1] * trueParent->DirectionEnd[1] +
                                trueDirection[2] * trueParent->DirectionEnd[2];
              } else {
                // Fallback: use direction from parent start to end
                Float_t trueParentDir[3] = {
                  trueParent->PositionEnd[0] - trueParent->Position[0],
                  trueParent->PositionEnd[1] - trueParent->Position[1],
                  trueParent->PositionEnd[2] - trueParent->Position[2]
                };
                Float_t trueParentNorm = sqrt(trueParentDir[0]*trueParentDir[0] + trueParentDir[1]*trueParentDir[1] + trueParentDir[2]*trueParentDir[2]);
                if (trueParentNorm > 0) {
                  trueParentDir[0] /= trueParentNorm;
                  trueParentDir[1] /= trueParentNorm;
                  trueParentDir[2] /= trueParentNorm;
                  trueDotProduct = trueDirection[0] * trueParentDir[0] + trueDirection[1] * trueParentDir[1] + trueDirection[2] * trueParentDir[2];
                }
              }

              // Clamp dot product to avoid numerical issues
              trueDotProduct = std::max(-1.0f, std::min(1.0f, trueDotProduct));
              trueAngle = acos(trueDotProduct) * 180.0 / M_PI; // Convert to degrees
            }
          }
          output.FillVectorVar(recovtx_par_dau_trueangles, trueAngle);
        }
      }
    }

    // Calculate distances between daughter start positions (excluding parent)
    if (reconVertex->NParticles > 1) {
      // Find daughters (exclude parent from distance calculations)
      std::vector<AnaParticlePD*> daughters;
      for (const auto& particle : reconVertex->Particles) {
        if (particle && particle != reconVertex->Parent) {
          daughters.push_back(particle);
        }
      }

      // Calculate distances between daughters only
      for (size_t i = 0; i < daughters.size(); i++) {
        for (size_t j = i + 1; j < daughters.size(); j++) {
          if (daughters[i] && daughters[j]) {
            // Reconstructed distances
            Float_t dx_ss = daughters[i]->PositionStart[0] - daughters[j]->PositionStart[0];
            Float_t dy_ss = daughters[i]->PositionStart[1] - daughters[j]->PositionStart[1];
            Float_t dz_ss = daughters[i]->PositionStart[2] - daughters[j]->PositionStart[2];
            Double_t distance_ss = sqrt(dx_ss*dx_ss + dy_ss*dy_ss + dz_ss*dz_ss);

            Float_t dx_se = daughters[i]->PositionStart[0] - daughters[j]->PositionEnd[0];
            Float_t dy_se = daughters[i]->PositionStart[1] - daughters[j]->PositionEnd[1];
            Float_t dz_se = daughters[i]->PositionStart[2] - daughters[j]->PositionEnd[2];
            Double_t distance_se = sqrt(dx_se*dx_se + dy_se*dy_se + dz_se*dz_se);

            Float_t dx_es = daughters[i]->PositionEnd[0] - daughters[j]->PositionStart[0];
            Float_t dy_es = daughters[i]->PositionEnd[1] - daughters[j]->PositionStart[1];
            Float_t dz_es = daughters[i]->PositionEnd[2] - daughters[j]->PositionStart[2];
            Double_t distance_es = sqrt(dx_es*dx_es + dy_es*dy_es + dz_es*dz_es);

            Float_t dx_ee = daughters[i]->PositionEnd[0] - daughters[j]->PositionEnd[0];
            Float_t dy_ee = daughters[i]->PositionEnd[1] - daughters[j]->PositionEnd[1];
            Float_t dz_ee = daughters[i]->PositionEnd[2] - daughters[j]->PositionEnd[2];
            Double_t distance_ee = sqrt(dx_ee*dx_ee + dy_ee*dy_ee + dz_ee*dz_ee);

            output.FillVectorVar(recovtx_dau_recodistances_ss, distance_ss);
            output.FillVectorVar(recovtx_dau_recodistances_se, distance_se);
            output.FillVectorVar(recovtx_dau_recodistances_es, distance_es);
            output.FillVectorVar(recovtx_dau_recodistances_ee, distance_ee);

            output.FillVectorVar(recovtx_dau_recodistances, std::min(distance_ss, std::min(distance_se, std::min(distance_es, distance_ee))));

            // True distances
            Double_t trueDistance = -999.0;
            AnaTrueParticleB* trueParticle1 = daughters[i]->GetTrueParticle();
            AnaTrueParticleB* trueParticle2 = daughters[j]->GetTrueParticle();
            if (trueParticle1 && trueParticle2) {
              Float_t trueDx = trueParticle1->Position[0] - trueParticle2->Position[0];
              Float_t trueDy = trueParticle1->Position[1] - trueParticle2->Position[1];
              Float_t trueDz = trueParticle1->Position[2] - trueParticle2->Position[2];
              trueDistance = sqrt(trueDx*trueDx + trueDy*trueDy + trueDz*trueDz);
            }
            output.FillVectorVar(recovtx_dau_truedistances, trueDistance);
          }
        }
      }

      // Calculate distances between pion daughters from the same K0
      std::vector<AnaParticlePD*> pionDaughters;

      // Find pion daughters from K0 parents
      for (const auto& daughter : daughters) {
        if (daughter) {
          AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
          if (trueDaughter) {
            // Check if this is a pion (PDG codes: 211, -211, 111)
            if (abs(trueDaughter->PDG) == 211 || trueDaughter->PDG == 111) {
              // Check if this pion has a K0 parent (PDG 310)
              if (trueDaughter->ParentPDG == 310) {
                pionDaughters.push_back(daughter);
              }
            }
          }
        }
      }

      // Calculate distances between pion daughters from the same K0
      Double_t pionRecoDistance = -999.0;
      Double_t pionTrueDistance = -999.0;

      // Check if we have at least 2 pion daughters and they have the same K0 parent
      if (pionDaughters.size() >= 2) {
        AnaTrueParticleB* truePion1 = pionDaughters[0]->GetTrueParticle();
        AnaTrueParticleB* truePion2 = pionDaughters[1]->GetTrueParticle();

        // Verify both pions have the same K0 parent
        if (truePion1 && truePion2 &&
            truePion1->ParentID == truePion2->ParentID &&
            truePion1->ParentID != 0) {

          // Calculate reconstructed distance between pion daughters
          Float_t dx = pionDaughters[0]->PositionStart[0] - pionDaughters[1]->PositionStart[0];
          Float_t dy = pionDaughters[0]->PositionStart[1] - pionDaughters[1]->PositionStart[1];
          Float_t dz = pionDaughters[0]->PositionStart[2] - pionDaughters[1]->PositionStart[2];
          pionRecoDistance = sqrt(dx*dx + dy*dy + dz*dz);

          // Calculate true distance between pion daughters
          Float_t trueDx = truePion1->Position[0] - truePion2->Position[0];
          Float_t trueDy = truePion1->Position[1] - truePion2->Position[1];
          Float_t trueDz = truePion1->Position[2] - truePion2->Position[2];
          pionTrueDistance = sqrt(trueDx*trueDx + trueDy*trueDy + trueDz*trueDz);
        }
      }

      output.FillVectorVar(recovtx_daupi_recodistance, pionRecoDistance);
      output.FillVectorVar(recovtx_daupi_truedistance, pionTrueDistance);
    }

    // Calculate distances between parent end position and daughter start positions
    if (reconVertex->Parent && reconVertex->NParticles > 0) {
      for (const auto& daughter : reconVertex->Particles) {
        if (daughter) {
          // Use parent end position if available, otherwise use vertex position
          Float_t parentPos[3];
          if (reconVertex->Parent->PositionEnd[0] != -999 &&
              reconVertex->Parent->PositionEnd[1] != -999 &&
              reconVertex->Parent->PositionEnd[2] != -999) {
            parentPos[0] = reconVertex->Parent->PositionEnd[0];
            parentPos[1] = reconVertex->Parent->PositionEnd[1];
            parentPos[2] = reconVertex->Parent->PositionEnd[2];
          } else {
            // Fallback to vertex position
            parentPos[0] = reconVertex->Position[0];
            parentPos[1] = reconVertex->Position[1];
            parentPos[2] = reconVertex->Position[2];
          }

          Float_t dx = daughter->PositionStart[0] - parentPos[0];
          Float_t dy = daughter->PositionStart[1] - parentPos[1];
          Float_t dz = daughter->PositionStart[2] - parentPos[2];
          Double_t distance = sqrt(dx*dx + dy*dy + dz*dz);
          output.FillVectorVar(recovtx_par_dau_recodistances, distance);

          // True distances between parent and daughter
          Double_t trueDistance = -999.0;
          AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
          if (trueDaughter && reconVertex->Parent) {
            AnaTrueParticleB* trueParent = reconVertex->Parent->GetTrueParticle();
            if (trueParent) {
              Float_t trueParentPos[3];
              if (trueParent->PositionEnd[0] != -999 &&
                  trueParent->PositionEnd[1] != -999 &&
                  trueParent->PositionEnd[2] != -999) {
                trueParentPos[0] = trueParent->PositionEnd[0];
                trueParentPos[1] = trueParent->PositionEnd[1];
                trueParentPos[2] = trueParent->PositionEnd[2];
              } else {
                trueParentPos[0] = trueParent->Position[0];
                trueParentPos[1] = trueParent->Position[1];
                trueParentPos[2] = trueParent->Position[2];
              }

              Float_t trueDx = trueDaughter->Position[0] - trueParentPos[0];
              Float_t trueDy = trueDaughter->Position[1] - trueParentPos[1];
              Float_t trueDz = trueDaughter->Position[2] - trueParentPos[2];
              trueDistance = sqrt(trueDx*trueDx + trueDy*trueDy + trueDz*trueDz);
            }
          }
          output.FillVectorVar(recovtx_par_dau_truedistances, trueDistance);
        }
      }
    }

    // Calculate angle between parent direction and daughter system direction
    if (reconVertex->Parent && reconVertex->NParticles > 0) {
      // Calculate daughter system direction (resultant direction of all daughters)
      Float_t daughterSystemDirection[3] = {0.0, 0.0, 0.0};
      for (const auto& daughter : reconVertex->Particles) {
        if (daughter) {
          daughterSystemDirection[0] += daughter->DirectionStart[0];
          daughterSystemDirection[1] += daughter->DirectionStart[1];
          daughterSystemDirection[2] += daughter->DirectionStart[2];
        }
      }

      // Normalize daughter system direction
      Float_t daughterNorm = sqrt(daughterSystemDirection[0]*daughterSystemDirection[0] +
                                 daughterSystemDirection[1]*daughterSystemDirection[1] +
                                 daughterSystemDirection[2]*daughterSystemDirection[2]);
      if (daughterNorm > 0) {
        daughterSystemDirection[0] /= daughterNorm;
        daughterSystemDirection[1] /= daughterNorm;
        daughterSystemDirection[2] /= daughterNorm;
      }

      // Get parent direction at end (if available)
      Float_t parentDirection[3];
      if (reconVertex->Parent->DirectionEnd[0] != -999 &&
          reconVertex->Parent->DirectionEnd[1] != -999 &&
          reconVertex->Parent->DirectionEnd[2] != -999) {
        parentDirection[0] = reconVertex->Parent->DirectionEnd[0];
        parentDirection[1] = reconVertex->Parent->DirectionEnd[1];
        parentDirection[2] = reconVertex->Parent->DirectionEnd[2];
      } else {
        // Fallback: calculate direction from parent start to end
        parentDirection[0] = reconVertex->Parent->PositionEnd[0] - reconVertex->Parent->PositionStart[0];
        parentDirection[1] = reconVertex->Parent->PositionEnd[1] - reconVertex->Parent->PositionStart[1];
        parentDirection[2] = reconVertex->Parent->PositionEnd[2] - reconVertex->Parent->PositionStart[2];

        // Normalize parent direction
        Float_t parentNorm = sqrt(parentDirection[0]*parentDirection[0] +
                                 parentDirection[1]*parentDirection[1] +
                                 parentDirection[2]*parentDirection[2]);
        if (parentNorm > 0) {
          parentDirection[0] /= parentNorm;
          parentDirection[1] /= parentNorm;
          parentDirection[2] /= parentNorm;
        }
      }

      // Calculate angle between parent direction and daughter system direction
      Float_t dotProduct = parentDirection[0] * daughterSystemDirection[0] +
                          parentDirection[1] * daughterSystemDirection[1] +
                          parentDirection[2] * daughterSystemDirection[2];

      // Clamp dot product to avoid numerical issues
      dotProduct = std::max(-1.0f, std::min(1.0f, dotProduct));
      Double_t angle = acos(dotProduct) * 180.0 / M_PI; // Convert to degrees
      output.FillVectorVar(recovtx_par_dau_recosyst_angle, angle);

      // Calculate true system angle
      Double_t trueSystemAngle = -999.0;
      if (reconVertex->Parent) {
        AnaTrueParticleB* trueParent = reconVertex->Parent->GetTrueParticle();
        if (trueParent) {
          // Calculate true daughter system direction
          Float_t trueDaughterSystemDirection[3] = {0.0, 0.0, 0.0};
          for (const auto& daughter : reconVertex->Particles) {
            if (daughter) {
              AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
              if (trueDaughter) {
                trueDaughterSystemDirection[0] += trueDaughter->Direction[0];
                trueDaughterSystemDirection[1] += trueDaughter->Direction[1];
                trueDaughterSystemDirection[2] += trueDaughter->Direction[2];
              }
            }
          }

          // Normalize true daughter system direction
          Float_t trueDaughterNorm = sqrt(trueDaughterSystemDirection[0]*trueDaughterSystemDirection[0] +
                                         trueDaughterSystemDirection[1]*trueDaughterSystemDirection[1] +
                                         trueDaughterSystemDirection[2]*trueDaughterSystemDirection[2]);
          if (trueDaughterNorm > 0) {
            trueDaughterSystemDirection[0] /= trueDaughterNorm;
            trueDaughterSystemDirection[1] /= trueDaughterNorm;
            trueDaughterSystemDirection[2] /= trueDaughterNorm;
          }

          // Get true parent direction at end
          Float_t trueParentDirection[3];
          if (trueParent->DirectionEnd[0] != -999 &&
              trueParent->DirectionEnd[1] != -999 &&
              trueParent->DirectionEnd[2] != -999) {
            trueParentDirection[0] = trueParent->DirectionEnd[0];
            trueParentDirection[1] = trueParent->DirectionEnd[1];
            trueParentDirection[2] = trueParent->DirectionEnd[2];
          } else {
            // Fallback: calculate direction from parent start to end
            trueParentDirection[0] = trueParent->PositionEnd[0] - trueParent->Position[0];
            trueParentDirection[1] = trueParent->PositionEnd[1] - trueParent->Position[1];
            trueParentDirection[2] = trueParent->PositionEnd[2] - trueParent->Position[2];

            // Normalize true parent direction
            Float_t trueParentNorm = sqrt(trueParentDirection[0]*trueParentDirection[0] +
                                         trueParentDirection[1]*trueParentDirection[1] +
                                         trueParentDirection[2]*trueParentDirection[2]);
            if (trueParentNorm > 0) {
              trueParentDirection[0] /= trueParentNorm;
              trueParentDirection[1] /= trueParentNorm;
              trueParentDirection[2] /= trueParentNorm;
            }
          }

          // Calculate angle between true parent direction and true daughter system direction
          Float_t trueDotProduct = trueParentDirection[0] * trueDaughterSystemDirection[0] +
                                  trueParentDirection[1] * trueDaughterSystemDirection[1] +
                                  trueParentDirection[2] * trueDaughterSystemDirection[2];

          // Clamp dot product to avoid numerical issues
          trueDotProduct = std::max(-1.0f, std::min(1.0f, trueDotProduct));
          trueSystemAngle = acos(trueDotProduct) * 180.0 / M_PI; // Convert to degrees
        }
      }
      output.FillVectorVar(recovtx_par_dau_truesyst_angle, trueSystemAngle);
    }

    // Fill daughter particle information (all variables in a single loop to avoid overwriting)
    if (reconVertex->NParticles > 0) {
      for (const auto& daughter : reconVertex->Particles) {
        if (daughter) {
          // Calculate reconstructed length
          Double_t length = -999.0;
          if (daughter->PositionStart[0] != -999 && daughter->PositionStart[1] != -999 && daughter->PositionStart[2] != -999 &&
              daughter->PositionEnd[0] != -999 && daughter->PositionEnd[1] != -999 && daughter->PositionEnd[2] != -999) {
            // Calculate 3D distance between start and end positions
            Float_t dx = daughter->PositionEnd[0] - daughter->PositionStart[0];
            Float_t dy = daughter->PositionEnd[1] - daughter->PositionStart[1];
            Float_t dz = daughter->PositionEnd[2] - daughter->PositionStart[2];
            length = sqrt(dx*dx + dy*dy + dz*dz);
          }
          output.FillVectorVar(recovtx_dau_recolengths, length);

          // Calculate true length
          Double_t trueLength = -999.0;
          AnaTrueParticleB* trueDaughter = daughter->GetTrueParticle();
          if (trueDaughter) {
            if (trueDaughter->Position[0] != -999 && trueDaughter->Position[1] != -999 && trueDaughter->Position[2] != -999 &&
                trueDaughter->PositionEnd[0] != -999 && trueDaughter->PositionEnd[1] != -999 && trueDaughter->PositionEnd[2] != -999) {
              Float_t trueDx = trueDaughter->PositionEnd[0] - trueDaughter->Position[0];
              Float_t trueDy = trueDaughter->PositionEnd[1] - trueDaughter->Position[1];
              Float_t trueDz = trueDaughter->PositionEnd[2] - trueDaughter->Position[2];
              trueLength = sqrt(trueDx*trueDx + trueDy*trueDy + trueDz*trueDz);
            }
          }
          output.FillVectorVar(recovtx_dau_truelengths, trueLength);

          // Fill daughter PDG (from true particle if available)
          Int_t daughterTruePDG = -999;
          Int_t daughterReconPDG = -999;
          if (trueDaughter) {
            daughterTruePDG = trueDaughter->PDG;
          }
          // Always try to get reconstructed PDG
          daughterReconPDG = daughter->ReconPDG[0];
          output.FillVectorVar(recovtx_dau_truepdgs, daughterTruePDG);
          output.FillVectorVar(recovtx_dau_recopdgs, daughterReconPDG);

          // Fill parent true PDG of daughter particle
          if (trueDaughter) {
            output.FillVectorVar(recovtx_dau_par_truepdg, trueDaughter->ParentPDG);
          }

          // Fill true process code of daughter particle
          Int_t daughterTrueProcess = -999;
          if (trueDaughter) {
            daughterTrueProcess = trueDaughter->ProcessStart;
          }
          output.FillVectorVar(recovtx_dau_trueproc, daughterTrueProcess);

          // Fill daughter end process (from true particle if available)
          Int_t trueEndProcess = -999;
          Int_t reconEndProcess = -999;
          if (trueDaughter) {
            trueEndProcess = trueDaughter->ProcessEnd;
          }
          // For reconstructed end process, we don't have this information directly
          // so we'll set it to -999 for now
          output.FillVectorVar(recovtx_dau_trueendproc, trueEndProcess);
          output.FillVectorVar(recovtx_dau_recoendproc, reconEndProcess);

          // Calculate angle between daughter direction and parent direction
          Double_t angleWithParent = -999.0;
          if (reconVertex->Parent && daughter->DirectionStart[0] != -999 &&
              daughter->DirectionStart[1] != -999 && daughter->DirectionStart[2] != -999) {

            // Get parent direction at end (if available)
            Float_t parentDirection[3];
            if (reconVertex->Parent->DirectionEnd[0] != -999 &&
                reconVertex->Parent->DirectionEnd[1] != -999 &&
                reconVertex->Parent->DirectionEnd[2] != -999) {
              parentDirection[0] = reconVertex->Parent->DirectionEnd[0];
              parentDirection[1] = reconVertex->Parent->DirectionEnd[1];
              parentDirection[2] = reconVertex->Parent->DirectionEnd[2];
            } else {
              // Fallback: calculate direction from parent start to end
              parentDirection[0] = reconVertex->Parent->PositionEnd[0] - reconVertex->Parent->PositionStart[0];
              parentDirection[1] = reconVertex->Parent->PositionEnd[1] - reconVertex->Parent->PositionStart[1];
              parentDirection[2] = reconVertex->Parent->PositionEnd[2] - reconVertex->Parent->PositionStart[2];

              // Normalize parent direction
              Float_t parentNorm = sqrt(parentDirection[0]*parentDirection[0] +
                                       parentDirection[1]*parentDirection[1] +
                                       parentDirection[2]*parentDirection[2]);
              if (parentNorm > 0) {
                parentDirection[0] /= parentNorm;
                parentDirection[1] /= parentNorm;
                parentDirection[2] /= parentNorm;
              }
            }

            // Calculate dot product between daughter and parent directions
            Float_t dotProduct = daughter->DirectionStart[0] * parentDirection[0] +
                                daughter->DirectionStart[1] * parentDirection[1] +
                                daughter->DirectionStart[2] * parentDirection[2];

            // Clamp dot product to avoid numerical issues
            dotProduct = std::max(-1.0f, std::min(1.0f, dotProduct));
            angleWithParent = acos(dotProduct) * 180.0 / M_PI; // Convert to degrees
          }
          output.FillVectorVar(recovtx_dau_par_recoangles, angleWithParent);

          // Calculate true angle between daughter and parent
          Double_t trueAngleWithParent = -999.0;
          if (trueDaughter && reconVertex->Parent) {
            AnaTrueParticleB* trueParent = reconVertex->Parent->GetTrueParticle();
            if (trueParent && trueDaughter->Direction[0] != -999 &&
                trueDaughter->Direction[1] != -999 && trueDaughter->Direction[2] != -999) {

              // Get true parent direction at end
              Float_t trueParentDirection[3];
              if (trueParent->DirectionEnd[0] != -999 &&
                  trueParent->DirectionEnd[1] != -999 &&
                  trueParent->DirectionEnd[2] != -999) {
                trueParentDirection[0] = trueParent->DirectionEnd[0];
                trueParentDirection[1] = trueParent->DirectionEnd[1];
                trueParentDirection[2] = trueParent->DirectionEnd[2];
              } else {
                // Fallback: calculate direction from parent start to end
                trueParentDirection[0] = trueParent->PositionEnd[0] - trueParent->Position[0];
                trueParentDirection[1] = trueParent->PositionEnd[1] - trueParent->Position[1];
                trueParentDirection[2] = trueParent->PositionEnd[2] - trueParent->Position[2];

                // Normalize true parent direction
                Float_t trueParentNorm = sqrt(trueParentDirection[0]*trueParentDirection[0] +
                                             trueParentDirection[1]*trueParentDirection[1] +
                                             trueParentDirection[2]*trueParentDirection[2]);
                if (trueParentNorm > 0) {
                  trueParentDirection[0] /= trueParentNorm;
                  trueParentDirection[1] /= trueParentNorm;
                  trueParentDirection[2] /= trueParentNorm;
                }
              }

              // Calculate dot product between true daughter and true parent directions
              Float_t trueDotProduct = trueDaughter->Direction[0] * trueParentDirection[0] +
                                      trueDaughter->Direction[1] * trueParentDirection[1] +
                                      trueDaughter->Direction[2] * trueParentDirection[2];

              // Clamp dot product to avoid numerical issues
              trueDotProduct = std::max(-1.0f, std::min(1.0f, trueDotProduct));
              trueAngleWithParent = acos(trueDotProduct) * 180.0 / M_PI; // Convert to degrees
            }
          }
          output.FillVectorVar(recovtx_dau_par_trueangles, trueAngleWithParent);

          // Fill daughter momentum using range-momentum calculation (pion hypothesis only)
          Double_t reconMomentum = -999.0;
          Double_t trueMomentum = -999.0;

          // Calculate reconstructed momentum from track length using pion range-momentum tables
          AnaParticlePD* daughterPD = static_cast<AnaParticlePD*>(daughter);
          if (daughterPD && daughterPD->Length > 0) {
            // Try pion hypothesis (PDG 211)
            reconMomentum = pdAnaUtils::ComputeRangeMomentum(daughterPD->Length, 211);
            // If it fails, keep -999
          }
          output.FillVectorVar(recovtx_dau_recomom, reconMomentum);

          // Calculate true momentum from true particle
          if (trueDaughter && trueDaughter->Momentum > 0) {
            trueMomentum = trueDaughter->Momentum / 1000.0; // Convert from MeV to GeV
          }
          output.FillVectorVar(recovtx_dau_truemom, trueMomentum);

          // Fill daughter dE/dx using pdAnaUtils methods
          Double_t reconDedx = -999.0;
          Double_t trueDedx = -999.0;

          // Calculate reconstructed dE/dx using pdAnaUtils
          if (daughterPD) {
            reconDedx = pdAnaUtils::ComputeAveragedEdxOverResRange(daughterPD);
            // If the above fails, try the truncLibo_dEdx as fallback
            if (reconDedx <= 0 && daughterPD->truncLibo_dEdx > 0) {
              reconDedx = daughterPD->truncLibo_dEdx;
            }
          }
          output.FillVectorVar(recovtx_dau_recodedx, reconDedx);

          // Calculate true dE/dx from true particle
          if (trueDaughter) {
            // For true particles, we can calculate dE/dx from energy and length
            // or use other available information
            if (trueDaughter->Momentum > 0 && daughter->Length > 0) {
              // This is a simplified calculation - in practice you might want to use
              // more sophisticated energy loss calculations
              trueDedx = trueDaughter->Momentum / daughter->Length;
            }
          }
          output.FillVectorVar(recovtx_dau_truededx, trueDedx);
        }
      }
    }

    // Calculate invariant mass of the parent particle from daughter momenta
    Double_t recoInvariantMass = -999.0;
    Double_t trueInvariantMass = -999.0;
    const Double_t m_pion = 0.13957; // GeV, pion mass

    // Only calculate if we have exactly 2 daughters (typical for K0 -> pi+ pi- decay)
    if (reconVertex->NParticles == 2) {
      // Get the two daughter particles
      AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(reconVertex->Particles[0]);
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(reconVertex->Particles[1]);

      if (daughter1 && daughter2) {
        // Calculate reconstructed invariant mass
        Double_t p1_reco = -999.0;
        Double_t p2_reco = -999.0;

        // Get reconstructed momenta using the same method as above
        if (daughter1->Length > 0) {
          p1_reco = pdAnaUtils::ComputeRangeMomentum(daughter1->Length, 211);
        }
        if (daughter2->Length > 0) {
          p2_reco = pdAnaUtils::ComputeRangeMomentum(daughter2->Length, 211);
        }

        // Calculate reconstructed invariant mass if both momenta are valid
        if (p1_reco > 0 && p2_reco > 0) {
          // For two pions, invariant mass = sqrt(2 * (E1*E2 - p1*p2*cos(theta)))
          // where E = sqrt(p^2 + m_pion^2) and theta is the angle between momenta
          Double_t E1 = sqrt(p1_reco*p1_reco + m_pion*m_pion);
          Double_t E2 = sqrt(p2_reco*p2_reco + m_pion*m_pion);

          // Calculate angle between daughter momenta
          Double_t cos_theta = 0.0;
          if (daughter1->DirectionStart[0] != -999 && daughter1->DirectionStart[1] != -999 && daughter1->DirectionStart[2] != -999 &&
              daughter2->DirectionStart[0] != -999 && daughter2->DirectionStart[1] != -999 && daughter2->DirectionStart[2] != -999) {
            cos_theta = daughter1->DirectionStart[0] * daughter2->DirectionStart[0] +
                       daughter1->DirectionStart[1] * daughter2->DirectionStart[1] +
                       daughter1->DirectionStart[2] * daughter2->DirectionStart[2];
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta)); // Clamp to avoid numerical issues
          }

          // Invariant mass formula: M^2 = 2 * (E1*E2 - p1*p2*cos(theta))
          Double_t M_squared = 2.0 * (E1*E2 - p1_reco*p2_reco*cos_theta);
          if (M_squared > 0) {
            recoInvariantMass = sqrt(M_squared);
          }
        }

        // Calculate true invariant mass independently
        AnaTrueParticleB* trueDaughter1 = daughter1->GetTrueParticle();
        AnaTrueParticleB* trueDaughter2 = daughter2->GetTrueParticle();

        if (trueDaughter1 && trueDaughter2 &&
            trueDaughter1->Momentum > 0 && trueDaughter2->Momentum > 0) {

          Double_t p1_true = trueDaughter1->Momentum / 1000.0; // Convert from MeV to GeV
          Double_t p2_true = trueDaughter2->Momentum / 1000.0; // Convert from MeV to GeV

          Double_t E1_true = sqrt(p1_true*p1_true + m_pion*m_pion);
          Double_t E2_true = sqrt(p2_true*p2_true + m_pion*m_pion);

          // Calculate angle between true daughter momenta
          Double_t cos_theta_true = 0.0;
          if (trueDaughter1->Direction[0] != -999 && trueDaughter1->Direction[1] != -999 && trueDaughter1->Direction[2] != -999 &&
              trueDaughter2->Direction[0] != -999 && trueDaughter2->Direction[1] != -999 && trueDaughter2->Direction[2] != -999) {
            cos_theta_true = trueDaughter1->Direction[0] * trueDaughter2->Direction[0] +
                            trueDaughter1->Direction[1] * trueDaughter2->Direction[1] +
                            trueDaughter1->Direction[2] * trueDaughter2->Direction[2];
            cos_theta_true = std::max(-1.0, std::min(1.0, cos_theta_true)); // Clamp to avoid numerical issues
          }

          // Invariant mass formula: M^2 = 2 * (E1*E2 - p1*p2*cos(theta))
          Double_t M_squared_true = 2.0 * (E1_true*E2_true - p1_true*p2_true*cos_theta_true);
          if (M_squared_true > 0) {
            trueInvariantMass = sqrt(M_squared_true);
          }
        }
      }
    }

    output.FillVectorVar(recovtx_recomass, recoInvariantMass);
    output.FillVectorVar(recovtx_truemass, trueInvariantMass);

    output.IncrementCounter(n_recovtx_candidates);
  }

  // Fill true vertex candidates
  for (const auto& trueVertex : trueCandidates) {
    if (!trueVertex) continue;

    output.FillVectorVar(true_vertex_nparticles, (Int_t)trueVertex->NTrueParticles);
    output.FillMatrixVarFromArray(true_vertex_position, trueVertex->Position, 3);
    output.FillVectorVar(true_vertex_quality, (Double_t)trueVertex->NTrueParticles); // Simple quality metric based on number of particles
    output.IncrementCounter(n_true_vertex_candidates);
  }

  // // Fill the count variables directly
  output.FillVar(n_recovtx_candidates, (Int_t)reconCandidates.size());
  // output.FillVar(n_true_vertex_candidates, (Int_t)trueCandidates.size());
}
