#include "pdNeutralUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdKalman.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"
#include "TMinuit.h"
#include "TVector3.h"
#include "TMath.h"
#include <algorithm>
#include <set>
#include <unordered_map>
#include <iostream>

// Namespace for internal structures and functions
namespace {
  // Structure to hold hit data for TMinuit FCN
  struct VertexFitData {
    std::vector<TVector3> hits1;  // Hits from first track
    std::vector<TVector3> hits2;  // Hits from second track
  };

  // Global pointer to fit data (needed for TMinuit callback)
  VertexFitData* gFitData = nullptr;

  //***************************************************************
  // TMinuit FCN function for vertex fitting
  //***************************************************************
  void VertexFitFCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // Suppress compiler warnings
    (void)npar;
    (void)gin;
    (void)iflag;

    if (!gFitData) {
      f = 1e10;
      return;
    }

    // Extract parameters
    // par[0-2]: vertex position (vx, vy, vz)
    // par[3-5]: direction 1 (dx1, dy1, dz1)
    // par[6-8]: direction 2 (dx2, dy2, dz2)
    TVector3 vertex(par[0], par[1], par[2]);
    TVector3 dir1(par[3], par[4], par[5]);
    TVector3 dir2(par[6], par[7], par[8]);

    // Normalize direction vectors
    if (dir1.Mag() > 0) dir1 = dir1.Unit();
    else {
      f = 1e10;
      return;
    }

    if (dir2.Mag() > 0) dir2 = dir2.Unit();
    else {
      f = 1e10;
      return;
    }

    // Calculate sum of squared distances along track direction
    double chi2 = 0.0;

    // Process hits from track 1
    for (const auto& hit : gFitData->hits1) {
      // Vector from vertex to hit
      TVector3 toHit = hit - vertex;

      // Project onto direction to find parameter t
      double t = toHit.Dot(dir1);

      // Calculate projection point on line
      TVector3 projection = vertex + t * dir1;

      // Calculate squared distance
      double dist2 = (hit - projection).Mag2();
      chi2 += dist2;
    }

    // Process hits from track 2
    for (const auto& hit : gFitData->hits2) {
      // Vector from vertex to hit
      TVector3 toHit = hit - vertex;

      // Project onto direction to find parameter t
      double t = toHit.Dot(dir2);

      // Calculate projection point on line
      TVector3 projection = vertex + t * dir2;

      // Calculate squared distance
      double dist2 = (hit - projection).Mag2();
      chi2 += dist2;
    }

    f = chi2;
  }
} // anonymous namespace

//***************************************************************
// Helper function to calculate Pandora vertex position from two particles
// Returns: 1 if used simple average (else branch), 0 if used line intersection
//***************************************************************
int CalculatePandoraVertexPosition(AnaParticlePD* particle1, AnaParticlePD* particle2, Float_t* position) {
    if (!particle1 || !particle2) {
        for (int i = 0; i < 3; i++) position[i] = -999.;
        return -999; // Default to "just average" for invalid inputs
    }

    TVector3 pos1(particle1->PositionStart[0], particle1->PositionStart[1], particle1->PositionStart[2]);
    TVector3 dir1(particle1->DirectionStart[0], particle1->DirectionStart[1], particle1->DirectionStart[2]);
    TVector3 pos2(particle2->PositionStart[0], particle2->PositionStart[1], particle2->PositionStart[2]);
    TVector3 dir2(particle2->DirectionStart[0], particle2->DirectionStart[1], particle2->DirectionStart[2]);

    // Normalize directions
    if (dir1.Mag() > 0) dir1 = dir1.Unit();
    if (dir2.Mag() > 0) dir2 = dir2.Unit();

    // Find closest point between two lines
    TVector3 w0 = pos1 - pos2;
    double a = dir1.Dot(dir1);
    double b = dir1.Dot(dir2);
    double c = dir2.Dot(dir2);
    double d = dir1.Dot(w0);
    double e = dir2.Dot(w0);
    double denom = a*c - b*b;

    TVector3 pandoraVertex;
    int isJustAverage;
    if (fabs(denom) > 1e-6) {
        double s = (b*e - c*d) / denom;
        double t = (a*e - b*d) / denom;
        TVector3 p1 = pos1 + s * dir1;
        TVector3 p2 = pos2 + t * dir2;
        pandoraVertex = 0.5 * (p1 + p2);
        isJustAverage = 0; // Used line intersection
    } else {
        pandoraVertex = 0.5 * (pos1 + pos2);
        isJustAverage = 1; // Used simple average (else branch)
    }

    position[0] = pandoraVertex.X();
    position[1] = pandoraVertex.Y();
    position[2] = pandoraVertex.Z();

    return isJustAverage;
}

//***************************************************************
// Find improved neutral particle start position by looking for nearby particles
//***************************************************************
TVector3 FindNeutralParticleStartPosition(
    AnaNeutralParticlePD* neutralParticle,
    AnaParticlePD* parentParticle,
    AnaEventB& event,
    double proximityThreshold) {

    // Default: use parent end position
    TVector3 defaultStart(parentParticle->PositionEnd[0],
                          parentParticle->PositionEnd[1],
                          parentParticle->PositionEnd[2]);

    if (!neutralParticle || !parentParticle || !neutralParticle->Vertex) {
        return defaultStart;
    }

    // Get parent particle hits
    std::vector<TVector3> parentHits;
    for (int plane = 0; plane < 3; plane++) {
        for (size_t h = 0; h < parentParticle->Hits[plane].size(); h++) {
            AnaHitPD hit = parentParticle->Hits[plane][h];
            if (hit.Position.X() > -900) {
                parentHits.push_back(hit.Position);
            }
        }
    }

    if (parentHits.empty()) {
        return defaultStart;
    }

    // Get vertex daughter particles to exclude
    std::set<AnaParticleB*> vertexParticles;
    for (int i = 0; i < neutralParticle->Vertex->NParticles; i++) {
        if (neutralParticle->Vertex->Particles[i]) {
            vertexParticles.insert(neutralParticle->Vertex->Particles[i]);
        }
    }

    // Find nearby particles
    std::vector<AnaParticlePD*> nearbyParticles;
    AnaParticleB** allParticles = event.Particles;
    int nParticles = event.nParticles;

    for (int i = 0; i < nParticles; i++) {
        AnaParticlePD* particle = static_cast<AnaParticlePD*>(allParticles[i]);
        if (!particle) continue;

        // Skip if this particle is in the vertex
        if (vertexParticles.find(particle) != vertexParticles.end()) {
            continue;
        }

        // Check if particle start position is close to any parent hit
        TVector3 particleStart(particle->PositionStart[0],
                               particle->PositionStart[1],
                               particle->PositionStart[2]);

        if (particleStart.X() < -900) continue;

        bool isNearby = false;
        for (const auto& hit : parentHits) {
            double distance = (particleStart - hit).Mag();
            if (distance < proximityThreshold) {
                isNearby = true;
                break;
            }
        }

        if (isNearby) {
            nearbyParticles.push_back(particle);
        }
    }

    if (nearbyParticles.empty()) {
        return defaultStart;
    }

    // Find the particle with minimum distance to parent track
    TVector3 parentPos(parentParticle->PositionStart[0],
                       parentParticle->PositionStart[1],
                       parentParticle->PositionStart[2]);
    TVector3 parentDir(parentParticle->DirectionStart[0],
                       parentParticle->DirectionStart[1],
                       parentParticle->DirectionStart[2]);
    if (parentDir.Mag() > 0) parentDir = parentDir.Unit();

    double minDistance = 999999.0;
    TVector3 bestNeutralStart = defaultStart;

    for (auto* particle : nearbyParticles) {
        TVector3 particlePos(particle->PositionStart[0],
                             particle->PositionStart[1],
                             particle->PositionStart[2]);
        TVector3 particleDir(particle->DirectionStart[0],
                             particle->DirectionStart[1],
                             particle->DirectionStart[2]);

        if (particleDir.Mag() > 0) particleDir = particleDir.Unit();

        // Find closest approach between two lines (same algorithm as vertex finding)
        TVector3 w0 = parentPos - particlePos;
        double a = parentDir.Dot(parentDir);
        double b = parentDir.Dot(particleDir);
        double c = particleDir.Dot(particleDir);
        double d = parentDir.Dot(w0);
        double e = particleDir.Dot(w0);
        double denom = a*c - b*b;

        if (fabs(denom) > 1e-6) {
            double s = (b*e - c*d) / denom;
            double t = (a*e - b*d) / denom;

            TVector3 p1 = parentPos + s * parentDir;
            TVector3 p2 = particlePos + t * particleDir;

            // Calculate distance between closest points
            double distance = (p1 - p2).Mag();

            // Keep track of particle with minimum distance
            if (distance < minDistance && distance < proximityThreshold) {
                minDistance = distance;
                // Use average of closest points as neutral start position
                bestNeutralStart = 0.5 * (p1 + p2);
            }
        }
    }

    return bestNeutralStart;
}

//***************************************************************
// Helper function to calculate degeneracy for all vertices in the list
//***************************************************************
void CalculateVertexDegeneracy(std::vector<AnaVertexPD*>& vertices, bool usePandora, double threshold) {
    for (size_t i = 0; i < vertices.size(); i++) {
        AnaVertexPD* vtx1 = vertices[i];
        if (!vtx1) continue;

        // Get position to compare (Pandora or fitted)
        TVector3 pos1;
        if (usePandora) {
            pos1.SetXYZ(vtx1->PositionPandora[0], vtx1->PositionPandora[1], vtx1->PositionPandora[2]);
        } else {
            pos1.SetXYZ(vtx1->Position[0], vtx1->Position[1], vtx1->Position[2]);
        }

        // Skip if position is invalid
        if (pos1.X() < -900) continue;

        int degeneracy = 0;
        std::set<AnaParticleB*> uniqueParticles;
        std::vector<double> distances; // Store all distances within threshold

        // Add this vertex's particles to the set
        for (int p = 0; p < vtx1->NParticles; p++) {
            if (vtx1->Particles[p]) {
                uniqueParticles.insert(vtx1->Particles[p]);
            }
        }

        // Compare with all other vertices
        for (size_t j = 0; j < vertices.size(); j++) {
            if (i == j) continue;

            AnaVertexPD* vtx2 = vertices[j];
            if (!vtx2) continue;

            TVector3 pos2;
            if (usePandora) {
                pos2.SetXYZ(vtx2->PositionPandora[0], vtx2->PositionPandora[1], vtx2->PositionPandora[2]);
            } else {
                pos2.SetXYZ(vtx2->Position[0], vtx2->Position[1], vtx2->Position[2]);
            }

            if (pos2.X() < -900) continue;

            // Check if within threshold distance
            double distance = (pos1 - pos2).Mag();
            if (distance < threshold) {
                degeneracy++;
                distances.push_back(distance); // Store the distance

                // Add this vertex's particles to unique set
                for (int p = 0; p < vtx2->NParticles; p++) {
                    if (vtx2->Particles[p]) {
                        uniqueParticles.insert(vtx2->Particles[p]);
                    }
                }
            }
        }

        vtx1->DegeneracyBeforeScoring = degeneracy;
        vtx1->NRecoParticles = uniqueParticles.size();

        // Sort distances and store the 5 minimum values
        vtx1->DegeneracyDistances.clear();
        if (!distances.empty()) {
            std::sort(distances.begin(), distances.end());
            size_t nToStore = std::min(distances.size(), size_t(5));
            for (size_t k = 0; k < nToStore; k++) {
                vtx1->DegeneracyDistances.push_back(static_cast<Float_t>(distances[k]));
            }
        }
    }
}

//***************************************************************
// Helper struct to store particle info with distances
//***************************************************************
struct IsolationParticleInfo {
    AnaParticlePD* particle;
    double pandoraDistance;
    double fitDistance;
    double startDistance;
};

//***************************************************************
// Helper struct to store creation vertex particle info
//***************************************************************
struct CreationVtxParticleInfo {
    AnaParticlePD* particle;
    double pandoraDistance;
    int truePDG;
    float chi2Proton;
};

//***************************************************************
// Helper function to count and characterize particles near neutral particle creation vertex
//***************************************************************
void CalculateNProtonInCreationVtx(AnaNeutralParticlePD* neutralParticle,
                                   AnaParticlePD* parentParticle,
                                   double protonCreationVtxDist) {
    // Initialize all attributes
    neutralParticle->NProtonInCreationVtx = 0;
    neutralParticle->NParticlesInCreationVtx = 0;
    neutralParticle->CreationVtxChi2Proton.clear();
    neutralParticle->CreationVtxDistances.clear();
    neutralParticle->CreationVtxTruePDG.clear();

    if (!parentParticle || !neutralParticle || neutralParticle->PositionStart[0] < -900) {
        return;
    }

    TVector3 neutralStart(neutralParticle->PositionStart[0],
                         neutralParticle->PositionStart[1],
                         neutralParticle->PositionStart[2]);

    std::vector<CreationVtxParticleInfo> particlesNearCreation;
    int nProtons = 0;

    // Iterate through parent's daughters
    for (size_t iDau = 0; iDau < parentParticle->Daughters.size(); iDau++) {
        AnaParticlePD* daughter = static_cast<AnaParticlePD*>(parentParticle->Daughters[iDau]);
        if (!daughter) continue;

        // Check if daughter has valid position and direction
        if (daughter->PositionStart[0] < -900 || daughter->DirectionStart[0] < -900) continue;

        // Calculate line-to-point distance from daughter's trajectory to neutral start position
        // Line: L(t) = L0 + t * L_dir (daughter trajectory)
        // Point: P0 (neutral start position)
        // Distance: d = ||(P0 - L0) × L_dir|| / ||L_dir||
        TVector3 L0(daughter->PositionStart[0], daughter->PositionStart[1], daughter->PositionStart[2]);
        TVector3 L_dir(daughter->DirectionStart[0], daughter->DirectionStart[1], daughter->DirectionStart[2]);

        double L_dir_mag = L_dir.Mag();
        if (L_dir_mag < 1e-6) continue; // Skip if direction is too small

        TVector3 diff = neutralStart - L0;
        TVector3 cross = diff.Cross(L_dir);
        double minDistance = cross.Mag() / L_dir_mag;

        // Check if distance is below threshold
        if (minDistance < protonCreationVtxDist) {
            // Count true protons
            int truePDG = 0;
            bool isTrueProton = false;
            if (daughter->TrueObject) {
                AnaTrueParticlePD* trueDau = static_cast<AnaTrueParticlePD*>(daughter->TrueObject);
                if (trueDau) {
                    truePDG = trueDau->PDG;
                    if (trueDau->PDG == 2212) {
                        isTrueProton = true;
                        nProtons++;
                    }
                }
            }

            // Calculate chi2/ndf under proton hypothesis
            std::pair<double, int> chi2Result = pdAnaUtils::Chi2PID(*daughter, 2212); // 2212 = proton PDG
            float chi2Proton = (chi2Result.second > 0) ?
                              static_cast<float>(chi2Result.first) / static_cast<float>(chi2Result.second) : -999.0;

            // Store particle info
            CreationVtxParticleInfo info;
            info.particle = daughter;
            info.pandoraDistance = minDistance;
            info.truePDG = truePDG;
            info.chi2Proton = chi2Proton;
            particlesNearCreation.push_back(info);
        }
    }

    // Sort particles by Pandora distance (ascending)
    std::sort(particlesNearCreation.begin(), particlesNearCreation.end(),
              [](const CreationVtxParticleInfo& a, const CreationVtxParticleInfo& b) {
                  return a.pandoraDistance < b.pandoraDistance;
              });

    // Store counts
    neutralParticle->NProtonInCreationVtx = nProtons;
    neutralParticle->NParticlesInCreationVtx = static_cast<int>(particlesNearCreation.size());

    // Store properties of up to 5 closest particles
    size_t nToStore = std::min(particlesNearCreation.size(), size_t(5));
    for (size_t k = 0; k < nToStore; k++) {
        neutralParticle->CreationVtxDistances.push_back(static_cast<Float_t>(particlesNearCreation[k].pandoraDistance));
        neutralParticle->CreationVtxChi2Proton.push_back(particlesNearCreation[k].chi2Proton);
        neutralParticle->CreationVtxTruePDG.push_back(particlesNearCreation[k].truePDG);
    }
}

//***************************************************************
// Helper function to calculate isolation distances for all vertices
//***************************************************************
void CalculateVertexIsolation(std::vector<AnaVertexPD*>& vertices, double threshold) {
    // Get track fit length from parameters
    double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

    for (size_t i = 0; i < vertices.size(); i++) {
        AnaVertexPD* vtx1 = vertices[i];
        if (!vtx1) continue;

        // Get vertex position (now using Pandora-based position as main Position)
        TVector3 vtxPos(vtx1->Position[0], vtx1->Position[1], vtx1->Position[2]);

        // Skip if position is invalid
        if (vtxPos.X() < -900) continue;

        // Get lengths of the two vertex particles for comparison
        double vtxPart1Length = -999.0;
        double vtxPart2Length = -999.0;
        if (vtx1->NParticles >= 2) {
            AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vtx1->Particles[0]);
            AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vtx1->Particles[1]);
            if (part1) vtxPart1Length = part1->Length;
            if (part2) vtxPart2Length = part2->Length;
        }

        std::vector<IsolationParticleInfo> isolationParticles;
        std::set<AnaParticleB*> vtx1Particles;

        // Store this vertex's particles to exclude them
        for (int p = 0; p < vtx1->NParticles; p++) {
            if (vtx1->Particles[p]) {
                vtx1Particles.insert(vtx1->Particles[p]);
            }
        }

        // Initialize counters
        int nProtons = 0;
        int nPions = 0;
        int isLongest = 0;

        // Find nearby vertices and collect their unique particles
        for (size_t j = 0; j < vertices.size(); j++) {
            if (i == j) continue;

            AnaVertexPD* vtx2 = vertices[j];
            if (!vtx2) continue;

            TVector3 pos2(vtx2->Position[0], vtx2->Position[1], vtx2->Position[2]);
            if (pos2.X() < -900) continue;

            // Check if within threshold distance
            double distance = (vtxPos - pos2).Mag();
            if (distance < threshold) {
                // Process particles from this nearby vertex
                for (int p = 0; p < vtx2->NParticles; p++) {
                    AnaParticleB* particle = vtx2->Particles[p];
                    if (!particle) continue;

                    // Skip if this particle belongs to vtx1
                    if (vtx1Particles.find(particle) != vtx1Particles.end()) continue;

                    AnaParticlePD* particlePD = static_cast<AnaParticlePD*>(particle);

                    // Check if particle has valid position and direction
                    if (particlePD->PositionStart[0] < -900 || particlePD->DirectionStart[0] < -900) continue;

                    TVector3 L0(particlePD->PositionStart[0], particlePD->PositionStart[1], particlePD->PositionStart[2]);

                    // === 1. Calculate Pandora-based line-to-point distance ===
                    TVector3 L_dir(particlePD->DirectionStart[0], particlePD->DirectionStart[1], particlePD->DirectionStart[2]);

                    double L_dir_mag = L_dir.Mag();
                    if (L_dir_mag < 1e-6) continue; // Skip if direction is too small

                    TVector3 diff = vtxPos - L0;
                    TVector3 cross = diff.Cross(L_dir);
                    double lineToPointDist = cross.Mag() / L_dir_mag;

                    // === 2. Calculate fitted track line-to-point distance ===
                    std::vector<double> fittedLineParams;
                    pdAnaUtils::ExtrapolateTrack(particlePD, fittedLineParams, trackFitLength, true);

                    double lineToPointDistFit = 9999.0;
                    if (fittedLineParams.size() == 6 && fittedLineParams[0] != -999.0) {
                        TVector3 L0_fit(fittedLineParams[0], fittedLineParams[1], fittedLineParams[2]);
                        TVector3 L_dir_fit(fittedLineParams[3], fittedLineParams[4], fittedLineParams[5]);

                        double L_dir_fit_mag = L_dir_fit.Mag();
                        if (L_dir_fit_mag > 1e-6) {
                            TVector3 diff_fit = vtxPos - L0_fit;
                            TVector3 cross_fit = diff_fit.Cross(L_dir_fit);
                            lineToPointDistFit = cross_fit.Mag() / L_dir_fit_mag;
                        }
                    }

                    // === 3. Calculate point-to-point distance to PositionStart ===
                    double pointToPointDist = (vtxPos - L0).Mag();

                    // Store particle info
                    IsolationParticleInfo info;
                    info.particle = particlePD;
                    info.pandoraDistance = lineToPointDist;
                    info.fitDistance = lineToPointDistFit;
                    info.startDistance = pointToPointDist;
                    isolationParticles.push_back(info);

                    // === Count protons and pions (ALL particles, not just 5) ===
                    if (particlePD->TrueObject) {
                        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(particlePD->TrueObject);
                        if (truePart) {
                            if (truePart->PDG == 2212) {
                                nProtons++;
                            }
                            if (abs(truePart->PDG) == 211) { // Both pi+ and pi-
                                nPions++;
                            }
                        }
                    }

                    // === Check if longer than BOTH vertex particles ===
                    if (particlePD->Length > vtxPart1Length && particlePD->Length > vtxPart2Length) {
                        isLongest = 1;
                    }
                }
            }
        }

        // Sort particles by Pandora distance
        std::sort(isolationParticles.begin(), isolationParticles.end(),
                  [](const IsolationParticleInfo& a, const IsolationParticleInfo& b) {
                      return a.pandoraDistance < b.pandoraDistance;
                  });

        // Clear all isolation vectors
        vtx1->IsolationDistances.clear();
        vtx1->IsolationDistancesFit.clear();
        vtx1->IsolationStartDistances.clear();
        vtx1->IsolationIsProton.clear();
        vtx1->IsolationChi2Proton.clear();
        vtx1->IsolationLength.clear();

        // Store the 5 closest particles' properties
        size_t nToStore = std::min(isolationParticles.size(), size_t(5));
        for (size_t k = 0; k < nToStore; k++) {
            AnaParticlePD* part = isolationParticles[k].particle;

            // Store distances
            vtx1->IsolationDistances.push_back(static_cast<Float_t>(isolationParticles[k].pandoraDistance));
            vtx1->IsolationDistancesFit.push_back(static_cast<Float_t>(isolationParticles[k].fitDistance));
            vtx1->IsolationStartDistances.push_back(static_cast<Float_t>(isolationParticles[k].startDistance));

            // Check if this particle is a true proton
            int isProton = 0;
            if (part->TrueObject) {
                AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
                if (truePart && truePart->PDG == 2212) {
                    isProton = 1;
                }
            }
            vtx1->IsolationIsProton.push_back(isProton);

            // Calculate chi2/ndf under proton hypothesis
            std::pair<double, int> chi2Result = pdAnaUtils::Chi2PID(*part, 2212); // 2212 = proton PDG
            Float_t chi2Proton = static_cast<Float_t>(chi2Result.first)/static_cast<Float_t>(chi2Result.second); // Extract chi2/ndf value
            vtx1->IsolationChi2Proton.push_back(chi2Proton);

            // Store length
            vtx1->IsolationLength.push_back(part->Length);
        }

        // Store total proton count, pion count, and longest flag
        vtx1->IsolationNProton = nProtons;
        vtx1->IsolationNPion = nPions;
        vtx1->IsolationIsLongest = isLongest;
    }
}

//***************************************************************
double pdNeutralUtils::FindVertexPositionWithFit(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!part1 || !part2) {
    return -999.0;
  }

  // Get reference positions for both particles
  TVector3 pos1 = pdAnaUtils::DefinePosition(part1, true);  // Use start position
  TVector3 pos2 = pdAnaUtils::DefinePosition(part2, true);  // Use start position

  // Check if positions are valid
  if (pos1.X() < -900 || pos2.X() < -900) {
    return -999.0;
  }

  // Collect hits within trackFitLength from each particle
  VertexFitData fitData;

  // Collect hits for particle 1
  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part1->Hits[plane].size(); i++) {
      AnaHitPD& hit = part1->Hits[plane][i];
      if (hit.Position.Z() != -999) {
        TVector3 hitPos = hit.Position;
        double distance = (hitPos - pos1).Mag();
        if (distance <= trackFitLength) {
          fitData.hits1.push_back(hitPos);
        }
      }
    }
  }

  // Collect hits for particle 2
  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part2->Hits[plane].size(); i++) {
      AnaHitPD& hit = part2->Hits[plane][i];
      if (hit.Position.Z() != -999) {
        TVector3 hitPos = hit.Position;
        double distance = (hitPos - pos2).Mag();
        if (distance <= trackFitLength) {
          fitData.hits2.push_back(hitPos);
        }
      }
    }
  }

  // Need sufficient hits for fitting
  if (fitData.hits1.size() < 2 || fitData.hits2.size() < 2) {
    return -999.0;
  }

  // Set global fit data pointer
  gFitData = &fitData;

  // Initialize TMinuit with 9 parameters
  TMinuit minuit(9);
  minuit.SetPrintLevel(-1);  // Suppress output

  // Set FCN function
  minuit.SetFCN(VertexFitFCN);

  // Calculate initial vertex position as midpoint
  double vx_init = 0.5 * (pos1.X() + pos2.X());
  double vy_init = 0.5 * (pos1.Y() + pos2.Y());
  double vz_init = 0.5 * (pos1.Z() + pos2.Z());

  // Get initial directions from particles
  TVector3 dir1_init(part1->DirectionStart[0], part1->DirectionStart[1], part1->DirectionStart[2]);
  TVector3 dir2_init(part2->DirectionStart[0], part2->DirectionStart[1], part2->DirectionStart[2]);

  // Normalize initial directions
  if (dir1_init.Mag() > 0) dir1_init = dir1_init.Unit();
  if (dir2_init.Mag() > 0) dir2_init = dir2_init.Unit();

  // Set parameters
  Double_t step = 0.1;
  minuit.DefineParameter(0, "vx", vx_init, step, 0, 0);
  minuit.DefineParameter(1, "vy", vy_init, step, 0, 0);
  minuit.DefineParameter(2, "vz", vz_init, step, 0, 0);
  minuit.DefineParameter(3, "dx1", dir1_init.X(), step, 0, 0);
  minuit.DefineParameter(4, "dy1", dir1_init.Y(), step, 0, 0);
  minuit.DefineParameter(5, "dz1", dir1_init.Z(), step, 0, 0);
  minuit.DefineParameter(6, "dx2", dir2_init.X(), step, 0, 0);
  minuit.DefineParameter(7, "dy2", dir2_init.Y(), step, 0, 0);
  minuit.DefineParameter(8, "dz2", dir2_init.Z(), step, 0, 0);

  // Run MIGRAD minimization
  minuit.Migrad();

  // Get fitted parameters
  Double_t vx, vy, vz, dx1, dy1, dz1, dx2, dy2, dz2;
  Double_t err;
  minuit.GetParameter(0, vx, err);
  minuit.GetParameter(1, vy, err);
  minuit.GetParameter(2, vz, err);
  minuit.GetParameter(3, dx1, err);
  minuit.GetParameter(4, dy1, err);
  minuit.GetParameter(5, dz1, err);
  minuit.GetParameter(6, dx2, err);
  minuit.GetParameter(7, dy2, err);
  minuit.GetParameter(8, dz2, err);

  // Get chi2 value
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  // Store fitted vertex position in PositionFit
  vertex->PositionFit[0] = vx;
  vertex->PositionFit[1] = vy;
  vertex->PositionFit[2] = vz;

  // Calculate Pandora-based vertex position and set flag
  vertex->IsJustAverage = CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

  // Copy Pandora position to main Position
  vertex->Position[0] = vertex->PositionPandora[0];
  vertex->Position[1] = vertex->PositionPandora[1];
  vertex->Position[2] = vertex->PositionPandora[2];

  // Calculate chi2/ndf
  int ndf = (fitData.hits1.size() + fitData.hits2.size()) * 3 - 9;  // 3 coords per hit, 9 parameters
  double chi2ndf = (ndf > 0) ? amin / ndf : amin;

  // Store score (chi2/ndf from minimization)
  vertex->Score = chi2ndf;

  // Calculate minimum distance between fitted lines
  TVector3 vertexPos(vx, vy, vz);
  TVector3 fittedDir1(dx1, dy1, dz1);
  TVector3 fittedDir2(dx2, dy2, dz2);

  // Normalize directions
  if (fittedDir1.Mag() > 0) fittedDir1 = fittedDir1.Unit();
  if (fittedDir2.Mag() > 0) fittedDir2 = fittedDir2.Unit();

  // Use cross product method to find minimum distance
  TVector3 w0 = vertexPos - vertexPos;  // Both lines start at vertex
  double a = fittedDir1.Dot(fittedDir1);
  double b = fittedDir1.Dot(fittedDir2);
  double c = fittedDir2.Dot(fittedDir2);
  double d = fittedDir1.Dot(w0);
  double e = fittedDir2.Dot(w0);

  double denom = a * c - b * b;
  vertex->MinimumDistance = 0.0;  // Both lines emanate from same vertex point

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = fittedDir1 + fittedDir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->DirectionFit[0] = fitDirSum.X();
  vertex->DirectionFit[1] = fitDirSum.Y();
  vertex->DirectionFit[2] = fitDirSum.Z();

  // Clean up global pointer
  gFitData = nullptr;

  return chi2ndf;
}

//***************************************************************
std::vector<AnaVertexPD*> pdNeutralUtils::FilterVerticesByScore(std::vector<AnaVertexPD*>& vertices) {
//***************************************************************

  // Sort vertices by Score (ascending - lower is better)
  std::sort(vertices.begin(), vertices.end(),
    [](const AnaVertexPD* a, const AnaVertexPD* b) {
      // Handle invalid scores
      if (a->Score < -900 && b->Score < -900) return false;
      if (a->Score < -900) return false;  // Invalid score goes to back
      if (b->Score < -900) return true;   // Invalid score goes to back
      return a->Score < b->Score;
    });

  // Track which particles have been used
  std::set<AnaParticlePD*> usedParticles;

  // Output vector for selected vertices
  std::vector<AnaVertexPD*> selectedVertices;
  std::set<AnaVertexPD*> selectedSet;  // For quick lookup

  // Iterate through sorted vertices
  for (AnaVertexPD* vertex : vertices) {
    if (!vertex || vertex->NParticles < 2) continue;

    AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
    AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

    if (!part1 || !part2) continue;

    // Check if either particle is already used
    if (usedParticles.find(part1) == usedParticles.end() &&
        usedParticles.find(part2) == usedParticles.end()) {
      // Neither particle is used, so we can keep this vertex
      selectedVertices.push_back(vertex);
      selectedSet.insert(vertex);

      // Mark particles as used
      usedParticles.insert(part1);
      usedParticles.insert(part2);
    }
  }

  // Delete vertices that were not selected
  for (AnaVertexPD* vertex : vertices) {
    if (selectedSet.find(vertex) == selectedSet.end()) {
      delete vertex;
    }
  }

  return selectedVertices;
}

//***************************************************************
double pdNeutralUtils::FindVertexPositionGeometric(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  (void)trackFitLength; // Not used in geometric method

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // Use existing geometric position finding (stores result in vertex->Position)
  pdAnaUtils::FindVertexPosition(vertex);

  // Store the geometric fit result in PositionFit
  vertex->PositionFit[0] = vertex->Position[0];
  vertex->PositionFit[1] = vertex->Position[1];
  vertex->PositionFit[2] = vertex->Position[2];

  // Set score to minimum distance (lower is better)
  vertex->Score = vertex->MinimumDistance;

  // Calculate fitted direction from the stored line parameters
  if(vertex->FittedLineParams.size() >= 2){
    // Extract line parameters for daughter 1
    std::vector<double> line1 = vertex->FittedLineParams[0];
    TVector3 dir1(line1[3], line1[4], line1[5]);
    if(dir1.Mag() > 0) dir1 = dir1.Unit();

    // Extract line parameters for daughter 2
    std::vector<double> line2 = vertex->FittedLineParams[1];
    TVector3 dir2(line2[3], line2[4], line2[5]);
    if(dir2.Mag() > 0) dir2 = dir2.Unit();

    // Calculate vertex fit direction (sum of daughter fit directions, normalized)
    TVector3 fitDirSum = dir1 + dir2;
    if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
    vertex->DirectionFit[0] = fitDirSum.X();
    vertex->DirectionFit[1] = fitDirSum.Y();
    vertex->DirectionFit[2] = fitDirSum.Z();
  }

  // Calculate Pandora-based vertex position and set flag
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
  vertex->IsJustAverage = CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

  // Copy Pandora position to main Position
  vertex->Position[0] = vertex->PositionPandora[0];
  vertex->Position[1] = vertex->PositionPandora[1];
  vertex->Position[2] = vertex->PositionPandora[2];

  return vertex->MinimumDistance;
}

//***************************************************************
double pdNeutralUtils::FindVertexPositionKalman(AnaVertexPD* vertex, double trackFitLength) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return -999.0;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!daughter1 || !daughter2) {
    return -999.0;
  }

  // Check if daughter1 and daughter2 are close enough
  TVector3 pos1 = pdAnaUtils::DefinePosition(daughter1);
  TVector3 pos2 = pdAnaUtils::DefinePosition(daughter2);

  // Check if positions are valid
  if (pos1.X() < -900 || pos2.X() < -900) {
    return -999.0;
  }

  // Get seed position from geometric approach (midpoint)
  TVector3 seedPosition = 0.5 * (pos1 + pos2);

  // Convert particles to track states
  pdKalman::TrackState track1 = pdKalman::ParticleToTrackState(daughter1, trackFitLength, true);
  pdKalman::TrackState track2 = pdKalman::ParticleToTrackState(daughter2, trackFitLength, true);

  // Perform Kalman vertex fit
  pdKalman::VertexState vertexState = pdKalman::FitVertex(track1, track2, seedPosition);

  // Check if fit quality is reasonable
  double maxKalmanChi2Ndf = ND::params().GetParameterD("neutralKaonAnalysis.MaxKalmanChi2Ndf");
  double chi2ndf = (vertexState.ndf > 0) ? vertexState.chi2 / vertexState.ndf : 9999.0;
  if (chi2ndf > maxKalmanChi2Ndf) {
    // Mark vertex as invalid by setting fitted position to -999
    vertex->PositionFit[0] = -999.0;
    vertex->PositionFit[1] = -999.0;
    vertex->PositionFit[2] = -999.0;
    // Also mark main Position as invalid
    vertex->Position[0] = -999.0;
    vertex->Position[1] = -999.0;
    vertex->Position[2] = -999.0;
    vertex->Score = 9999.0;
    return 9999.0;
  }

  // Set fitted vertex position in PositionFit
  vertex->PositionFit[0] = vertexState.position.X();
  vertex->PositionFit[1] = vertexState.position.Y();
  vertex->PositionFit[2] = vertexState.position.Z();

  // Calculate Pandora-based vertex position and set flag
  vertex->IsJustAverage = CalculatePandoraVertexPosition(daughter1, daughter2, vertex->PositionPandora);

  // Copy Pandora position to main Position
  vertex->Position[0] = vertex->PositionPandora[0];
  vertex->Position[1] = vertex->PositionPandora[1];
  vertex->Position[2] = vertex->PositionPandora[2];

  // Set score to chi2 (lower is better)
  vertex->Score = vertexState.chi2;

  // Calculate minimum distance between fitted tracks at vertex
  TVector3 vtxPos = vertexState.position;
  TVector3 track1ToVtx = vtxPos - track1.position;
  TVector3 track2ToVtx = vtxPos - track2.position;
  double dist1 = (track1ToVtx - track1.direction * track1ToVtx.Dot(track1.direction)).Mag();
  double dist2 = (track2ToVtx - track2.direction * track2ToVtx.Dot(track2.direction)).Mag();
  vertex->MinimumDistance = sqrt(dist1*dist1 + dist2*dist2);

  // Calculate vertex fit direction from Kalman track states
  TVector3 dir1 = track1.direction.Unit();
  TVector3 dir2 = track2.direction.Unit();

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = dir1 + dir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->DirectionFit[0] = fitDirSum.X();
  vertex->DirectionFit[1] = fitDirSum.Y();
  vertex->DirectionFit[2] = fitDirSum.Z();


  return vertexState.chi2;
}

//***************************************************************
bool pdNeutralUtils::ValidateVertex(AnaVertexPD* vertex) {
//***************************************************************

  if (!vertex) return false;

  // Check valid position (> -900 for all coordinates)
  if (vertex->Position[0] < -900 || vertex->Position[1] < -900 || vertex->Position[2] < -900) {
    return false;
  }

  // Check valid score
  if (vertex->Score < -900 || vertex->Score > 1e6) {
    return false;
  }

  return true;
}

//***************************************************************
AnaTrueEquivalentVertexPD* pdNeutralUtils::FillTrueEquivalentVertex(AnaVertexPD* vertex) {
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return nullptr;
  }

  AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!daughter1 || !daughter2) {
    return nullptr;
  }

  // Create true equivalent vertex
  AnaTrueEquivalentVertexPD* trueEquivalentVertex = new AnaTrueEquivalentVertexPD();

  AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);

  trueEquivalentVertex->TrueParticles.push_back(trueParticle1);
  trueEquivalentVertex->TrueParticles.push_back(trueParticle2);

  // Calculate original distance between true particles
  if (trueParticle1 && trueParticle2) {
    float trueDistance = sqrt(pow(trueParticle1->Position[0] - trueParticle2->Position[0], 2) +
                              pow(trueParticle1->Position[1] - trueParticle2->Position[1], 2) +
                              pow(trueParticle1->Position[2] - trueParticle2->Position[2], 2));
    trueEquivalentVertex->OriginalDistance = trueDistance;
  } else {
    trueEquivalentVertex->OriginalDistance = -999.0;
  }

  // Initialize position and distance
  trueEquivalentVertex->Position[0] = -999.0;
  trueEquivalentVertex->Position[1] = -999.0;
  trueEquivalentVertex->Position[2] = -999.0;
  trueEquivalentVertex->MinimumDistance = -999.0;

  // Calculate true vertex position
  pdAnaUtils::FindVertexPosition(trueEquivalentVertex);

  // Calculate true direction (sum of daughter directions, normalized)
  if (trueParticle1 && trueParticle2) {
    TVector3 trueDir1(trueParticle1->Direction[0], trueParticle1->Direction[1], trueParticle1->Direction[2]);
    TVector3 trueDir2(trueParticle2->Direction[0], trueParticle2->Direction[1], trueParticle2->Direction[2]);
    TVector3 trueDirSum = trueDir1 + trueDir2;
    if(trueDirSum.Mag() > 0) trueDirSum = trueDirSum.Unit();
    trueEquivalentVertex->Direction[0] = trueDirSum.X();
    trueEquivalentVertex->Direction[1] = trueDirSum.Y();
    trueEquivalentVertex->Direction[2] = trueDirSum.Z();
  } else {
    trueEquivalentVertex->Direction[0] = -999.0;
    trueEquivalentVertex->Direction[1] = -999.0;
    trueEquivalentVertex->Direction[2] = -999.0;
  }

  // Copy fitted direction from reconstructed vertex
  trueEquivalentVertex->DirectionFit[0] = vertex->DirectionFit[0];
  trueEquivalentVertex->DirectionFit[1] = vertex->DirectionFit[1];
  trueEquivalentVertex->DirectionFit[2] = vertex->DirectionFit[2];

  // Initialize Pandora position and degeneracy
  for (int i = 0; i < 3; i++) {
    trueEquivalentVertex->PositionPandora[i] = -999;
  }
  trueEquivalentVertex->DegeneracyBeforeScoring = 0;
  trueEquivalentVertex->DegeneracyAfterScoring = 0;
  trueEquivalentVertex->NRecoParticles = 0;

  return trueEquivalentVertex;
}

//***************************************************************
std::vector<AnaVertexPD*> pdNeutralUtils::CreateVerticesCommon(
    AnaEventB& event,
    double maxVertexRadius,
    double maxDaughterDistance,
    double (*positionFinder)(AnaVertexPD*, double)) {
//***************************************************************

  // Note: maxVertexRadius parameter is not currently used in this implementation
  (void)maxVertexRadius; // Suppress unused parameter warning

  // Get the array of particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // OPTIMIZATION: Build hash map for O(1) particle lookups by UniqueID
  std::unordered_map<Int_t, AnaParticlePD*> particleByUniqueID;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part){
      particleByUniqueID[part->UniqueID] = part;
    }
  }

  // Create reconstructed vertices
  std::vector<AnaVertexPD*> reconstructedVertices;
  int vertexID = 0; // Counter for unique vertex IDs

  // Set to track which particle pairs have already been used to create vertices
  std::set<std::pair<AnaParticlePD*, AnaParticlePD*>> usedPairs;

  // Get track fit length from parameters
  double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

  // Get PID chi2/ndf thresholds from parameters
  double maxChi2NdfProton = ND::params().GetParameterD("neutralKaonAnalysis.MaxChi2NdfProton");
  double maxChi2NdfKaon = ND::params().GetParameterD("neutralKaonAnalysis.MaxChi2NdfKaon");

  // Get beam particle direction for angle calculations
  TVector3 beamDir(0, 0, 1); // Default to z-axis if beam not available
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(event.Beam);
  if (beam && beam->BeamParticle) {
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    if (beamPart) {
      // Use beam particle direction (preferably end direction at TPC entrance)
      if (beamPart->DirectionEnd[0] != -999 && beamPart->DirectionEnd[1] != -999 && beamPart->DirectionEnd[2] != -999) {
        beamDir.SetXYZ(beamPart->DirectionEnd[0], beamPart->DirectionEnd[1], beamPart->DirectionEnd[2]);
      } else if (beamPart->DirectionStart[0] != -999 && beamPart->DirectionStart[1] != -999 && beamPart->DirectionStart[2] != -999) {
        beamDir.SetXYZ(beamPart->DirectionStart[0], beamPart->DirectionStart[1], beamPart->DirectionStart[2]);
      }
      beamDir = beamDir.Unit(); // Ensure normalized
    }
  }

  for(int i = 0; i < nParts; i++){
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if(!daughter1) continue;

    // Check if daughter1 has valid start position
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) {
      continue; // Skip daughter1 with invalid start positions
    }
    // Check if daughter1 has valid end position
    if (daughter1->PositionEnd[0] < -900 || daughter1->PositionEnd[1] < -900 || daughter1->PositionEnd[2] < -900) {
      continue; // Skip daughter1 with invalid end positions
    }

    // Check PID compatibility for daughter1 - skip if compatible with proton or kaon
    std::pair<double, int> chi2Proton1 = pdAnaUtils::Chi2PID(*daughter1, 2212);  // Proton PDG
    std::pair<double, int> chi2Kaon1 = pdAnaUtils::Chi2PID(*daughter1, 321);    // Kaon PDG
    double chi2ndfProton1 = (chi2Proton1.second > 0) ? chi2Proton1.first / chi2Proton1.second : 9999.0;
    double chi2ndfKaon1 = (chi2Kaon1.second > 0) ? chi2Kaon1.first / chi2Kaon1.second : 9999.0;

    if (chi2ndfProton1 < maxChi2NdfProton || chi2ndfKaon1 < maxChi2NdfKaon) {
      continue; // Skip particles compatible with proton or kaon
    }

    for(int j = i + 1; j < nParts; j++){
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if(!daughter2) continue;

      // Skip if daughter1 and daughter2 are the same particle
      if (daughter1 == daughter2) {
        continue;
      }

      // Check if both particles have the same ParentID (EARLY CHECK - avoid useless computation)
      if (daughter1->ParentID != daughter2->ParentID) {
        continue; // Skip if parents are different
      }

      // Check if daughter2 has valid start position
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) {
        continue; // Skip daughter2 with invalid start positions
      }
      // Check if daughter2 has valid end position
      if (daughter2->PositionEnd[0] < -900 || daughter2->PositionEnd[1] < -900 || daughter2->PositionEnd[2] < -900) {
        continue; // Skip daughter2 with invalid end positions
      }

      // Check PID compatibility for daughter2 - skip if compatible with proton or kaon
      std::pair<double, int> chi2Proton2 = pdAnaUtils::Chi2PID(*daughter2, 2212);  // Proton PDG
      std::pair<double, int> chi2Kaon2 = pdAnaUtils::Chi2PID(*daughter2, 321);    // Kaon PDG
      double chi2ndfProton2 = (chi2Proton2.second > 0) ? chi2Proton2.first / chi2Proton2.second : 9999.0;
      double chi2ndfKaon2 = (chi2Kaon2.second > 0) ? chi2Kaon2.first / chi2Kaon2.second : 9999.0;

      if (chi2ndfProton2 < maxChi2NdfProton || chi2ndfKaon2 < maxChi2NdfKaon) {
        continue; // Skip particles compatible with proton or kaon
      }

      // Create a pair to check for duplicates (order-independent)
      std::pair<AnaParticlePD*, AnaParticlePD*> currentPair;
      if (daughter1 < daughter2) {
        currentPair = std::make_pair(daughter1, daughter2);
      } else {
        currentPair = std::make_pair(daughter2, daughter1);
      }

      // Check if this pair has already been used to create a vertex
      if (usedPairs.find(currentPair) != usedPairs.end()) {
        continue; // Skip if this pair has already been used
      }

      // Check if daughter1 and daughter2 are close enough
      TVector3 pos1 = pdAnaUtils::DefinePosition(daughter1);
      TVector3 pos2 = pdAnaUtils::DefinePosition(daughter2);

      // Check if positions are valid
      if (pos1.X() < -900 || pos2.X() < -900) {
        continue; // Skip if positions are invalid
      }

      // Convert TVector3 to arrays for GetSeparationSquared
      Float_t pos1_array[3] = {static_cast<Float_t>(pos1.X()), static_cast<Float_t>(pos1.Y()), static_cast<Float_t>(pos1.Z())};
      Float_t pos2_array[3] = {static_cast<Float_t>(pos2.X()), static_cast<Float_t>(pos2.Y()), static_cast<Float_t>(pos2.Z())};

      float distance = sqrt(anaUtils::GetSeparationSquared(pos1_array, pos2_array));

      if(distance > maxDaughterDistance){
        continue; // Skip daughter1 and daughter2 if they are not close enough
      }

      // Create reconstructed vertex
      AnaVertexPD* reconstructedVertex = new AnaVertexPD();
      reconstructedVertex->OriginalDistance = distance;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->NParticles = 2;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);
      reconstructedVertex->ParentID = daughter1->ParentID;

      // Initialize vertex position and minimum distance to invalid values
      reconstructedVertex->Position[0] = -999.0;
      reconstructedVertex->Position[1] = -999.0;
      reconstructedVertex->Position[2] = -999.0;
      reconstructedVertex->MinimumDistance = -999.0;
      reconstructedVertex->Score = -999.0;

      // Call the position finder function to calculate vertex position and score
      positionFinder(reconstructedVertex, trackFitLength);

      // Calculate vertex direction (average of daughter directions)
      Float_t direction[3] = {daughter1->DirectionStart[0] + daughter2->DirectionStart[0],
                              daughter1->DirectionStart[1] + daughter2->DirectionStart[1],
                              daughter1->DirectionStart[2] + daughter2->DirectionStart[2]};

      Float_t norm = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
      direction[0] = direction[0] / norm;
      direction[1] = direction[1] / norm;
      direction[2] = direction[2] / norm;

      reconstructedVertex->Direction[0] = direction[0];
      reconstructedVertex->Direction[1] = direction[1];
      reconstructedVertex->Direction[2] = direction[2];

      // Compute the angle between the start directions of the two particles in degrees
      Float_t cosAngle = daughter1->DirectionStart[0]*daughter2->DirectionStart[0] +
                         daughter1->DirectionStart[1]*daughter2->DirectionStart[1] +
                         daughter1->DirectionStart[2]*daughter2->DirectionStart[2];
      // Clamp cosAngle to [-1, 1] to avoid numerical errors in acos
      if (cosAngle > 1.0) cosAngle = 1.0;
      if (cosAngle < -1.0) cosAngle = -1.0;
      reconstructedVertex->OpeningAngle = TMath::ACos(cosAngle) * 180.0 / TMath::Pi();

      // Calculate vertex momentum, energy, and angle with beam
      // Ensure particles have momentum calculated first
      reconstructedVertex->EnsureParticleMomentum();

      const double pionMass = 0.13957; // GeV
      TVector3 p1vec(daughter1->DirectionStart[0] * daughter1->Momentum,
                    daughter1->DirectionStart[1] * daughter1->Momentum,
                    daughter1->DirectionStart[2] * daughter1->Momentum);
      TVector3 p2vec(daughter2->DirectionStart[0] * daughter2->Momentum,
                    daughter2->DirectionStart[1] * daughter2->Momentum,
                    daughter2->DirectionStart[2] * daughter2->Momentum);
      TVector3 pTot = p1vec + p2vec;

      reconstructedVertex->Momentum[0] = pTot.X();
      reconstructedVertex->Momentum[1] = pTot.Y();
      reconstructedVertex->Momentum[2] = pTot.Z();

      double E1 = sqrt(daughter1->Momentum*daughter1->Momentum + pionMass*pionMass);
      double E2 = sqrt(daughter2->Momentum*daughter2->Momentum + pionMass*pionMass);
      reconstructedVertex->Energy = (E1 + E2) * 1000; // Convert to MeV

      // Calculate angle with beam direction
      TVector3 vtxDir(reconstructedVertex->Direction[0],
                     reconstructedVertex->Direction[1],
                     reconstructedVertex->Direction[2]);
      reconstructedVertex->AngleWithBeam = vtxDir.Angle(beamDir);

      // Create and fill true equivalent vertex using helper function
      reconstructedVertex->TrueEquivalentVertex = FillTrueEquivalentVertex(reconstructedVertex);

      // Set the generation and process
      reconstructedVertex->Generation = -999;
      reconstructedVertex->Process = -999;

      // Validate vertex before adding
      if (ValidateVertex(reconstructedVertex)) {
        reconstructedVertices.push_back(reconstructedVertex);
        // Mark this pair as used to prevent duplicate vertices
        usedPairs.insert(currentPair);
      } else {
        // Delete the vertex if validation failed
        delete reconstructedVertex;
      }
    }
  }

  // Calculate degeneracy before scoring
  double degThreshold = ND::params().GetParameterD("neutralKaonAnalysis.DegeneracyDistance");
  bool usePandora = (bool)ND::params().GetParameterI("neutralKaonAnalysis.UsePandoraForDegeneracy");
  CalculateVertexDegeneracy(reconstructedVertices, usePandora, degThreshold);

  // Calculate isolation distances for all vertices
  CalculateVertexIsolation(reconstructedVertices, degThreshold);

  // Filter vertices to ensure each particle belongs to at most one vertex
  std::vector<AnaVertexPD*> filteredVertices = FilterVerticesByScore(reconstructedVertices);

  // Calculate degeneracy after scoring
  for (auto* vtx : filteredVertices) {
    if (!vtx) continue;

    TVector3 pos1;
    if (usePandora) {
      pos1.SetXYZ(vtx->PositionPandora[0], vtx->PositionPandora[1], vtx->PositionPandora[2]);
    } else {
      pos1.SetXYZ(vtx->Position[0], vtx->Position[1], vtx->Position[2]);
    }

    if (pos1.X() < -900) continue;

    int degeneracy = 0;
    for (auto* vtx2 : filteredVertices) {
      if (vtx == vtx2 || !vtx2) continue;

      TVector3 pos2;
      if (usePandora) {
        pos2.SetXYZ(vtx2->PositionPandora[0], vtx2->PositionPandora[1], vtx2->PositionPandora[2]);
      } else {
        pos2.SetXYZ(vtx2->Position[0], vtx2->Position[1], vtx2->Position[2]);
      }

      if (pos2.X() < -900) continue;

      if ((pos1 - pos2).Mag() < degThreshold) {
        degeneracy++;
      }
    }

    vtx->DegeneracyAfterScoring = degeneracy;
  }

  return filteredVertices;
}

//***************************************************************
std::pair<Int_t, Int_t> pdNeutralUtils::CalculateNeutralScore(
    AnaNeutralParticlePD* neutralParticle,
    AnaVertexPD* vertex,
    AnaParticlePD* parentParticle,
    AnaEventB& event,
    const std::unordered_map<Int_t, AnaParticlePD*>& particleByUniqueID) {
//***************************************************************

  // Initialize return values
  Int_t NPotentialParents = 0;
  Int_t NRecoHitsInVertex = 0;

  if (!neutralParticle || !vertex || !parentParticle) {
    return std::make_pair(NPotentialParents, NRecoHitsInVertex);
  }

  // Get array of particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // Get parameters
  double daughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
  double cylinderRadius = ND::params().GetParameterD("neutralKaonAnalysis.CylinderRadius");

  TVector3 vertexPos(vertex->Position[0], vertex->Position[1], vertex->Position[2]);

  // === ENHANCED NEUTRAL PARTICLE SCORING ALGORITHM ===
  // LOWER score = more neutral-like (few hits, far from path, not aligned)
  // HIGHER score = more charged-like (many hits, close to path, aligned)

  // Define neutral particle trajectory
  TVector3 neutralStart(neutralParticle->PositionStart[0],
                        neutralParticle->PositionStart[1],
                        neutralParticle->PositionStart[2]);
  TVector3 neutralEnd(neutralParticle->PositionEnd[0],
                      neutralParticle->PositionEnd[1],
                      neutralParticle->PositionEnd[2]);
  TVector3 neutralDirection = (neutralEnd - neutralStart).Unit();
  double neutralLength = neutralParticle->Length;

  // Collect all hits within the cylinder for analysis
  std::vector<TVector3> hitsInCylinder;
  std::vector<double> perpDistances;
  std::vector<double> longitudinalProjections;
  double totalDistance = 0.0;
  int nHitsInCylinder = 0;

  for(int j = 0; j < nParts; j++){
    AnaParticlePD* potentialParent = static_cast<AnaParticlePD*>(parts[j]);
    if(!potentialParent) continue;

    // Skip vertex particles and current parent
    if(vertex->Particles.size() >= 2){
      if(potentialParent->UniqueID == vertex->Particles[0]->UniqueID ||
         potentialParent->UniqueID == vertex->Particles[1]->UniqueID ||
         potentialParent->UniqueID == parentParticle->UniqueID) continue;
    }

    TVector3 potentialParentEnd(potentialParent->PositionEnd[0],
                                 potentialParent->PositionEnd[1],
                                 potentialParent->PositionEnd[2]);
    double potentialParentDistance = (potentialParentEnd - vertexPos).Mag();

    if(potentialParentDistance < daughterDistance){
      // OPTIMIZATION: Count potential parents (merged from CreateNeutral loop)
      NPotentialParents++;

      // Check hits from this particle
      for(int plane = 0; plane < 3; plane++){
        for(size_t k = 0; k < potentialParent->Hits[plane].size(); k++){
          AnaHitPD& hit = potentialParent->Hits[plane][k];
          if(hit.Position.Z() == -999) continue;

          TVector3 hitPos = hit.Position;

          // Check if hit is within cylinder
          TVector3 startToHit = hitPos - neutralStart;
          double projection = startToHit.Dot(neutralDirection);

          // Only consider hits along the neutral particle path
          if(projection >= 0 && projection <= neutralLength){
            TVector3 projectionPoint = neutralStart + projection * neutralDirection;
            double perpDistance = (hitPos - projectionPoint).Mag();

            if(perpDistance < cylinderRadius){
              nHitsInCylinder++;
              hitsInCylinder.push_back(hitPos);
              perpDistances.push_back(perpDistance);
              longitudinalProjections.push_back(projection);
              totalDistance += perpDistance;
              // OPTIMIZATION: Count hits in vertex (merged from CreateNeutral loop)
              NRecoHitsInVertex++;
            }
          }
        }
      }
    }
  }

  // === 1. ENHANCED DISTANCE METRICS ===

  // Calculate average perpendicular distance
  double avgDistance = (nHitsInCylinder > 0) ? totalDistance / nHitsInCylinder : 0.0;

  // Calculate RMS of perpendicular distances
  double rmsDistance = 0.0;
  if(nHitsInCylinder > 1){
    double sumSquaredDiff = 0.0;
    for(double dist : perpDistances){
      sumSquaredDiff += (dist - avgDistance) * (dist - avgDistance);
    }
    rmsDistance = TMath::Sqrt(sumSquaredDiff / nHitsInCylinder);
  }

  // === 2. LONGITUDINAL SPAN CALCULATION ===

  double longitudinalSpan = 0.0;
  if(nHitsInCylinder > 0 && neutralLength > 0){
    double minProj = *std::min_element(longitudinalProjections.begin(), longitudinalProjections.end());
    double maxProj = *std::max_element(longitudinalProjections.begin(), longitudinalProjections.end());
    double span = maxProj - minProj;
    longitudinalSpan = span / neutralLength;  // Fraction of path covered by hits
  }

  // === 3. ROBUST ALIGNMENT CALCULATION ===

  double hitsAlignment = 0.0;  // Default: no alignment

  if(nHitsInCylinder == 1){
    // Single hit: check if vector from start to hit is aligned with neutral direction
    TVector3 hitVector = (hitsInCylinder[0] - neutralStart);
    if(hitVector.Mag() > 0){
      hitVector = hitVector.Unit();
      hitsAlignment = TMath::Abs(hitVector.Dot(neutralDirection));
    }
  }
  else if(nHitsInCylinder == 2){
    // Two hits: fit line between them
    TVector3 hitLine = (hitsInCylinder[1] - hitsInCylinder[0]);
    if(hitLine.Mag() > 0){
      hitLine = hitLine.Unit();
      hitsAlignment = TMath::Abs(hitLine.Dot(neutralDirection));
    }
  }
  else if(nHitsInCylinder >= 3){
    // Three or more hits: use PCA
    // Calculate centroid of hits
    TVector3 centroid(0, 0, 0);
    for(const auto& hit : hitsInCylinder){
      centroid += hit;
    }
    centroid *= (1.0 / hitsInCylinder.size());

    // Build covariance matrix
    double cov_xx = 0, cov_yy = 0, cov_zz = 0;
    double cov_xy = 0, cov_xz = 0, cov_yz = 0;

    for(const auto& hit : hitsInCylinder){
      TVector3 diff = hit - centroid;
      cov_xx += diff.X() * diff.X();
      cov_yy += diff.Y() * diff.Y();
      cov_zz += diff.Z() * diff.Z();
      cov_xy += diff.X() * diff.Y();
      cov_xz += diff.X() * diff.Z();
      cov_yz += diff.Y() * diff.Z();
    }

    // Normalize
    double norm = 1.0 / hitsInCylinder.size();
    cov_xx *= norm;
    cov_yy *= norm;
    cov_zz *= norm;
    cov_xy *= norm;
    cov_xz *= norm;
    cov_yz *= norm;

    // Find principal eigenvector (direction of maximum variance)
    // Use power iteration method
    TVector3 fittedDirection(1, 1, 1);  // Initial guess
    for(int iter = 0; iter < 20; iter++){
      TVector3 newDir;
      newDir.SetX(cov_xx * fittedDirection.X() + cov_xy * fittedDirection.Y() + cov_xz * fittedDirection.Z());
      newDir.SetY(cov_xy * fittedDirection.X() + cov_yy * fittedDirection.Y() + cov_yz * fittedDirection.Z());
      newDir.SetZ(cov_xz * fittedDirection.X() + cov_yz * fittedDirection.Y() + cov_zz * fittedDirection.Z());

      double mag = newDir.Mag();
      if(mag > 0){
        fittedDirection = newDir * (1.0 / mag);
      }
    }

    // Calculate alignment as absolute value of dot product
    hitsAlignment = TMath::Abs(fittedDirection.Dot(neutralDirection));
  }

  // === 4. CALCULATE SCORE COMPONENTS ===

  // A) Hit density (hits per unit distance)
  // This will be a multiplicative factor for hit-dependent scores
  // Ensures particles with 0 hits get 0 score for hit-dependent component
  double hitDensity = 0.0;
  if(neutralLength > 0){
    hitDensity = nHitsInCylinder / neutralLength;
  }

  // B) Enhanced distance score - combines average and RMS
  // Scattered hits (high RMS) are more neutral-like
  // Concentrated hits (low RMS) are more charged-like
  double distanceScore = 0.0;
  if(nHitsInCylinder > 0 && cylinderRadius > 0){
    double avgFactor = 1.0 - TMath::Min(avgDistance / cylinderRadius, 1.0);
    double rmsFactor = (nHitsInCylinder > 1) ? (1.0 - TMath::Min(rmsDistance / cylinderRadius, 1.0)) : avgFactor;
    // Lower RMS (concentrated) → higher rmsFactor → higher score (charged-like)
    distanceScore = avgFactor * rmsFactor * 100.0;
  }

  // C) Alignment score - already calculated (0 to 1, higher = more aligned)
  double alignmentScore = hitsAlignment * 100.0;

  // D) Clustering score - longitudinal span
  // High span (hits spread along path) → charged-like
  // Low span (localized hits) → could be neutral or random
  double clusterScore = longitudinalSpan * 100.0;

  // === 5. CALCULATE SIBLING ALIGNMENT SCORE ===
  // Look at alignment between all parent's daughters (including vertex particles) and neutral direction

  double siblingAlignmentScore = 0.0;
  int nSiblings = 0;
  double totalSiblingAlignment = 0.0;

  if(parentParticle && parentParticle->Daughters.size() > 0){
    // Loop through parent's daughters (including vertex particles)
    for(size_t i = 0; i < parentParticle->Daughters.size(); i++){
      AnaParticlePD* sibling = static_cast<AnaParticlePD*>(parentParticle->Daughters[i]);
      if(!sibling) continue;

      // Check if sibling has valid direction
      if(sibling->DirectionStart[0] < -900) continue;

      // Calculate alignment between sibling direction and neutral direction
      TVector3 siblingDir(sibling->DirectionStart[0],
                         sibling->DirectionStart[1],
                         sibling->DirectionStart[2]);
      if(siblingDir.Mag() > 0){
        siblingDir = siblingDir.Unit();
        double alignment = TMath::Abs(siblingDir.Dot(neutralDirection));
        totalSiblingAlignment += alignment;
        nSiblings++;
      }
    }

    // Average sibling alignment
    if(nSiblings > 0){
      siblingAlignmentScore = (totalSiblingAlignment / nSiblings) * 100.0;
    }
  }

  // === 6. COMBINED SCORING ===
  // neutralScore = (nhits/length) * (alignmentScore + distanceScore + clusterScore) + siblingAlignmentScore
  // Equal weights for hit-dependent components, multiplied by hit density
  // Sibling alignment is independent of hit density

  double hitDependentScore = hitDensity * (alignmentScore + distanceScore + clusterScore);
  double neutralScore = hitDependentScore + siblingAlignmentScore;

  // Store all metrics in neutral particle
  neutralParticle->NeutralScore = neutralScore;
  neutralParticle->HitsAlignment = hitsAlignment;
  neutralParticle->NHitsInCylinder = nHitsInCylinder;
  neutralParticle->HitsAvgDistance = avgDistance;
  neutralParticle->HitsRMSDistance = rmsDistance;
  neutralParticle->HitsLongitudinalSpan = longitudinalSpan;

  // OPTIMIZATION: Return counts that were previously calculated in a separate loop
  return std::make_pair(NPotentialParents, NRecoHitsInVertex);
}

//***************************************************************
AnaTrueEquivalentNeutralParticlePD* pdNeutralUtils::FillTrueEquivalentNeutralParticle(
    AnaVertexPD* vertex,
    AnaParticlePD* parentParticle) {
//***************************************************************

  if (!vertex || !parentParticle) {
    return nullptr;
  }

  // Create true equivalent neutral particle
  AnaTrueEquivalentNeutralParticlePD* trueEquivalentNeutralParticle = new AnaTrueEquivalentNeutralParticlePD();
  trueEquivalentNeutralParticle->TrueEquivalentVertex = static_cast<AnaTrueEquivalentVertexPD*>(vertex->TrueEquivalentVertex);
  trueEquivalentNeutralParticle->TrueParent = static_cast<AnaTrueParticlePD*>(parentParticle->TrueObject);

  // Add the position and direction of the true particle
  if (parentParticle->TrueObject) {
    AnaTrueParticlePD* trueParticle = static_cast<AnaTrueParticlePD*>(parentParticle->TrueObject);
    trueEquivalentNeutralParticle->Position[0] = trueParticle->PositionEnd[0];
    trueEquivalentNeutralParticle->Position[1] = trueParticle->PositionEnd[1];
    trueEquivalentNeutralParticle->Position[2] = trueParticle->PositionEnd[2];
  } else {
    trueEquivalentNeutralParticle->Position[0] = -999.0;
    trueEquivalentNeutralParticle->Position[1] = -999.0;
    trueEquivalentNeutralParticle->Position[2] = -999.0;
  }

  // Add the end position of the true vertex
  if (vertex->TrueEquivalentVertex) {
    trueEquivalentNeutralParticle->PositionEnd[0] = vertex->TrueEquivalentVertex->Position[0];
    trueEquivalentNeutralParticle->PositionEnd[1] = vertex->TrueEquivalentVertex->Position[1];
    trueEquivalentNeutralParticle->PositionEnd[2] = vertex->TrueEquivalentVertex->Position[2];
  } else {
    trueEquivalentNeutralParticle->PositionEnd[0] = -999.0;
    trueEquivalentNeutralParticle->PositionEnd[1] = -999.0;
    trueEquivalentNeutralParticle->PositionEnd[2] = -999.0;
  }

  // Calculate the direction of the true particle
  Float_t trueDirectionStart[3];
  trueDirectionStart[0] = trueEquivalentNeutralParticle->PositionEnd[0] - trueEquivalentNeutralParticle->Position[0];
  trueDirectionStart[1] = trueEquivalentNeutralParticle->PositionEnd[1] - trueEquivalentNeutralParticle->Position[1];
  trueDirectionStart[2] = trueEquivalentNeutralParticle->PositionEnd[2] - trueEquivalentNeutralParticle->Position[2];

  Float_t trueNorm = sqrt(trueDirectionStart[0]*trueDirectionStart[0] +
                          trueDirectionStart[1]*trueDirectionStart[1] +
                          trueDirectionStart[2]*trueDirectionStart[2]);
  if (trueNorm > 0) {
    trueDirectionStart[0] /= trueNorm;
    trueDirectionStart[1] /= trueNorm;
    trueDirectionStart[2] /= trueNorm;
  } else {
    trueDirectionStart[0] = -999.0;
    trueDirectionStart[1] = -999.0;
    trueDirectionStart[2] = -999.0;
  }

  trueEquivalentNeutralParticle->Direction[0] = trueDirectionStart[0];
  trueEquivalentNeutralParticle->Direction[1] = trueDirectionStart[1];
  trueEquivalentNeutralParticle->Direction[2] = trueDirectionStart[2];

  Float_t trueDirectionEnd[3];
  trueDirectionEnd[0] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[0] : -999;
  trueDirectionEnd[1] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[1] : -999;
  trueDirectionEnd[2] = vertex->TrueEquivalentVertex ? vertex->TrueEquivalentVertex->Direction[2] : -999;

  Float_t trueNormEnd = sqrt(trueDirectionEnd[0]*trueDirectionEnd[0] +
                             trueDirectionEnd[1]*trueDirectionEnd[1] +
                             trueDirectionEnd[2]*trueDirectionEnd[2]);
  if (trueNormEnd > 0) {
    trueDirectionEnd[0] /= trueNormEnd;
    trueDirectionEnd[1] /= trueNormEnd;
    trueDirectionEnd[2] /= trueNormEnd;
  } else {
    trueDirectionEnd[0] = -999.0;
    trueDirectionEnd[1] = -999.0;
    trueDirectionEnd[2] = -999.0;
  }

  trueEquivalentNeutralParticle->DirectionEnd[0] = trueDirectionEnd[0];
  trueEquivalentNeutralParticle->DirectionEnd[1] = trueDirectionEnd[1];
  trueEquivalentNeutralParticle->DirectionEnd[2] = trueDirectionEnd[2];

  Float_t trueLength = sqrt(pow(trueEquivalentNeutralParticle->PositionEnd[0]-trueEquivalentNeutralParticle->Position[0],2)+
                            pow(trueEquivalentNeutralParticle->PositionEnd[1]-trueEquivalentNeutralParticle->Position[1],2)+
                            pow(trueEquivalentNeutralParticle->PositionEnd[2]-trueEquivalentNeutralParticle->Position[2],2));
  trueEquivalentNeutralParticle->Length = trueLength;

  // Calculate invariant mass and momentum for true equivalent neutral particle
  Float_t trueInvariantMass = -999;
  Float_t trueMomentumEnd = -999;

  if (vertex->TrueEquivalentVertex && vertex->TrueEquivalentVertex->TrueParticles.size() >= 2) {
    AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[0]);
    AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[1]);

    if (trueParticle1 && trueParticle2 &&
        trueParticle1->Momentum > 0 && trueParticle2->Momentum > 0 &&
        trueParticle1->Momentum != -999 && trueParticle2->Momentum != -999) {
      const Float_t pionMass = 0.13957;
      trueInvariantMass = pdAnaUtils::ComputeTrueInvariantMass(*trueParticle1, *trueParticle2, pionMass, pionMass);

      // Calculate total momentum (sum of daughter momenta)
      trueMomentumEnd = trueParticle1->Momentum + trueParticle2->Momentum;
    }
  }

  trueEquivalentNeutralParticle->Mass = trueInvariantMass;
  trueEquivalentNeutralParticle->MomentumEnd = trueMomentumEnd;

  return trueEquivalentNeutralParticle;
}

//***************************************************************
AnaNeutralParticlePD* pdNeutralUtils::CreateNeutral(
    AnaEventB& event,
    AnaVertexPD* vertex,
    int neutralParticleID,
    const std::unordered_map<Int_t, AnaParticlePD*>& particleByUniqueID){
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return nullptr;
  }

  // Find the parent particle using the vertex's ParentID
  AnaParticlePD* parentParticle = nullptr;

  // OPTIMIZATION: O(1) hash map lookup instead of O(n) linear search
  auto it = particleByUniqueID.find(vertex->ParentID);
  if(it != particleByUniqueID.end()){
    parentParticle = it->second;
  }

  // If not found in regular particles, check beam particle
  if(!parentParticle){
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(event.Beam);
    if(beam && beam->BeamParticle){
      AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
      if(beamPart && beamPart->UniqueID == vertex->ParentID){
        parentParticle = beamPart;
      }
    }
  }

  // If parent not found, return nullptr
  if(!parentParticle){
    return nullptr;
  }

  // Create neutral particle
  AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();

  // Set vertex and parent
  neutralParticle->Vertex = vertex;
  neutralParticle->Parent = parentParticle;
  neutralParticle->UniqueID = neutralParticleID;

  // Calculate neutral particle start position
  bool useRecalculation = (bool)ND::params().GetParameterI("neutralKaonAnalysis.UseNeutralStartRecalculation");

  if (useRecalculation) {
    // Use improved start position calculation
    double proximityThreshold = ND::params().GetParameterD("neutralKaonAnalysis.NeutralStartProximityThreshold");
    TVector3 neutralStart = FindNeutralParticleStartPosition(neutralParticle, parentParticle, event, proximityThreshold);

    neutralParticle->PositionStart[0] = neutralStart.X();
    neutralParticle->PositionStart[1] = neutralStart.Y();
    neutralParticle->PositionStart[2] = neutralStart.Z();
    neutralParticle->PositionStart[3] = parentParticle->PositionEnd[3];
  } else {
    // Use parent end position directly
    neutralParticle->PositionStart[0] = parentParticle->PositionEnd[0];
    neutralParticle->PositionStart[1] = parentParticle->PositionEnd[1];
    neutralParticle->PositionStart[2] = parentParticle->PositionEnd[2];
    neutralParticle->PositionStart[3] = parentParticle->PositionEnd[3];
  }

  neutralParticle->PositionEnd[0] = vertex->Position[0];
  neutralParticle->PositionEnd[1] = vertex->Position[1];
  neutralParticle->PositionEnd[2] = vertex->Position[2];
  neutralParticle->PositionEnd[3] = -999; // Time component

  // Calculate start direction as the vector from PositionStart to PositionEnd
  Float_t directionStart[3];
  directionStart[0] = neutralParticle->PositionEnd[0] - neutralParticle->PositionStart[0];
  directionStart[1] = neutralParticle->PositionEnd[1] - neutralParticle->PositionStart[1];
  directionStart[2] = neutralParticle->PositionEnd[2] - neutralParticle->PositionStart[2];

  // Normalize the direction vector
  Float_t norm = sqrt(directionStart[0]*directionStart[0] +
                      directionStart[1]*directionStart[1] +
                      directionStart[2]*directionStart[2]);
  if (norm > 0) {
    directionStart[0] /= norm;
    directionStart[1] /= norm;
    directionStart[2] /= norm;
  }

  neutralParticle->DirectionStart[0] = directionStart[0];
  neutralParticle->DirectionStart[1] = directionStart[1];
  neutralParticle->DirectionStart[2] = directionStart[2];

  // Set the end direction from the vertex direction
  neutralParticle->DirectionEnd[0] = vertex->Direction[0];
  neutralParticle->DirectionEnd[1] = vertex->Direction[1];
  neutralParticle->DirectionEnd[2] = vertex->Direction[2];

  // Calculate length
  neutralParticle->Length = sqrt(pow(neutralParticle->PositionEnd[0]-neutralParticle->PositionStart[0],2)+
                                 pow(neutralParticle->PositionEnd[1]-neutralParticle->PositionStart[1],2)+
                                 pow(neutralParticle->PositionEnd[2]-neutralParticle->PositionStart[2],2));

  // Calculate impact parameter and create FitParent from parent track extrapolation
  neutralParticle->ImpactParameter = -999;
  if(parentParticle){
    // Get track fit length from parameters
    double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

    // Extrapolate parent track from end position
    std::vector<double> parentLineParams;
    pdAnaUtils::ExtrapolateTrack(parentParticle, parentLineParams, trackFitLength, false); // Use end position for parent

    // Calculate impact parameter if extrapolation was successful
    if(parentLineParams.size() == 6 && parentLineParams[0] != -999.0){
      TVector3 vertexPos(vertex->Position[0], vertex->Position[1], vertex->Position[2]);
      neutralParticle->ImpactParameter = pdAnaUtils::CalculateImpactParameter(parentLineParams, vertexPos);
    }
  }

  // Calculate creation vertex particle quantities: proton count, total count, chi2, distances, PDGs
  double protonCreationVtxDist = ND::params().GetParameterD("neutralKaonAnalysis.ProtonCreationVtxDistance");
  CalculateNProtonInCreationVtx(neutralParticle, parentParticle, protonCreationVtxDist);

  // OPTIMIZATION: Calculate neutral particle score and get counts in a single pass
  // This replaces the duplicate loop that was previously here (lines 1903-1973)
  auto [NPotentialParents, NRecoHitsInVertex] = CalculateNeutralScore(
      neutralParticle, vertex, parentParticle, event, particleByUniqueID);

  neutralParticle->Vertex->NPotentialParents = NPotentialParents;
  neutralParticle->NRecoHitsInVertex = NRecoHitsInVertex;

  // Create and fill true equivalent neutral particle using helper function
  neutralParticle->TrueEquivalentNeutralParticle = FillTrueEquivalentNeutralParticle(vertex, parentParticle);

  // Logic to assign TrueObject to the neutral particle
  if(vertex->TrueEquivalentVertex && vertex->TrueEquivalentVertex->TrueParticles.size() >= 2){
    AnaTrueParticlePD* trueDaughter1 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[0]);
    AnaTrueParticlePD* trueDaughter2 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[1]);
    AnaTrueParticlePD* trueParent = static_cast<AnaTrueParticlePD*>(parentParticle->TrueObject);

    if(trueDaughter1 && trueDaughter2 && trueParent && trueParent->Daughters.size() > 0){
      // Both daughters have the same parent
      if(trueDaughter1->ParentID == trueDaughter2->ParentID){
        // Find if there is a true particle that is the parent of the two daughters
        for(size_t i = 0; i < trueParent->Daughters.size(); i++){
          Int_t daughterID = trueParent->Daughters[i];
          AnaTrueParticlePD* trueDaughter = pdAnaUtils::GetTrueParticle(&event, daughterID);
          if(trueDaughter && (trueDaughter->ID == trueDaughter1->ParentID || trueDaughter->ID == trueDaughter2->ParentID)){
            // If found, assign it to the true object of the neutral particle
            neutralParticle->TrueObject = trueDaughter;
            // Check if the true object has a reconstructed particle associated to it
            if(trueDaughter->ReconParticles.size() > 0){
              neutralParticle->RecoParticle = static_cast<AnaParticlePD*>(trueDaughter->ReconParticles[0]);
            }
            break;
          }
        }
      }
      else{
        neutralParticle->TrueObject = nullptr;
        neutralParticle->RecoParticle = nullptr;
      }
    }
    else{
      neutralParticle->TrueObject = nullptr;
      neutralParticle->RecoParticle = nullptr;
    }
  }
  else{
    neutralParticle->TrueObject = nullptr;
    neutralParticle->RecoParticle = nullptr;
  }

  // Ensure particles have reliable momentum values
  vertex->EnsureParticleMomentum();

  // Calculate invariant mass from the two particles in the vertex
  Float_t invariantMass = -999;
  if (vertex->Particles.size() >= 2) {
    AnaParticlePD* particle1 = vertex->Particles[0];
    AnaParticlePD* particle2 = vertex->Particles[1];

    if (particle1 && particle2 &&
        particle1->Momentum > 0 && particle2->Momentum > 0 &&
        particle1->Momentum != -999 && particle2->Momentum != -999) {
      const Float_t pionMass = 0.13957;
      invariantMass = anaUtils::ComputeInvariantMass(*particle1, *particle2, pionMass, pionMass);
    }
  }

  // Set other properties
  neutralParticle->Mass = invariantMass;
  neutralParticle->Momentum = -999;
  neutralParticle->PDG = -999;
  neutralParticle->Lifetime = -999;
  neutralParticle->DecayLength = -999;

  return neutralParticle;
}

//***************************************************************
std::vector<AnaNeutralParticlePD*> pdNeutralUtils::CreateNeutrals(AnaEventB& event, const std::vector<AnaVertexPD*>& vertices){
//***************************************************************

  std::vector<AnaNeutralParticlePD*> neutralParticles;

  int neutralParticleID = 0;

  // OPTIMIZATION: Build hash map once for all neutral particles
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  std::unordered_map<Int_t, AnaParticlePD*> particleByUniqueID;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part){
      particleByUniqueID[part->UniqueID] = part;
    }
  }

  // Loop over all vertices and create one neutral particle per vertex
  for(size_t v = 0; v < vertices.size(); v++){
    AnaVertexPD* vertex = vertices[v];
    if(!vertex) continue;

    AnaNeutralParticlePD* neutralParticle = CreateNeutral(event, vertex, neutralParticleID, particleByUniqueID);
    if(neutralParticle){
      neutralParticles.push_back(neutralParticle);
      neutralParticleID++;
    }
  }

  return neutralParticles;
}

//***************************************************************
std::vector<AnaVertexPD*> pdNeutralUtils::CreateVertices(AnaEventB& event, double maxVertexRadius, double maxDaughterDistance) {
//***************************************************************

  // Check parameter to decide which algorithm to use
  int algorithmChoice = ND::params().GetParameterI("neutralKaonAnalysis.VertexFindingMethod");

  // Select the appropriate position finder function
  double (*positionFinder)(AnaVertexPD*, double) = nullptr;

  // 0 = Geometric, 1 = Fitted (TMinuit), 2 = Kalman
  switch(algorithmChoice) {
    case 0:
      positionFinder = FindVertexPositionGeometric;
      break;
    case 1:
      positionFinder = FindVertexPositionWithFit;
      break;
    case 2:
      positionFinder = FindVertexPositionKalman;
      break;
    default:
      std::cout << "WARNING: Unknown vertex finding method " << algorithmChoice
                << ", using Fitted method (1)" << std::endl;
      positionFinder = FindVertexPositionWithFit;
      break;
  }

  // Call common vertex creation function with selected position finder
  return CreateVerticesCommon(event, maxVertexRadius, maxDaughterDistance, positionFinder);
}

