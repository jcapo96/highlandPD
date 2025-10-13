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
//***************************************************************
void CalculatePandoraVertexPosition(AnaParticlePD* particle1, AnaParticlePD* particle2, Float_t* position) {
    if (!particle1 || !particle2) {
        for (int i = 0; i < 3; i++) position[i] = -999.;
        return;
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
    if (fabs(denom) > 1e-6) {
        double s = (b*e - c*d) / denom;
        double t = (a*e - b*d) / denom;
        TVector3 p1 = pos1 + s * dir1;
        TVector3 p2 = pos2 + t * dir2;
        pandoraVertex = 0.5 * (p1 + p2);
    } else {
        pandoraVertex = 0.5 * (pos1 + pos2);
    }

    position[0] = pandoraVertex.X();
    position[1] = pandoraVertex.Y();
    position[2] = pandoraVertex.Z();
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

  // Store fitted vertex position
  vertex->Position[0] = vx;
  vertex->Position[1] = vy;
  vertex->Position[2] = vz;

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

  // Create fitted particle objects to store fit results
  AnaVertexFittedParticlePD* fitPart1 = new AnaVertexFittedParticlePD();
  fitPart1->FitPosition[0] = vx;
  fitPart1->FitPosition[1] = vy;
  fitPart1->FitPosition[2] = vz;
  fitPart1->FitDirection[0] = fittedDir1.X();
  fitPart1->FitDirection[1] = fittedDir1.Y();
  fitPart1->FitDirection[2] = fittedDir1.Z();

  AnaVertexFittedParticlePD* fitPart2 = new AnaVertexFittedParticlePD();
  fitPart2->FitPosition[0] = vx;
  fitPart2->FitPosition[1] = vy;
  fitPart2->FitPosition[2] = vz;
  fitPart2->FitDirection[0] = fittedDir2.X();
  fitPart2->FitDirection[1] = fittedDir2.Y();
  fitPart2->FitDirection[2] = fittedDir2.Z();

  vertex->FitParticles.push_back(fitPart1);
  vertex->FitParticles.push_back(fitPart2);

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = fittedDir1 + fittedDir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->FitDirection[0] = fitDirSum.X();
  vertex->FitDirection[1] = fitDirSum.Y();
  vertex->FitDirection[2] = fitDirSum.Z();

  // Calculate Pandora-based vertex position
  CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

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

  // Use existing geometric position finding
  pdAnaUtils::FindVertexPosition(vertex);

  // Set score to minimum distance (lower is better)
  vertex->Score = vertex->MinimumDistance;

  // Create fitted particle objects from the stored line parameters
  if(vertex->FittedLineParams.size() >= 2){
    // Extract line parameters for daughter 1
    std::vector<double> line1 = vertex->FittedLineParams[0];
    AnaVertexFittedParticlePD* fitPart1 = new AnaVertexFittedParticlePD();
    fitPart1->FitPosition[0] = line1[0];
    fitPart1->FitPosition[1] = line1[1];
    fitPart1->FitPosition[2] = line1[2];
    TVector3 dir1(line1[3], line1[4], line1[5]);
    if(dir1.Mag() > 0) dir1 = dir1.Unit();
    fitPart1->FitDirection[0] = dir1.X();
    fitPart1->FitDirection[1] = dir1.Y();
    fitPart1->FitDirection[2] = dir1.Z();

    // Extract line parameters for daughter 2
    std::vector<double> line2 = vertex->FittedLineParams[1];
    AnaVertexFittedParticlePD* fitPart2 = new AnaVertexFittedParticlePD();
    fitPart2->FitPosition[0] = line2[0];
    fitPart2->FitPosition[1] = line2[1];
    fitPart2->FitPosition[2] = line2[2];
    TVector3 dir2(line2[3], line2[4], line2[5]);
    if(dir2.Mag() > 0) dir2 = dir2.Unit();
    fitPart2->FitDirection[0] = dir2.X();
    fitPart2->FitDirection[1] = dir2.Y();
    fitPart2->FitDirection[2] = dir2.Z();

    vertex->FitParticles.push_back(fitPart1);
    vertex->FitParticles.push_back(fitPart2);

    // Calculate vertex fit direction (sum of daughter fit directions, normalized)
    TVector3 fitDirSum = dir1 + dir2;
    if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
    vertex->FitDirection[0] = fitDirSum.X();
    vertex->FitDirection[1] = fitDirSum.Y();
    vertex->FitDirection[2] = fitDirSum.Z();
  }

  // Calculate Pandora-based vertex position
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);
  CalculatePandoraVertexPosition(part1, part2, vertex->PositionPandora);

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
    // Mark vertex as invalid by setting position to -999
    vertex->Position[0] = -999.0;
    vertex->Position[1] = -999.0;
    vertex->Position[2] = -999.0;
    vertex->Score = 9999.0;
    return 9999.0;
  }

  // Set vertex position from Kalman fit
  vertex->Position[0] = vertexState.position.X();
  vertex->Position[1] = vertexState.position.Y();
  vertex->Position[2] = vertexState.position.Z();

  // Set score to chi2 (lower is better)
  vertex->Score = vertexState.chi2;

  // Calculate minimum distance between fitted tracks at vertex
  TVector3 vtxPos = vertexState.position;
  TVector3 track1ToVtx = vtxPos - track1.position;
  TVector3 track2ToVtx = vtxPos - track2.position;
  double dist1 = (track1ToVtx - track1.direction * track1ToVtx.Dot(track1.direction)).Mag();
  double dist2 = (track2ToVtx - track2.direction * track2ToVtx.Dot(track2.direction)).Mag();
  vertex->MinimumDistance = sqrt(dist1*dist1 + dist2*dist2);

  // Create fitted particle objects from Kalman track states
  AnaVertexFittedParticlePD* fitPart1 = new AnaVertexFittedParticlePD();
  fitPart1->FitPosition[0] = vertexState.position.X();
  fitPart1->FitPosition[1] = vertexState.position.Y();
  fitPart1->FitPosition[2] = vertexState.position.Z();
  TVector3 dir1 = track1.direction.Unit();
  fitPart1->FitDirection[0] = dir1.X();
  fitPart1->FitDirection[1] = dir1.Y();
  fitPart1->FitDirection[2] = dir1.Z();

  AnaVertexFittedParticlePD* fitPart2 = new AnaVertexFittedParticlePD();
  fitPart2->FitPosition[0] = vertexState.position.X();
  fitPart2->FitPosition[1] = vertexState.position.Y();
  fitPart2->FitPosition[2] = vertexState.position.Z();
  TVector3 dir2 = track2.direction.Unit();
  fitPart2->FitDirection[0] = dir2.X();
  fitPart2->FitDirection[1] = dir2.Y();
  fitPart2->FitDirection[2] = dir2.Z();

  vertex->FitParticles.push_back(fitPart1);
  vertex->FitParticles.push_back(fitPart2);

  // Calculate vertex fit direction (sum of daughter fit directions, normalized)
  TVector3 fitDirSum = dir1 + dir2;
  if(fitDirSum.Mag() > 0) fitDirSum = fitDirSum.Unit();
  vertex->FitDirection[0] = fitDirSum.X();
  vertex->FitDirection[1] = fitDirSum.Y();
  vertex->FitDirection[2] = fitDirSum.Z();

  // Calculate Pandora-based vertex position
  CalculatePandoraVertexPosition(daughter1, daughter2, vertex->PositionPandora);

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

  // Copy fit results from reconstructed vertex
  trueEquivalentVertex->FitParticles = vertex->FitParticles;
  trueEquivalentVertex->FitDirection[0] = vertex->FitDirection[0];
  trueEquivalentVertex->FitDirection[1] = vertex->FitDirection[1];
  trueEquivalentVertex->FitDirection[2] = vertex->FitDirection[2];

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
void pdNeutralUtils::CalculateNeutralScore(
    AnaNeutralParticlePD* neutralParticle,
    AnaVertexPD* vertex,
    AnaParticlePD* parentParticle,
    AnaEventB& event) {
//***************************************************************

  if (!neutralParticle || !vertex || !parentParticle) {
    return;
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

  // A) Hit count score - normalized by expected hits for charged particle
  // Expected hits = (length / pitch) * (radius / wire_spacing)
  // Typical pitch ~0.3 cm, wire spacing ~0.5 cm for collection plane
  double wireSpacing = 0.5; // cm
  double pitch = 0.3; // cm
  double expectedHits = (neutralLength / pitch) * (cylinderRadius / wireSpacing);
  if(expectedHits < 1.0) expectedHits = 1.0; // Avoid division by zero

  double hitScore = TMath::Min(nHitsInCylinder / expectedHits, 1.0) * 100.0;

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

  // === 5. OPTIMIZED WEIGHTING ===

  // Physics-motivated weights:
  // - Hit count: most discriminating (50%)
  // - Alignment: important for charged tracks (25%)
  // - Distance: secondary (geometry-dependent) (20%)
  // - Clustering: additional info (5%)

  double neutralScore = 0.50 * hitScore
                      + 0.25 * alignmentScore
                      + 0.20 * distanceScore
                      + 0.05 * clusterScore;

  // Store all metrics in neutral particle
  neutralParticle->NeutralScore = neutralScore;
  neutralParticle->HitsAlignment = hitsAlignment;
  neutralParticle->NHitsInCylinder = nHitsInCylinder;
  neutralParticle->HitsAvgDistance = avgDistance;
  neutralParticle->HitsRMSDistance = rmsDistance;
  neutralParticle->HitsLongitudinalSpan = longitudinalSpan;
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

  // Calculate invariant mass for true equivalent neutral particle
  Float_t trueInvariantMass = -999;
  if (vertex->TrueEquivalentVertex && vertex->TrueEquivalentVertex->TrueParticles.size() >= 2) {
    AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[0]);
    AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(vertex->TrueEquivalentVertex->TrueParticles[1]);

    if (trueParticle1 && trueParticle2 &&
        trueParticle1->Momentum > 0 && trueParticle2->Momentum > 0 &&
        trueParticle1->Momentum != -999 && trueParticle2->Momentum != -999) {
      const Float_t pionMass = 0.13957;
      trueInvariantMass = pdAnaUtils::ComputeTrueInvariantMass(*trueParticle1, *trueParticle2, pionMass, pionMass);
    }
  }

  trueEquivalentNeutralParticle->Mass = trueInvariantMass;

  return trueEquivalentNeutralParticle;
}

//***************************************************************
AnaNeutralParticlePD* pdNeutralUtils::CreateNeutral(AnaEventB& event, AnaVertexPD* vertex, int neutralParticleID){
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return nullptr;
  }

  // Find the parent particle using the vertex's ParentID
  AnaParticlePD* parentParticle = nullptr;

  // First search in regular particles
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part && part->UniqueID == vertex->ParentID){
      parentParticle = part;
      break;
    }
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

  // Set positions
  neutralParticle->PositionStart[0] = parentParticle->PositionEnd[0];
  neutralParticle->PositionStart[1] = parentParticle->PositionEnd[1];
  neutralParticle->PositionStart[2] = parentParticle->PositionEnd[2];
  neutralParticle->PositionStart[3] = parentParticle->PositionEnd[3];

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

      // Create FitParent from extrapolated parent track
      AnaVertexFittedParticlePD* fitParent = new AnaVertexFittedParticlePD();
      fitParent->FitPosition[0] = parentLineParams[0];
      fitParent->FitPosition[1] = parentLineParams[1];
      fitParent->FitPosition[2] = parentLineParams[2];
      TVector3 parentDir(parentLineParams[3], parentLineParams[4], parentLineParams[5]);
      if(parentDir.Mag() > 0) parentDir = parentDir.Unit();
      fitParent->FitDirection[0] = parentDir.X();
      fitParent->FitDirection[1] = parentDir.Y();
      fitParent->FitDirection[2] = parentDir.Z();

      neutralParticle->FitParent = fitParent;
    }
  }

  // Calculate NPotentialParents and NRecoHitsInVertex
  // These count particles and hits within a cylinder from parent end to vertex
  Int_t NPotentialParents = 0;
  Int_t NRecoHitsInVertex = 0;

  // Get parameters for the cylinder
  double daughterDistance = ND::params().GetParameterD("neutralKaonAnalysis.DaughterDistance");
  double cylinderRadius = ND::params().GetParameterD("neutralKaonAnalysis.CylinderRadius");

  // Reuse parts and nParts variables already declared at the beginning of this function
  TVector3 vertexPos(vertex->Position[0], vertex->Position[1], vertex->Position[2]);

  // Loop over all particles to find potential parents
  for(int j = 0; j < nParts; j++){
    AnaParticlePD* potentialParent = static_cast<AnaParticlePD*>(parts[j]);
    if(!potentialParent) continue;

    // Skip if the potential parent is one of the particles in the vertex or the current parent particle
    if(vertex->Particles.size() >= 2){
      if(potentialParent->UniqueID == vertex->Particles[0]->UniqueID ||
         potentialParent->UniqueID == vertex->Particles[1]->UniqueID ||
         potentialParent->UniqueID == parentParticle->UniqueID) continue;
    }

    TVector3 potentialParentEnd(potentialParent->PositionEnd[0], potentialParent->PositionEnd[1], potentialParent->PositionEnd[2]);
    double potentialParentDistance = (potentialParentEnd - vertexPos).Mag();

    if(potentialParentDistance < daughterDistance){
      NPotentialParents++;

      // Count hits from this particle within the cylinder
      Int_t nHits = potentialParent->Hits[2].size();
      for(size_t k = 0; k < nHits; k++){
        AnaHitPD* hit = &potentialParent->Hits[2][k];

        // Define cylinder from parent->PositionEnd to vertex->Position
        TVector3 cylinderStart(parentParticle->PositionEnd[0], parentParticle->PositionEnd[1], parentParticle->PositionEnd[2]);
        TVector3 cylinderEnd(vertex->Position[0], vertex->Position[1], vertex->Position[2]);
        TVector3 cylinderAxis = cylinderEnd - cylinderStart;
        double cylinderLength = cylinderAxis.Mag();

        if(cylinderLength <= 0) continue; // Skip if cylinder has no length

        // Get hit position
        TVector3 hitPos = hit->Position;

        // Vector from cylinder start to hit
        TVector3 startToHit = hitPos - cylinderStart;

        // Project onto cylinder axis to check if hit is within cylinder length
        double projection = startToHit.Dot(cylinderAxis) / cylinderLength;

        // Check if projection is within cylinder bounds (0 to cylinderLength)
        if(projection >= 0 && projection <= cylinderLength){
          // Calculate perpendicular distance from hit to cylinder axis
          TVector3 axisDirection = cylinderAxis.Unit();
          TVector3 projectionVector = axisDirection * startToHit.Dot(axisDirection);
          TVector3 perpendicularVector = startToHit - projectionVector;
          double perpendicularDistance = perpendicularVector.Mag();

          // Check if hit is within cylinder radius
          if(perpendicularDistance < cylinderRadius){
            NRecoHitsInVertex++;
          }
        }
      }
    }
  }

  neutralParticle->Vertex->NPotentialParents = NPotentialParents;
  neutralParticle->NRecoHitsInVertex = NRecoHitsInVertex;

  // Calculate neutral particle score using helper function
  CalculateNeutralScore(neutralParticle, vertex, parentParticle, event);

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

  // Loop over all vertices and create one neutral particle per vertex
  for(size_t v = 0; v < vertices.size(); v++){
    AnaVertexPD* vertex = vertices[v];
    if(!vertex) continue;

    AnaNeutralParticlePD* neutralParticle = CreateNeutral(event, vertex, neutralParticleID);
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

