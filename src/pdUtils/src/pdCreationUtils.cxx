#include "pdCreationUtils.hxx"
#include "pdNeutralHelpers.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "BasicUtils.hxx"
#include "TVector3.h"
#include <cmath>
#include <set>
#include <vector>
#include <map>
#include <algorithm>

namespace pdCreationUtils {

//***************************************************************
void CalculateCreationVertexScores(AnaCreationVertexPD* creationVtx) {
//***************************************************************

  if (!creationVtx || !creationVtx->BeamParticle) {
    return;
  }

  AnaParticlePD* beam = creationVtx->BeamParticle;
  AnaParticlePD* second = creationVtx->SecondParticle;

  // Handle case with no second particle (beam-only vertex)
  if (!second) {
    // Set position to beam end position
    creationVtx->Position[0] = beam->PositionEnd[0];
    creationVtx->Position[1] = beam->PositionEnd[1];
    creationVtx->Position[2] = beam->PositionEnd[2];

    // Set default scores
    creationVtx->ProtonScore = 999.0; // No proton candidate
    creationVtx->DistanceScore = 0.0; // No distance (vertex at beam end)
    creationVtx->MinDistanceScore = 0.0;
    // Set Score to 999.0 so these vertices are filtered out
    creationVtx->Score = 999.0;

    return;
  }

  // 1. ProtonScore: Chi2/ndf under proton hypothesis for secondary particle
  std::pair<double, int> chi2Result = pdAnaUtils::Chi2PID(*second, 2212); // 2212 = proton PDG
  creationVtx->ProtonScore = (chi2Result.second > 0) ?
                             static_cast<float>(chi2Result.first) / static_cast<float>(chi2Result.second) : 999.0;

  // 2. DistanceScore: Distance from beam end to secondary particle start
  TVector3 beamEnd(beam->PositionEnd[0], beam->PositionEnd[1], beam->PositionEnd[2]);
  TVector3 secondStart(second->PositionStart[0], second->PositionStart[1], second->PositionStart[2]);
  creationVtx->DistanceScore = (beamEnd - secondStart).Mag();

  // 3 & 4. Calculate creation vertex position using Pandora directions (consistent with annihilation vertex)
  // Use Pandora PositionEnd and DirectionEnd for beam (similar to how annihilation uses PositionStart/DirectionStart)
  TVector3 beamPandoraPos(beam->PositionEnd[0], beam->PositionEnd[1], beam->PositionEnd[2]);
  TVector3 beamPandoraDir(beam->DirectionEnd[0], beam->DirectionEnd[1], beam->DirectionEnd[2]);

  TVector3 secondPandoraPos(second->PositionStart[0], second->PositionStart[1], second->PositionStart[2]);
  TVector3 secondPandoraDir(second->DirectionStart[0], second->DirectionStart[1], second->DirectionStart[2]);

  // Normalize directions
  if (beamPandoraDir.Mag() > 0) beamPandoraDir = beamPandoraDir.Unit();
  if (secondPandoraDir.Mag() > 0) secondPandoraDir = secondPandoraDir.Unit();

  // Use helper function to calculate midpoint of minimum distance
  TVector3 vertexPosition = pdNeutralHelpers::CalculateLineMinDistanceMidpoint(
      beamPandoraPos, beamPandoraDir,
      secondPandoraPos, secondPandoraDir);

  creationVtx->Position[0] = vertexPosition.X();
  creationVtx->Position[1] = vertexPosition.Y();
  creationVtx->Position[2] = vertexPosition.Z();

  // Calculate closest points manually for storage
  TVector3 w0 = beamPandoraPos - secondPandoraPos;
  double a = beamPandoraDir.Dot(beamPandoraDir);
  double b = beamPandoraDir.Dot(secondPandoraDir);
  double c = secondPandoraDir.Dot(secondPandoraDir);
  double d = beamPandoraDir.Dot(w0);
  double e = secondPandoraDir.Dot(w0);
  double denom = a * c - b * b;

  if (fabs(denom) > 1e-6) {
    // Lines are not parallel - calculate closest points
    double s = (b * e - c * d) / denom;
    double t = (a * e - b * d) / denom;

    TVector3 closestPoint1 = beamPandoraPos + s * beamPandoraDir;
    TVector3 closestPoint2 = secondPandoraPos + t * secondPandoraDir;

    // Store closest points for downstream consumers
    creationVtx->ClosestPointBeam[0] = closestPoint1.X();
    creationVtx->ClosestPointBeam[1] = closestPoint1.Y();
    creationVtx->ClosestPointBeam[2] = closestPoint1.Z();
    creationVtx->ClosestPointSecond[0] = closestPoint2.X();
    creationVtx->ClosestPointSecond[1] = closestPoint2.Y();
    creationVtx->ClosestPointSecond[2] = closestPoint2.Z();

    // Store the minimum distance as MinDistanceScore
    double minDistance = (closestPoint1 - closestPoint2).Mag();
    creationVtx->MinDistanceScore = static_cast<Float_t>(minDistance);
    // Set Score to MinimumDistance (same as annihilation vertices)
    creationVtx->Score = static_cast<Float_t>(minDistance);

    // Store Pandora-based line parameters for event display (consistent with visualization)
    std::vector<double> beamLineParams = {
      beam->PositionEnd[0], beam->PositionEnd[1], beam->PositionEnd[2],
      beam->DirectionEnd[0], beam->DirectionEnd[1], beam->DirectionEnd[2]
    };
    std::vector<double> secondLineParams = {
      second->PositionStart[0], second->PositionStart[1], second->PositionStart[2],
      second->DirectionStart[0], second->DirectionStart[1], second->DirectionStart[2]
    };

    creationVtx->FittedLineParams.clear();
    creationVtx->FittedLineParams.push_back(beamLineParams);    // Line 0: beam particle (uses endPos/endDir)
    creationVtx->FittedLineParams.push_back(secondLineParams);   // Line 1: second particle (uses startPos/startDir)
  } else {
    // Lines are parallel - fallback to simple midpoint
    TVector3 midpoint = 0.5 * (beamEnd + secondStart);
    creationVtx->Position[0] = midpoint.X();
    creationVtx->Position[1] = midpoint.Y();
    creationVtx->Position[2] = midpoint.Z();
    creationVtx->MinDistanceScore = creationVtx->DistanceScore;
    creationVtx->ClosestPointBeam[0] = midpoint.X();
    creationVtx->ClosestPointBeam[1] = midpoint.Y();
    creationVtx->ClosestPointBeam[2] = midpoint.Z();
    creationVtx->ClosestPointSecond[0] = midpoint.X();
    creationVtx->ClosestPointSecond[1] = midpoint.Y();
    creationVtx->ClosestPointSecond[2] = midpoint.Z();
    // Set Score to DistanceScore (fallback when lines are parallel)
    creationVtx->Score = creationVtx->DistanceScore;
  }
}

//***************************************************************
std::vector<AnaCreationVertexPD*> CreateCreationVertices(
    AnaEventB& event,
    double creationVertexRadius,
    const std::set<Int_t>& excludeParticleIDs) {
//***************************************************************

  std::vector<AnaCreationVertexPD*> creationVertices;

  // Get all particles from event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // Loop through all particles to create creation vertices for each
  for (int i = 0; i < nParts; i++) {
    AnaParticlePD* parentParticle = static_cast<AnaParticlePD*>(parts[i]);
    if (!parentParticle) continue;

    // Only seed creation vertices from the Pandora beam particle.
    if (!parentParticle->isPandora) continue;

    // Skip excluded particles
    if (excludeParticleIDs.find(parentParticle->UniqueID) != excludeParticleIDs.end()) {
      continue;
    }

    // Check particle has valid end position
    if (parentParticle->PositionEnd[0] < -900) continue;

    TVector3 parentEnd(parentParticle->PositionEnd[0],
                       parentParticle->PositionEnd[1],
                       parentParticle->PositionEnd[2]);

    // Loop through all particles to find potential secondary particles for this parent
    bool foundSecondParticle = false;
    for (int j = 0; j < nParts; j++) {
      AnaParticlePD* secondaryParticle = static_cast<AnaParticlePD*>(parts[j]);
      if (!secondaryParticle) continue;

      // Skip the parent particle itself
      if (secondaryParticle->UniqueID == parentParticle->UniqueID) continue;

      // Skip excluded particles
      if (excludeParticleIDs.find(secondaryParticle->UniqueID) != excludeParticleIDs.end()) {
        continue;
      }

      // Check secondary particle has valid start position
      if (secondaryParticle->PositionStart[0] < -900) continue;

      TVector3 secondaryStart(secondaryParticle->PositionStart[0],
                              secondaryParticle->PositionStart[1],
                              secondaryParticle->PositionStart[2]);

      // Check if secondary particle start is within radius of parent end
      double distance = (parentEnd - secondaryStart).Mag();

      if (distance < creationVertexRadius) {
        // Create creation vertex candidate
        AnaCreationVertexPD* creationVtx = new AnaCreationVertexPD();
        creationVtx->BeamParticle = parentParticle;
        creationVtx->SecondParticle = secondaryParticle;

        // Calculate all scores
        CalculateCreationVertexScores(creationVtx);

        creationVertices.push_back(creationVtx);
        foundSecondParticle = true;
      }
    }

    // If no secondary particles were found, create a parent-only creation vertex
    if (!foundSecondParticle) {
      AnaCreationVertexPD* creationVtx = new AnaCreationVertexPD();
      creationVtx->BeamParticle = parentParticle;
      creationVtx->SecondParticle = nullptr; // No second particle

      // Calculate scores (will use parent end position for vertex)
      CalculateCreationVertexScores(creationVtx);

      creationVertices.push_back(creationVtx);
    }
  }

  // Minimal branch: disable score/uniqueness filtering at creation-vertex stage.
  return creationVertices;
}

//***************************************************************
std::vector<AnaCreationVertexPD*> FilterCreationVerticesByScore(std::vector<AnaCreationVertexPD*>& vertices) {
//***************************************************************

  // Sort vertices by MinDistanceScore first (ascending - lower is better)
  std::sort(vertices.begin(), vertices.end(),
    [](const AnaCreationVertexPD* a, const AnaCreationVertexPD* b) {
      // Handle invalid MinDistanceScores
      if (a->MinDistanceScore < -900 && b->MinDistanceScore < -900) return false;
      if (a->MinDistanceScore < -900) return false;
      if (b->MinDistanceScore < -900) return true;
      return a->MinDistanceScore < b->MinDistanceScore;
    });

  // Track which vertex each particle is associated with
  std::map<AnaParticlePD*, AnaCreationVertexPD*> particleToVertex;

  // Output vector for selected vertices
  std::vector<AnaCreationVertexPD*> selectedVertices;
  std::set<AnaCreationVertexPD*> selectedSet;  // For quick lookup

  // Threshold for considering two MinDistanceScores as equivalent (1cm = 1.0)
  const double minDistanceThreshold = 2.5;

  // Iterate through sorted vertices (best MinDistanceScore first)
  for (AnaCreationVertexPD* vertex : vertices) {
    if (!vertex || !vertex->BeamParticle) continue;

    AnaParticlePD* beamParticle = vertex->BeamParticle;
    AnaParticlePD* secondParticle = vertex->SecondParticle;

    // Check if beam particle is already associated with a vertex
    auto beamIt = particleToVertex.find(beamParticle);
    bool beamUsed = (beamIt != particleToVertex.end());

    // Check if second particle is already associated with a vertex (if it exists)
    auto secondIt = particleToVertex.end();
    bool secondUsed = false;
    if (secondParticle) {
      secondIt = particleToVertex.find(secondParticle);
      secondUsed = (secondIt != particleToVertex.end());
    }

    // If neither particle is used, keep this vertex
    if (!beamUsed && !secondUsed) {
      selectedVertices.push_back(vertex);
      selectedSet.insert(vertex);
      particleToVertex[beamParticle] = vertex;
      if (secondParticle) {
        particleToVertex[secondParticle] = vertex;
      }
    }
    // If a particle is already used, check if we should replace the existing vertex
    else {
      // Check if we should replace based on MinDistanceScore difference and ProtonScore
      bool shouldReplace = false;
      AnaCreationVertexPD* vertexToReplace = nullptr;

      if (beamUsed) {
        AnaCreationVertexPD* existingVertex = beamIt->second;
        double minDistDiff = std::abs(vertex->MinDistanceScore - existingVertex->MinDistanceScore);

        // If MinDistanceScores are within 1cm, prefer lower ProtonScore (more proton-like)
        if (minDistDiff < minDistanceThreshold) {
          // Both ProtonScores should be valid to compare
          if (vertex->ProtonScore >= -900 && existingVertex->ProtonScore >= -900) {
            if (vertex->ProtonScore < existingVertex->ProtonScore) {
              shouldReplace = true;
              vertexToReplace = existingVertex;
            }
          }
        }
      }

      if (secondUsed && !shouldReplace) {
        AnaCreationVertexPD* existingVertex = secondIt->second;
        double minDistDiff = std::abs(vertex->MinDistanceScore - existingVertex->MinDistanceScore);

        // If MinDistanceScores are within 1cm, prefer lower ProtonScore (more proton-like)
        if (minDistDiff < minDistanceThreshold) {
          // Both ProtonScores should be valid to compare
          if (vertex->ProtonScore >= -900 && existingVertex->ProtonScore >= -900) {
            if (vertex->ProtonScore < existingVertex->ProtonScore) {
              shouldReplace = true;
              vertexToReplace = existingVertex;
            }
          }
        }
      }

      if (shouldReplace && vertexToReplace) {
        // Remove the old vertex from selected set
        selectedSet.erase(vertexToReplace);
        selectedVertices.erase(
          std::remove(selectedVertices.begin(), selectedVertices.end(), vertexToReplace),
          selectedVertices.end()
        );

        // Remove old particle associations
        if (vertexToReplace->BeamParticle) {
          particleToVertex.erase(vertexToReplace->BeamParticle);
        }
        if (vertexToReplace->SecondParticle) {
          particleToVertex.erase(vertexToReplace->SecondParticle);
        }

        // Add the new vertex
        selectedVertices.push_back(vertex);
        selectedSet.insert(vertex);
        particleToVertex[beamParticle] = vertex;
        if (secondParticle) {
          particleToVertex[secondParticle] = vertex;
        }
      }
    }
  }

  // Delete vertices that were not selected
  for (AnaCreationVertexPD* vertex : vertices) {
    if (selectedSet.find(vertex) == selectedSet.end()) {
      delete vertex;
    }
  }

  return selectedVertices;
}


} // namespace pdCreationUtils

