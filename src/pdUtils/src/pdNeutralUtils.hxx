#ifndef pdNeutralUtils_h
#define pdNeutralUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include <vector>
#include <unordered_map>
#include <utility>
#include <set>
#include "TVector3.h"

namespace pdNeutralUtils {

  // Create vertices - checks parameter to decide which algorithm to use
  std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxVertexRadius = 30.0, double maxDaughterDistance = 5.0);

  // Common vertex creation logic with selectable position finder
  std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxVertexRadius, double maxDaughterDistance, double (*positionFinder)(AnaVertexPD*, double));

  // Position finder functions - return score value
  double FindVertexPositionWithFit(AnaVertexPD* vertex, double trackFitLength);
  double FindVertexPositionGeometric(AnaVertexPD* vertex, double trackFitLength);
  double FindVertexPositionKalman(AnaVertexPD* vertex, double trackFitLength);

  // Validate vertex (position and score checks)
  bool ValidateVertex(AnaVertexPD* vertex);

  // Filter vertices ensuring each particle belongs to at most one vertex
  std::vector<AnaAnnihilationVertexPD*> FilterVerticesByScore(std::vector<AnaAnnihilationVertexPD*>& vertices);

  // Helper function to create and fill true equivalent vertex
  AnaTrueEquivalentVertexPD* FillTrueEquivalentVertex(AnaVertexPD* vertex);

  // Helper function to create and fill true equivalent neutral particle
  AnaTrueEquivalentNeutralParticlePD* FillTrueEquivalentNeutralParticle(
      AnaVertexPD* vertex,
      AnaParticlePD* parentParticle);

  // Calculate neutral particle score and metrics
  // Returns: {NPotentialParents, NRecoHitsInVertex}
  std::pair<Int_t, Int_t> CalculateNeutralScore(
      AnaNeutralParticlePD* neutralParticle,
      AnaVertexPD* vertex,
      AnaParticlePD* parentParticle,
      AnaEventB& event,
      const std::unordered_map<Int_t, AnaParticlePD*>& particleByUniqueID);

  // Create single neutral particle from a vertex
  AnaNeutralParticlePD* CreateNeutral(
      AnaEventB& event,
      AnaVertexPD* vertex,
      int neutralParticleID,
      const std::unordered_map<Int_t, AnaParticlePD*>& particleByUniqueID);

  // Create neutral particles from all vertices
  std::vector<AnaNeutralParticlePD*> CreateNeutrals(AnaEventB& event, const std::vector<AnaAnnihilationVertexPD*>& vertices);

  // Calculate midpoint of minimum distance between two 3D lines
  TVector3 CalculateLineMinDistanceMidpoint(
      const TVector3& pos1, const TVector3& dir1,
      const TVector3& pos2, const TVector3& dir2);

  // Calculate scores for a creation vertex (ProtonScore, DistanceScore, MinDistanceScore, Position)
  void CalculateCreationVertexScores(AnaCreationVertexPD* creationVtx);

  // Create all valid creation vertices for an event
  // excludeParticleIDs: particle IDs to exclude (e.g., annihilation vertex particles)
  std::vector<AnaCreationVertexPD*> CreateCreationVertices(
      AnaEventB& event,
      double creationVertexRadius,
      const std::set<Int_t>& excludeParticleIDs = {});

  // Calculate degeneracies for both creation and annihilation vertices
  void CalculateVertexDegeneracies(AnaEventB& event,
                                     AnaCreationVertexPD* creationVertex,
                                     AnaAnnihilationVertexPD* annihilationVertex);

}

#endif

