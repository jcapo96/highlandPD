#ifndef pdCreationUtils_h
#define pdCreationUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include "pdNeutralHelpers.hxx"
#include <vector>
#include <set>

namespace pdCreationUtils {

  // Calculate scores for a creation vertex (ProtonScore, DistanceScore, MinDistanceScore, Position)
  void CalculateCreationVertexScores(AnaCreationVertexPD* creationVtx);

  // Create all valid creation vertices for an event
  // excludeParticleIDs: particle IDs to exclude (e.g., annihilation vertex particles)
  std::vector<AnaCreationVertexPD*> CreateCreationVertices(
      AnaEventB& event,
      double creationVertexRadius,
      const std::set<Int_t>& excludeParticleIDs = {});

  // Filter creation vertices ensuring each particle belongs to at most one vertex
  // Sorts by Score (MinimumDistance) and keeps only the best vertex for each particle
  std::vector<AnaCreationVertexPD*> FilterCreationVerticesByScore(std::vector<AnaCreationVertexPD*>& vertices);

}

#endif

