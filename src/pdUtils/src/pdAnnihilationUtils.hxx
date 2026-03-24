#ifndef pdAnnihilationUtils_h
#define pdAnnihilationUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include "pdNeutralHelpers.hxx"
#include <vector>

namespace pdAnnihilationUtils {

  // Create vertices - checks parameter to decide which algorithm to use
  std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance = 5.0);

  // Common vertex creation logic with selectable position finder
  std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxDaughterDistance, double (*positionFinder)(AnaVertexPD*, double));

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

  // Find vertex position by fitting lines to daughter particles and finding closest points
  void FindVertexPositionFit(AnaVertexPD* vertex);

  // Find vertex position by fitting lines to true daughter particles and finding closest points
  void FindVertexPositionFit(AnaTrueEquivalentVertexPD* vertex);

}

#endif

