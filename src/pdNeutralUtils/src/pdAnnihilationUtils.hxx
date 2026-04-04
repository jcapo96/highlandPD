#ifndef pdAnnihilationUtils_h
#define pdAnnihilationUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include <vector>

namespace pdAnnihilationUtils {

  // Create vertices with simplified geometric position assignment.
  // Optional outputs return counts before/after overlap filtering.
  std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance = 5.0,
                                                       Int_t* nBeforeFiltering = nullptr,
                                                       Int_t* nAfterFiltering = nullptr);

  // Common vertex creation logic based only on pair conditions.
  std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxDaughterDistance,
                                                             Int_t* nBeforeFiltering = nullptr,
                                                             Int_t* nAfterFiltering = nullptr);

  // Fill PositionPandora from the two daughter start-point/direction lines.
  void FillPositionPandora(AnaAnnihilationVertexPD* vertex);

  // Fill PositionFit from line fits in [trackFitDistanceFromStart, trackFitDistanceFromStart + trackFitLength].
  void FillPositionFit(AnaAnnihilationVertexPD* vertex, double trackFitLength, double trackFitDistanceFromStart = 0.0);

  // Keep vertices with smallest fit minimum distance first, with each particle used at most once.
  std::vector<AnaAnnihilationVertexPD*> FilterVerticesByMinimumDistanceFit(std::vector<AnaAnnihilationVertexPD*>& vertices);

  // Angle (rad) between neutral creation→annihilation axis and reconstructed vertex momentum (Pandora / fit).
  void FillNeutralParticleAlignment(AnaNeutralParticlePD* neutral, const AnaEventB& event, double trackFitLength,
                                    double trackFitDistanceFromStart);

}

#endif