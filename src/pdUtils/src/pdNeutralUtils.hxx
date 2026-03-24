#ifndef pdNeutralUtils_h
#define pdNeutralUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include "pdAnnihilationUtils.hxx"
#include "pdCreationUtils.hxx"
#include "pdNeutralHelpers.hxx"
#include <vector>
#include <unordered_map>
#include <utility>
#include <set>
#include "TVector3.h"

namespace pdNeutralUtils {

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

  // Create neutral particles from all vertices
  std::vector<AnaNeutralParticlePD*> CreateNeutrals(AnaEventB& event,
                                                     const std::vector<AnaCreationVertexPD*>& creationVertices,
                                                     const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices);

  // Filter neutral particles ensuring each annihilation vertex belongs to at most one neutral particle
  std::vector<AnaNeutralParticlePD*> FilterNeutralsByScore(std::vector<AnaNeutralParticlePD*>& neutralParticles);

  // Calculate raw score components for a neutral particle (ProtonScore, alignment, NHitsInCylinder)
  // Returns: {protonScore, alignmentScore, hitsScore}
  std::tuple<double, double, double> CalculateRawScoreComponents(AnaNeutralParticlePD* neutralParticle);

  // Normalize scores using percentile-based normalization across all neutral particles in the event
  // Each component is normalized to [0,1] where best (lowest) gets 1.0 and worst (highest) gets 0.0
  void NormalizeNeutralParticleScores(std::vector<AnaNeutralParticlePD*>& neutralParticles);

  // Calculate degeneracies for both creation and annihilation vertices
  void CalculateVertexDegeneracies(AnaEventB& event,
                                     AnaCreationVertexPD* creationVertex,
                                     AnaAnnihilationVertexPD* annihilationVertex);

}

#endif

