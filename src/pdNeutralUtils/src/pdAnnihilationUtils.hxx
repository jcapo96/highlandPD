#ifndef pdAnnihilationUtils_h
#define pdAnnihilationUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include <vector>

namespace pdAnnihilationUtils {

  // Create annihilation-vertex candidates from daughter pairs passing the endpoint-distance cut
  // (minimum of SS/SE/ES/EE within maxDaughterDistance). Reversal: SS none; SE/ES/EE flip each daughter
  // if PositionStart(Z)>PositionEnd(Z), swapping endpoints/directions and reversing TrjPoints and Hits order.
  // Optional outputs return counts before/after overlap filtering by smallest MinimumDistanceFit.
  std::vector<AnaAnnihilationVertexPD*> CreateVertices(AnaEventB& event, double maxDaughterDistance = 5.0,
                                                       Int_t* nBeforeFiltering = nullptr,
                                                       Int_t* nAfterFiltering = nullptr);

  // Common vertex creation logic using the same pairing cut and fit-distance overlap filtering.
  std::vector<AnaAnnihilationVertexPD*> CreateVerticesCommon(AnaEventB& event, double maxDaughterDistance,
                                                             Int_t* nBeforeFiltering = nullptr,
                                                             Int_t* nAfterFiltering = nullptr);

  // Plane used to judge a vertex daughter: the collection plane; with inductionFallback, a daughter with no more than
  // minCollectionHits collection hits is judged on its induction plane with more hits (filled only from PDSPAnalyzer
  // ntuples with pdIO.PDSPAnalyzerFillAllPlanes = 1). Hits counted by the hit requirement, and Chi2PID(pdg)/npts.
  Int_t VertexDaughterHits(const AnaParticlePD* p, bool inductionFallback, int minCollectionHits);
  Float_t VertexDaughterChi2PerHit(const AnaParticlePD* p, int pdg, bool inductionFallback, int minCollectionHits);

  // Fill PositionPandora from the two daughter start-point/direction lines.
  void FillPositionPandora(AnaAnnihilationVertexPD* vertex);

  // Fill PositionFit from line fits in [trackFitDistanceFromStart, trackFitDistanceFromStart + trackFitLength],
  // measured from the start (useStartN = true) or from the end of daughter N.
  void FillPositionFit(AnaAnnihilationVertexPD* vertex, double trackFitLength, double trackFitDistanceFromStart = 0.0,
                       bool useStart1 = true, bool useStart2 = true);

  // Keep vertices with smallest MinimumDistanceFit first. overlapMode 0: each particle used at most once;
  // 1: a particle may be shared by several vertices as long as all of them pair it by the same endpoint.
  std::vector<AnaAnnihilationVertexPD*> FilterVerticesByMinimumDistanceFit(std::vector<AnaAnnihilationVertexPD*>& vertices,
                                                                           int overlapMode = 0);

  // Angle (rad) between neutral creation→annihilation axis and reconstructed vertex momentum (Pandora / fit).
  void FillNeutralParticleAlignment(AnaNeutralParticlePD* neutral, const AnaEventB& event, double trackFitLength,
                                    double trackFitDistanceFromStart);

  // Assign particle momentum with the same free-range dE/dx fit used for annihilation daughters,
  // but under a proton hypothesis (PDG 2212). Returns method code as in DaughterMomentumMethod.
  Int_t AssignProtonMomentumFromResidualRange(AnaParticlePD* particle);

  // Compute annihilation-vertex degeneracy using configured parameters while optionally
  // excluding one reco particle by UniqueID (e.g., neutral parent candidate).
  Int_t ComputeAnnihilationVertexDegeneracyWithExclusion(const AnaEventB& event,
                                                         const AnaAnnihilationVertexPD* vertex,
                                                         Int_t excludedParticleUniqueID = -1,
                                                         const AnaCreationVertexPD* creationVertex = nullptr);

  // Compute creation-vertex degeneracy with the same geometric algorithm used for
  // annihilation vertices. The annihilation daughters are excluded from counting.
  Int_t ComputeCreationVertexDegeneracy(const AnaEventB& event,
                                        const AnaCreationVertexPD* creationVertex,
                                        const AnaAnnihilationVertexPD* annihilationVertex,
                                        Int_t excludedParticleUniqueID = -1);

}

#endif