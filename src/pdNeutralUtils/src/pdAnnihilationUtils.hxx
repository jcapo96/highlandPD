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

  // Geometry of the two daughter lines at the given endpoints (line fits as in FillPositionFit, the Pandora endpoint
  // and direction as fallback): distance of closest approach, cosine of the opening angle (both lines pointing into
  // their tracks) and, per daughter, the gap: the distance along its line from the closest-approach point to its
  // endpoint, positive when the vertex lies upstream of the endpoint (a track whose first segment is missing).
  struct PairLineGeometry {
    Float_t doca = -999.f;
    Float_t gap1 = -999.f;
    Float_t gap2 = -999.f;
    Float_t cosOpening = -999.f;
    bool valid = false;
  };
  PairLineGeometry ComputePairLineGeometry(const AnaParticlePD* d1, const AnaParticlePD* d2, bool useStart1,
                                           bool useStart2, double trackFitLength, double trackFitDistanceFromStart);

  // Geometric pairing (AnnihilationVertexPairingMode = 1): closest approach <= AnnihilationVertexGeomMaxDoca, both
  // gaps in [AnnihilationVertexGeomMinGap, AnnihilationVertexGeomMaxGap] and cosine of the opening angle below
  // AnnihilationVertexGeomMaxCosOpening (the closest-approach point is undefined for nearly parallel lines).
  struct GeometricPairingCuts {
    double maxDoca = 4.0;
    double minGap = -3.0;
    double maxGap = 60.0;
    double maxCosOpening = 0.998;
  };
  GeometricPairingCuts ReadGeometricPairingCuts();
  bool PassesGeometricPairing(const PairLineGeometry& geometry, const GeometricPairingCuts& cuts);

  // Endpoint combination (0=SS, 1=SE, 2=ES, 3=EE; SS only without endpoint combinatorics) that passes the geometric
  // pairing with the smallest gap sum, or -1. geometry is that combination's, or the one of the combination with the
  // smallest endpoint distance when none passes.
  Int_t BestGeometricPairing(const AnaParticlePD* d1, const AnaParticlePD* d2, bool useEndpointCombinatorics,
                             double trackFitLength, double trackFitDistanceFromStart, const GeometricPairingCuts& cuts,
                             PairLineGeometry& geometry);

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