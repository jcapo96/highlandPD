#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax);

  void FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                Int_t nVerticesBeforeFiltering = -1, Int_t nVerticesAfterFiltering = -1,
                                AnaBeamB* beam = NULL);
  void FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaAnnihilationVertexPD* vertex);
  void WriteHitDistanceProfiles(OutputManager& output);

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Candidate info
    nk0 = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    k0nvtxbeforefiltering, // Number of annihilation vertices before overlap filtering
    k0nvtxafterfiltering,  // Number of annihilation vertices after overlap filtering
    // Active microtree variables (vertex-only set)
    k0vtxtruepos, //K0 vertex true position
    k0vtxoriginaldistance, //K0 vertex original distance
    k0vtxtrueoriginaldistance, //K0 vertex true original distance
    k0vtxpandorapos, //Pandora-based vertex position
    k0vtxfitpos, //Algorithm-specific fitted vertex position
    k0vtxpandorax,
    k0vtxpandoray,
    k0vtxpandoraz,
    k0vtxfitx,
    k0vtxfity,
    k0vtxfitz,
    k0vtxpandoraresidual, //|Pandora vertex position - true vertex position| [cm]
    k0vtxfitresidual, //|Fit vertex position - true vertex position| [cm]
    k0vtxpandoraresidualx, //Pandora x_reco - x_true [cm]
    k0vtxpandoraresidualy, //Pandora y_reco - y_true [cm]
    k0vtxpandoraresidualz, //Pandora z_reco - z_true [cm]
    k0vtxfitresidualx, //Fit x_reco - x_true [cm]
    k0vtxfitresidualy, //Fit y_reco - y_true [cm]
    k0vtxfitresidualz, //Fit z_reco - z_true [cm]
  enumNeutralKaonMicroTreesLast
  };
}

#endif