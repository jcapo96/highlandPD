#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_K0Particle(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0VtxDaughters(OutputManager& output, UInt_t nmax);

  void FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                Int_t nVerticesBeforeFiltering = -1, Int_t nVerticesAfterFiltering = -1,
                                AnaBeamB* beam = NULL, size_t neutralCandidateIndex = 0);
  void FillNeutralKaonVariables_K0Particle(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
  void FillNeutralKaonVariables_K0vtxDaughters(OutputManager& output, AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
  void WriteHitDistanceProfiles(OutputManager& output);

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Candidate info
    nk0 = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    k0nvtxbeforefiltering, // Number of annihilation vertices before overlap filtering
    k0nvtxafterfiltering,  // Number of annihilation vertices after overlap filtering
    k0lengthpandora, // Neutral length from creation to annihilation Pandora position [cm]
    k0lengthfit, // Neutral length from creation to annihilation Fit position [cm]
    k0truelength, // True K0 length from true creation to true decay vertex [cm]
    k0alignmentpandora, // Angle (rad): neutral axis vs vertex momentum (Pandora dirs)
    k0alignmentfit, // Angle (rad): neutral axis vs vertex momentum (fit dirs)
    // Active microtree variables (vertex-only set)
    k0vtxtruepos, //K0 vertex true position
    k0vtxoriginaldistance, //K0 vertex original distance
    k0vtxtrueoriginaldistance, //K0 vertex true original distance
    k0vtxdegeneracy, //K0 annihilation-vertex degeneracy (Reco)
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
    k0vtxmomentum, //|p1 + p2| using daughter pion momenta [GeV/c]
    k0vtxinvariantmass, //Invariant mass from daughter pion hypothesis [GeV/c^2]
    k0vtxmomentumpandora, //|p1 + p2| using Pandora daughter directions [GeV/c]
    k0vtxinvariantmasspandora, //Invariant mass with Pandora directions [GeV/c^2]
    k0vtxmomentumfit, //|p1 + p2| using fit daughter directions [GeV/c]
    k0vtxinvariantmassfit, //Invariant mass with fit directions [GeV/c^2]
    k0vtxresultantmomentumreco, // Reco resultant momentum magnitude at vertex [GeV/c]
    k0vtxresultantmomentumtrue, // True resultant momentum magnitude at vertex [GeV/c]
    // Daughter-level variables (filled in FillNeutralKaonVariables_K0vtxDaughters)
    k0vtxdau1momentumreco, // Reco momentum magnitude of daughter 1 [GeV/c]
    k0vtxdau2momentumreco, // Reco momentum magnitude of daughter 2 [GeV/c]
    k0vtxdau1momentumtrue, // True momentum magnitude of daughter 1 [GeV/c]
    k0vtxdau2momentumtrue, // True momentum magnitude of daughter 2 [GeV/c]
    k0vtxdau1trueendproc, // True end process enum for daughter 1
    k0vtxdau2trueendproc, // True end process enum for daughter 2
    k0vtxdau1trueendmom, // True end momentum for daughter 1 [GeV/c]
    k0vtxdau2trueendmom, // True end momentum for daughter 2 [GeV/c]
    k0vtxdau1truestartmom, // True start momentum for daughter 1 [GeV/c]
    k0vtxdau2truestartmom, // True start momentum for daughter 2 [GeV/c]
    k0vtxdau1truendau, // True number of daughters for daughter 1
    k0vtxdau2truendau, // True number of daughters for daughter 2
    k0vtxdau1truelength, // True length for daughter 1 [cm]
    k0vtxdau2truelength, // True length for daughter 2 [cm]
    k0vtxdau1recolength, // Reco length for daughter 1 [cm]
    k0vtxdau2recolength, // Reco length for daughter 2 [cm]
    k0vtxdau1nhitsreco, // Reco number of collection-plane hits for daughter 1
    k0vtxdau2nhitsreco, // Reco number of collection-plane hits for daughter 2
    k0vtxdau1ndaughtersreco, // Reco number of daughters of vertex daughter 1
    k0vtxdau2ndaughtersreco, // Reco number of daughters of vertex daughter 2
    k0vtxdau1nrecodau, // Number of true daughters of daughter 1 with at least one reco object
    k0vtxdau2nrecodau, // Number of true daughters of daughter 2 with at least one reco object
  enumNeutralKaonMicroTreesLast
  };
}

#endif