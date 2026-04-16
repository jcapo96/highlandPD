#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_K0Particle(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Parent(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0CreationVtx(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0VtxDaughters(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Kinematics(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0AlignVariants(OutputManager& output, UInt_t nmax);

  void FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                Int_t nVerticesBeforeFiltering = -1, Int_t nVerticesAfterFiltering = -1,
                                AnaBeamB* beam = NULL, size_t neutralCandidateIndex = 0);
  void FillNeutralKaonVariables_K0Particle(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0Parent(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0CreationVtx(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
  void FillNeutralKaonVariables_K0vtxDaughters(OutputManager& output, AnaAnnihilationVertexPD* vertex, const AnaEventB& event);
  void FillNeutralKaonVariables_K0Kinematics(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0AlignVariants(OutputManager& output, AnaNeutralParticlePD* candidate);
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
    k0al_alltrue, // Alignment with all true quantities
    k0al_allreco, // Alignment with all reco quantities (fit reco geometry/directions)
    k0al_cvreco, // Alignment with reco creation only (all else true)
    k0al_avreco, // Alignment with reco annihilation only (all else true)
    k0al_vtxreco, // Alignment with reco creation+annihilation only (all else true)
    k0al_d1magreco, // Alignment with daughter1 reco momentum magnitude only
    k0al_d2magreco, // Alignment with daughter2 reco momentum magnitude only
    k0al_d1dirreco, // Alignment with daughter1 reco momentum direction only
    k0al_d2dirreco, // Alignment with daughter2 reco momentum direction only
    k0al_d1preco, // Alignment with daughter1 reco momentum (direction+magnitude)
    k0al_d2preco, // Alignment with daughter2 reco momentum (direction+magnitude)
    k0al_d12magreco, // Alignment with both daughters reco momentum magnitudes only
    k0al_d12dirreco, // Alignment with both daughters reco momentum directions only
    k0al_d12preco, // Alignment with both daughters reco momentum (direction+magnitude)
    k0al_allreco_cvtrue, // Alignment with all reco except creation vertex true
    k0al_allreco_avtrue, // Alignment with all reco except annihilation vertex true
    k0al_allreco_d1true, // Alignment with all reco except daughter1 momentum true
    k0al_allreco_d2true, // Alignment with all reco except daughter2 momentum true
    // K0 kinematics: true K0 positions, momentum, direction
    k0truecreationx, // True K0 start position x [cm]
    k0truecreationy, // True K0 start position y [cm]
    k0truecreationz, // True K0 start position z [cm]
    k0trueannihilationx, // True K0 end position (decay) x [cm]
    k0trueannihilationy, // True K0 end position (decay) y [cm]
    k0trueannihilationz, // True K0 end position (decay) z [cm]
    k0truemomentum, // True K0 momentum magnitude [GeV/c]
    k0truedirectionx, // True K0 direction (normalized) x component
    k0truedirectiony, // True K0 direction (normalized) y component
    k0truedirectionz, // True K0 direction (normalized) z component
    // Reco K0 creation position: Pandora variant (raw parent end)
    k0creationpandorax, // Reco K0 creation Pandora x from raw parent end [cm]
    k0creationpandoray, // Reco K0 creation Pandora y from raw parent end [cm]
    k0creationpandoraz, // Reco K0 creation Pandora z from raw parent end [cm]
    // Reco K0 creation position: Fit variant (projected parent end)
    k0creationfitx, // Reco K0 creation Fit x from projected parent end [cm]
    k0creationfity, // Reco K0 creation Fit y from projected parent end [cm]
    k0creationfitz, // Reco K0 creation Fit z from projected parent end [cm]
    // Reco K0 direction: Pandora variant
    k0directionpandorax, // Reco K0 direction Pandora (normalized) x component
    k0directionpandoray, // Reco K0 direction Pandora (normalized) y component
    k0directionpandoraz, // Reco K0 direction Pandora (normalized) z component
    // Reco K0 direction: Fit variant
    k0directionfitx, // Reco K0 direction Fit (normalized) x component
    k0directionfity, // Reco K0 direction Fit (normalized) y component
    k0directionfitz, // Reco K0 direction Fit (normalized) z component
    // True daughter directions at annihilation vertex
    k0vtxdau1truedirectionx, // True daughter 1 direction (normalized) x component
    k0vtxdau1truedirectiony, // True daughter 1 direction (normalized) y component
    k0vtxdau1truedirectionz, // True daughter 1 direction (normalized) z component
    k0vtxdau2truedirectionx, // True daughter 2 direction (normalized) x component
    k0vtxdau2truedirectiony, // True daughter 2 direction (normalized) y component
    k0vtxdau2truedirectionz, // True daughter 2 direction (normalized) z component
    // Reco daughter directions: Pandora variant
    k0vtxdau1directionpandorax, // Reco daughter 1 Pandora direction x component
    k0vtxdau1directionpandoray, // Reco daughter 1 Pandora direction y component
    k0vtxdau1directionpandoraz, // Reco daughter 1 Pandora direction z component
    k0vtxdau2directionpandorax, // Reco daughter 2 Pandora direction x component
    k0vtxdau2directionpandoray, // Reco daughter 2 Pandora direction y component
    k0vtxdau2directionpandoraz, // Reco daughter 2 Pandora direction z component
    // Reco daughter directions: Fit variant
    k0vtxdau1directionfitx, // Reco daughter 1 Fit direction x component
    k0vtxdau1directionfity, // Reco daughter 1 Fit direction y component
    k0vtxdau1directionfitz, // Reco daughter 1 Fit direction z component
    k0vtxdau2directionfitx, // Reco daughter 2 Fit direction x component
    k0vtxdau2directionfity, // Reco daughter 2 Fit direction y component
    k0vtxdau2directionfitz, // Reco daughter 2 Fit direction z component
    // Parent-level variables (parent of reconstructed neutral candidate)
    k0partruepdg, // True PDG of neutral parent
    k0partrueendmom, // True end momentum of neutral parent [GeV/c]
    k0partruelength, // True track length of neutral parent [cm]
    k0parrecolength, // Reco track length of neutral parent [cm]
    k0parndau, // Reco number of daughters of neutral parent
    k0parisbeam, // 1 if reco parent has isPandora==true (beam), 0 otherwise
    // Creation-vertex residuals: reco - true (true creation from true K0 start when available)
    k0cvtxpandoraresidual, // |creation Pandora(raw parent end) - true creation| [cm]
    k0cvtxfitresidual, // |creation Fit(projected parent end) - true creation| [cm]
    k0cvtxpandoraresidualx, // creation Pandora x_reco - x_true [cm]
    k0cvtxpandoraresidualy, // creation Pandora y_reco - y_true [cm]
    k0cvtxpandoraresidualz, // creation Pandora z_reco - z_true [cm]
    k0cvtxfitresidualx, // creation Fit x_reco - x_true [cm]
    k0cvtxfitresidualy, // creation Fit y_reco - y_true [cm]
    k0cvtxfitresidualz, // creation Fit z_reco - z_true [cm]
    // Active microtree variables (vertex-only set)
    k0vtxtruepos, //K0 vertex true position
    k0vtxoriginaldistance, //K0 vertex original distance
    k0vtxtrueoriginaldistance, //K0 vertex true original distance
    k0vtxdegeneracy, // K0 annihilation-vertex degeneracy count within the configured degeneracy radius
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
    k0vtxopeninganglepandora, // Opening angle between daughter directions using Pandora dirs [rad]
    k0vtxmomentumfit, //|p1 + p2| using fit daughter directions [GeV/c]
    k0vtxinvariantmassfit, //Invariant mass with fit directions [GeV/c^2]
    k0vtxopeninganglefit, // Opening angle between daughter directions using fit dirs [rad]
    k0vtxresultantmomentumreco, // Reco resultant momentum magnitude at vertex [GeV/c]
    k0vtxresultantmomentumtrue, // True resultant momentum magnitude at vertex [GeV/c]
    // Daughter-level variables (filled in FillNeutralKaonVariables_K0vtxDaughters)
    k0vtxdau1momentumreco, // Reco momentum magnitude of daughter 1 [GeV/c]
    k0vtxdau2momentumreco, // Reco momentum magnitude of daughter 2 [GeV/c]
    k0vtxdau1mommethod, // Momentum assignment method enum for daughter 1
    k0vtxdau2mommethod, // Momentum assignment method enum for daughter 2
    k0vtxdau1extchi2ndf, // Free-range fit log-likelihood used by momentum assignment for daughter 1
    k0vtxdau2extchi2ndf, // Free-range fit log-likelihood used by momentum assignment for daughter 2
    k0vtxdau1dedxdrift, // Mean dEdx bias (actual - expected from fit PDF mode) for daughter 1 [MeV/cm]
    k0vtxdau2dedxdrift, // Mean dEdx bias (actual - expected from fit PDF mode) for daughter 2 [MeV/cm]
    k0vtxdau1dedxsigma, // Sigma of dEdx bias Gaussian fit for daughter 1 [MeV/cm]
    k0vtxdau2dedxsigma, // Sigma of dEdx bias Gaussian fit for daughter 2 [MeV/cm]
    k0vtxdau1dedxfitok, // 1 if Gaussian dEdx bias fit succeeded for daughter 1, 0 otherwise
    k0vtxdau2dedxfitok, // 1 if Gaussian dEdx bias fit succeeded for daughter 2, 0 otherwise
    k0vtxdau1momentumtrue, // True momentum magnitude of daughter 1 [GeV/c]
    k0vtxdau2momentumtrue, // True momentum magnitude of daughter 2 [GeV/c]
    k0vtxdau1trueendproc, // True end process enum for daughter 1
    k0vtxdau2trueendproc, // True end process enum for daughter 2
    k0vtxdau1trueendmom, // True end momentum for daughter 1 [GeV/c]
    k0vtxdau2trueendmom, // True end momentum for daughter 2 [GeV/c]
    k0vtxdau1truestartmom, // True start momentum for daughter 1 [GeV/c]
    k0vtxdau2truestartmom, // True start momentum for daughter 2 [GeV/c]
    k0vtxdau1truestartx, // True start x for daughter 1 [cm]
    k0vtxdau1truestarty, // True start y for daughter 1 [cm]
    k0vtxdau1truestartz, // True start z for daughter 1 [cm]
    k0vtxdau2truestartx, // True start x for daughter 2 [cm]
    k0vtxdau2truestarty, // True start y for daughter 2 [cm]
    k0vtxdau2truestartz, // True start z for daughter 2 [cm]
    k0vtxdau1truendau, // True number of daughters for daughter 1
    k0vtxdau2truendau, // True number of daughters for daughter 2
    k0vtxdau1truelength, // True length for daughter 1 [cm]
    k0vtxdau2truelength, // True length for daughter 2 [cm]
    k0vtxdau1recolength, // Reco length for daughter 1 [cm]
    k0vtxdau2recolength, // Reco length for daughter 2 [cm]
    k0vtxdau1recostartx, // Reco start x for daughter 1 [cm]
    k0vtxdau1recostarty, // Reco start y for daughter 1 [cm]
    k0vtxdau1recostartz, // Reco start z for daughter 1 [cm]
    k0vtxdau2recostartx, // Reco start x for daughter 2 [cm]
    k0vtxdau2recostarty, // Reco start y for daughter 2 [cm]
    k0vtxdau2recostartz, // Reco start z for daughter 2 [cm]
    k0vtxdau1nhitsreco, // Reco number of collection-plane hits for daughter 1
    k0vtxdau2nhitsreco, // Reco number of collection-plane hits for daughter 2
    k0vtxdau1protonchi2ndf, // Reco proton PID chi2/ndf for daughter 1
    k0vtxdau2protonchi2ndf, // Reco proton PID chi2/ndf for daughter 2
    k0vtxdau1pionchi2ndf, // Reco pion PID chi2/ndf for daughter 1
    k0vtxdau2pionchi2ndf, // Reco pion PID chi2/ndf for daughter 2
    k0vtxdau1truepdg, // True PDG of daughter 1
    k0vtxdau2truepdg, // True PDG of daughter 2
    k0vtxdau1ndaughtersreco, // Reco number of daughters of vertex daughter 1
    k0vtxdau2ndaughtersreco, // Reco number of daughters of vertex daughter 2
    k0vtxdau1nrecodau, // Number of true daughters of daughter 1 with at least one reco object
    k0vtxdau2nrecodau, // Number of true daughters of daughter 2 with at least one reco object
  enumNeutralKaonMicroTreesLast
  };
}

#endif