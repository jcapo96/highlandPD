#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_K0(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Par(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0vtxDaughter1(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0vtxDaughter2(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Brother(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonVariables_K0Vtx(OutputManager& output, UInt_t nmax);
  // void AddNeutralKaonVariables_K0VtxPions(OutputManager& output, UInt_t nmax);

  void FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event);
  void FillNeutralKaonVariables_K0(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0Par(OutputManager& output, AnaParticlePD* parentCandidate);
  void FillNeutralKaonVariables_K0Brother(OutputManager& output, AnaParticlePD* parentCandidate);
  void FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaVertexPD* vertex);
  void FillNeutralKaonVariables_K0vtxDaughter1(OutputManager& output, AnaParticlePD* daughterCandidate);
  void FillNeutralKaonVariables_K0vtxDaughter2(OutputManager& output, AnaParticlePD* daughterCandidate);

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Candidate info
    nk0 = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    // Variables about K0
    k0id, //unique ID of neutral particle candidates
    k0recostartpos, //reconstructed start position: from parent end to vertex start
    k0truestartpos, //true start position: k0hastrueobject=1: true object start, k0hastrueobject=0: true parent end
    k0recostartdir, //reconstructed start direction
    k0truestartdir, //true start direction: k0hastrueobject=1: true object start direction, k0hastrueobject=0: vector from true parent end to true vertex start
    k0recoendpos, //reconstructed end position
    k0trueendpos, //true end position: k0hastrueobject=1: true object end, k0hastrueobject=0: true vertex start
    k0recoenddir, //reconstructed end direction
    k0trueenddir, //true end direction: k0hastrueobject=1: true object end direction, k0hastrueobject=0: true vertex start direction
    k0truestartenddir, //true start-end direction scalar product
    k0recostartenddir, //reconstructed start-end direction scalar product
    k0recolength, //reconstructed length
    k0truelength, //true length
    k0truestartmom, //true start momentum
    k0trueendmom, //true end momentum
    k0recomass, //reconstructed mass from vertex particles under pion hypothesis
    k0truemass, //true mass from vertex particles under pion hypothesis
    k0truendau, //true number of daughters
    k0trueproc, //true creation process
    k0trueendproc, // true annihilation process
    k0hastrueobject, //1: if we have a true object, 0: if we don't have a true object
    k0hasrecoparticle, //1: if true object has a reconstructed particle, 0: if true object doesn't have a reconstructed particle
    k0hasequivalenttrueobject, //1: if true object has an equivalent true object, 0: if true object doesn't have an equivalent true object (should be always 1)
    k0truepdg, //true PDG
    k0recopdg, //reconstructed PDG
    k0truegeneration, //true generation
    k0nrecohits, //number of reconstructed hits in cylinder of radius X from parent end to vertex start

    k0truerecodist, //True-vertex-position minus reconstructed-vertex-position
    k0impactparameter, //Impact parameter

    // Variables about K0 daughter1
    k0dau1recostartpos, //reconstructed start position
    k0dau1truestartpos, //true start position
    k0dau1recoendpos, //reconstructed end position
    k0dau1trueendpos, //true end position
    k0dau1recostartdir,
    k0dau1truestartdir,
    k0dau1recomom,
    k0dau1recoenddir,
    k0dau1trueenddir,
    k0dau1recolength,
    k0dau1truelength,
    k0dau1truestartmom,
    k0dau1trueendmom,
    k0dau1recondau,
    k0dau1truendau,
    k0dau1trueproc,
    k0dau1trueendproc,
    k0dau1truepdg,
    k0dau1parid,
    k0dau1istrack, //0: Unknown, 1: Shower, 2: Track
    k0dau1chi2pion, //chi2 for pion PID hypothesis
    k0dau1nptchi2pion, //number of points used in chi2 pion calculation
    k0dau1avgdedx15cm, //average dEdx for residual range < 15cm

    //Variables about K0 daughter2
    k0dau2recostartpos,
    k0dau2truestartpos,
    k0dau2recoendpos,
    k0dau2trueendpos,
    k0dau2recostartdir,
    k0dau2truestartdir,
    k0dau2recomom,
    k0dau2recoenddir,
    k0dau2trueenddir,
    k0dau2recolength,
    k0dau2truelength,
    k0dau2truestartmom,
    k0dau2trueendmom,
    k0dau2recondau,
    k0dau2truendau,
    k0dau2trueproc,
    k0dau2trueendproc,
    k0dau2truepdg,
    k0dau2parid,
    k0dau2istrack, //0: Unknown, 1: Shower, 2: Track
    k0dau2chi2pion, //chi2 for pion PID hypothesis
    k0dau2nptchi2pion, //number of points used in chi2 pion calculation
    k0dau2avgdedx15cm, //average dEdx for residual range < 15cm

    // Variables about K0 parent
    k0parrecostartpos,
    k0partruestartpos,
    k0parrecoendpos,
    k0partrueendpos,
    k0parrecostartdir,
    k0partruestartdir,
    k0parrecoenddir,
    k0partrueenddir,
    k0parrecostartenddir,
    k0partruestartenddir,
    k0parrecolength,
    k0partruelength,
    k0partruestartmom,
    k0partrueendmom,
    k0partruendau,
    k0parrecondau,
    k0partruepdg,
    k0parrecopdg,
    k0partruepdgdau,
    k0parrecopdgdau,
    k0parisbeam,
    k0partrueproc,
    k0partrueendproc,
    k0partruegeneration,

    // Variables about K0 brothers (siblings)
    k0reconbrother,
    k0truenbrother,
    k0brothreconprot,
    k0brothreconpiplus,
    k0brothreconpiminus,

     // Variables about the vertex system (two particles)
    k0vtxrecopos,
    k0vtxtruepos,
    k0vtxrecomom,
    k0vtxtruemom,
    k0vtxrecoenergy,
    k0vtxtrueenergy,
    k0vtxrecomass,
    k0vtxtruemass,
    k0vtxrecodir,
    k0vtxtruedir,
    k0vtxrecoopening,
    k0vtxtrueopening,
    k0vtxrecoangle,
    k0vtxtrueangle,
    k0vtxnpotpar,
    k0vtxoriginaldistance,
    k0vtxminimumdistance,
    k0vtxtrueoriginaldistance,
    k0vtxtrueminimumdistance,

  enumNeutralKaonMicroTreesLast
  };
}

#endif