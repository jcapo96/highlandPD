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

  void FillNeutralKaonVariables(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event, AnaBeamB* beam = NULL);
  void FillNeutralKaonVariables_K0(OutputManager& output, AnaNeutralParticlePD* candidate);
  void FillNeutralKaonVariables_K0Par(OutputManager& output, AnaNeutralParticlePD* neutralCandidate, const AnaEventB& event, AnaBeamB* beam = NULL);
  void FillNeutralKaonVariables_K0Brother(OutputManager& output, AnaNeutralParticlePD* neutralCandidate, AnaParticlePD* parentCandidate, const AnaEventB& event);
  void FillNeutralKaonVariables_K0vtx(OutputManager& output, AnaVertexPD* vertex);
  void FillNeutralKaonVariables_K0vtxDaughter1(OutputManager& output, AnaParticlePD* daughterCandidate, AnaVertexPD* vertex);
  void FillNeutralKaonVariables_K0vtxDaughter2(OutputManager& output, AnaParticlePD* daughterCandidate, AnaVertexPD* vertex);

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
    k0neutralscore, //neutral particle quality score (lower = more neutral-like)
    k0hitsalignment, //alignment between hits in cylinder and neutral particle direction
    k0nhitsincylinder, //number of hits in the cylinder around neutral particle
    k0hitsavgdistance, //average perpendicular distance of hits to neutral particle path
    k0hitsrmsdistance, //RMS of perpendicular distances of hits
    k0hitslongspan, //longitudinal span fraction (span along path / total length)

    k0truerecodist, //True-vertex-position minus reconstructed-vertex-position
    k0impactparameter, //Impact parameter
    k0nprotonincreationvtx, //Number of parent's daughters that are true protons near neutral start position
    k0nparticlesincreationvtx, //Total number of parent's daughters near neutral start position
    k0creationvtxchi2proton, //Chi2/ndf under proton hypothesis for particles near creation vertex (5 closest)
    k0creationvtxdistances, //Pandora-based distances to neutral start (5 closest)
    k0creationvtxtruepdg, //True PDG codes of particles near creation vertex (5 closest)

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
    k0dau1chi2proton, //chi2 for proton PID hypothesis
    k0dau1nptchi2proton, //number of points used in chi2 proton calculation
    k0dau1chi2kaon, //chi2 for kaon PID hypothesis
    k0dau1nptchi2kaon, //number of points used in chi2 kaon calculation
    k0dau1avgdedx5cm, //average dEdx for residual range < 5cm
    k0dau1avgdedx10cm, //average dEdx for residual range < 10cm
    k0dau1avgdedx15cm, //average dEdx for residual range < 15cm
    k0dau1avgdedx25cm, //average dEdx for residual range < 25cm
    k0dau1avgdedx50cm, //average dEdx for residual range < 50cm
    k0dau1nhits, //number of hits in collection plane

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
    k0dau2chi2proton, //chi2 for proton PID hypothesis
    k0dau2nptchi2proton, //number of points used in chi2 proton calculation
    k0dau2chi2kaon, //chi2 for kaon PID hypothesis
    k0dau2nptchi2kaon, //number of points used in chi2 kaon calculation
    k0dau2avgdedx5cm, //average dEdx for residual range < 5cm
    k0dau2avgdedx10cm, //average dEdx for residual range < 10cm
    k0dau2avgdedx15cm, //average dEdx for residual range < 15cm
    k0dau2avgdedx25cm, //average dEdx for residual range < 25cm
    k0dau2avgdedx50cm, //average dEdx for residual range < 50cm
    k0dau2nhits, //number of hits in collection plane

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
    k0partruestartenergy,
    k0partrueendenergy,
    k0partruendau,
    k0parrecondau,
    k0partruepdg,
    k0parrecopdg,
    k0partruepdgdau,
    k0parrecopdgdau,

    //Variables about K0 parent beam
    k0parisbeam,
    k0beaminstpdg,
    k0partrueproc,
    k0partrueendproc,
    k0partruegeneration,

    // Variables about K0 brothers (siblings)
    // Summary counters (kept for backward compatibility)
    k0reconbrother, // Number of reconstructed candidate brothers reconstructed
    k0truenbrother, // Number of true candidate brothers reconstructed
    k0brothreconprot, // Number of reconstructed candidate brothers compatible with proton
    k0brothrecoprotenergy, // Energy of reconstructed candidate brothers compatible with proton with the highest momentum
    k0brothrecoprotmom, // Momentum of reconstructed candidate brothers compatible with proton with the highest momentum
    k0brothreconpiplus, // Number of reconstructed candidate brothers compatible with pi+
    k0brothreconpiminus, // Number of reconstructed candidate brothers compatible with pi-
    k0brothreconprotchi2, //number of brothers compatible with proton (chi2/ndf<60)

    // True brothers - ALL from trueParent->Daughters (including K0)
    k0ntruebroth,
    k0truebrothpdg,
    k0truebrothmomentum,
    k0truebrothenergy,
    k0truebrothprocessstart,
    k0truebrothlength,
    k0truebrothtotalmom,
    k0truebrothtotalenergy,
    k0truebrothtruetotaldir,
    k0truebrothprotonmaxenergy,
    k0truebrothprotonmaxmomentum,
    k0truebrothprotonmaxdir,
    k0truebrothalign,

    // True brothers with reco - subset with reconstructed objects
    k0ntruebrothreco,
    k0truebrothrecotruepdg,
    k0truebrothrecotruemom,
    k0truebrothrecotrueenergy,
    k0truebrothrecotrueprocessstart,
    k0truebrothrecopdg,
    k0truebrothrecomom,
    k0truebrothrecoenergy,
    k0truebrothrecolength,
    k0truebrothrecochi2prot,
    k0truebrothrecotruetotalmom,
    k0truebrothrecotruetotalenergy,
    k0truebrothrecotruetotaldir,
    k0truebrothrecoprotonmaxenergy,
    k0truebrothrecoprotonmaxmomentum,
    k0truebrothrecoprotonmaxrecoenergy,
    k0truebrothrecoprotonmaxrecomom,
    k0truebrothrecoprotonmaxdir,
    k0truebrothrecoalign,

    // Reco brothers - ALL from Parent->Daughters
    k0nrecobroth,
    k0recobrothpdg,
    k0recobrothmom,
    k0recobrothenergy,
    k0recobrothlength,
    k0recobrothchi2prot,
    k0recobrothtruepdg,
    k0recobrothtrueenergy,
    k0recobrothtruemom,
    k0recobrothtrueprocessstart,
    k0recobrothtruetotaldir,
    k0recobrothrecototaldir,
    k0recobrothprotonmaxenergy,
    k0recobrothprotonmaxmomentum,
    k0recobrothprotonmaxtrueenergy,
    k0recobrothprotonmaxtruemom,
    k0recobrothprotonmaxdir,
    k0recobrothprotonmaxtruedir,
    k0recobrothalign,
    k0recobrothtruealign,

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
    k0vtxscore,

    k0vtxpandorapos, //Pandora-based vertex position
    k0vtxfitpos, //Algorithm-specific fitted vertex position (geometric/TMinuit/Kalman)
    k0vtxfitdir, //Algorithm-specific fitted vertex direction (geometric/TMinuit/Kalman)
    k0vtxisjustavg, //Flag: 1 if Pandora used simple average, 0 if line intersection, -999 if invalid
    k0vtxdegbefore, //Degeneracy count before scoring
    k0vtxdegafter, //Degeneracy count after scoring
    k0vtxnrecopart, //Number of unique reco particles in degenerate vertices
    k0vtxdegdistances, //5 minimum distances to vertices within degeneracy threshold
    k0vtxisolationdistances, //5 minimum line-to-point isolation distances (Pandora-based)
    k0vtxisolationdistancesfit, //5 minimum line-to-point isolation distances (fitted track)
    k0vtxisolationstartdistances, //5 minimum point-to-point distances to particle PositionStart
    k0vtxisolnproton, //Total number of isolation particles that are true protons
    k0vtxisolnpion, //Total number of isolation particles that are true pions (both charges)
    k0vtxisolisproton, //Flags (1=proton, 0=not) for 5 closest isolation particles
    k0vtxisolchi2proton, //Chi2/ndf under proton hypothesis for 5 closest isolation particles
    k0vtxisollength, //Lengths of 5 closest isolation particles
    k0vtxisolislongest, //Flag: 1 if any isolation particle longer than both vertex particles

  enumNeutralKaonMicroTreesLast
  };
}

#endif