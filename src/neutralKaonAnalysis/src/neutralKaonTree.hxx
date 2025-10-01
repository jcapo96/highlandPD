#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_Candidates(OutputManager& output, UInt_t nmax);

  // Methods to fill the neutralKaonAnalysis sets of variables in the output tree
  void FillNeutralKaonVariables_Candidates(OutputManager& output, const std::vector<AnaNeutralParticlePD*>& candidates, const AnaEventB& event);
  void FillNeutralKaonVariables_SingleCandidate(OutputManager& output, AnaNeutralParticlePD* candidate, const AnaEventB& event);

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Candidate info
    nk0 = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,

    // Variables about K0
    k0id,
    k0recostartpos,
    k0truestartpos,
    k0recostartdir,
    k0truestartdir,
    k0recoendpos,
    k0trueendpos,
    k0recoenddir,
    k0trueenddir,
    k0truestartenddir,
    k0recostartenddir,
    k0recolength,
    k0truelength,
    k0truestartmom,
    k0trueendmom,
    k0recomass,
    k0truemass,
    k0truendau,
    k0truenbrothers,
    k0truebrotherspdg,
    k0trueproc,
    k0trueendproc,
    k0hastrueobject,
    k0truepdg,

    k0truerecodist,
    k0impactparameter,
    k0vtxoriginaldistance,
    k0vtxminimumdistance,

    // Variables about K0 daughters
    k0daurecostartpos,
    k0dautruestartpos,
    k0daurecoendpos,
    k0dautrueendpos,
    k0daurecostartdir,
    k0dautruestartdir,
    k0daurecoenddir,
    k0dautrueenddir,
    k0daurecolength,
    k0dautruelength,
    k0dautruestartmom,
    k0dautrueendmom,
    k0dautruendau,
    k0dautrueproc,
    k0dautrueendproc,
    k0dautruepdg,

    // Variables about K0 parent
    k0parrecostartpos,
    k0partruestartpos,
    k0parrecoendpos,
    k0partrueendpos,
    k0parrecostartdir,
    k0partruestartdir,
    k0parrecoenddir,
    k0partrueenddir,
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
    k0beaminstpdg,
    k0partrueproc,
    k0partrueendproc,

    // Variables about the vertex system (two particles)
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
    k0vtxrecoinvariantmass,
    k0vtxtrueinvariantmass,
    k0vtxrecoctau,
    k0vtxtruectau,
    k0vtxrecoangle,
    k0vtxtrueangle,
    k0vtxnpotpar,

    // Variables about the vertex pions only
    k0vtxpplength,
    k0vtxpmlength,
    k0vtxpptruepdg,
    k0vtxpmtruepdg,
    k0vtxpptruendau,
    k0vtxpmtruendau,
    k0vtxpprecondau,
    k0vtxpmrecondau,
    k0vtxpptruepdgdau,
    k0vtxpmtruepdgdau,
    k0vtxpprecopdgdau,
    k0vtxpmrecopdgdau,
    k0vtxpptruestartdir,
    k0vtxpmtruestartdir,
    k0vtxpptrueenddir,
    k0vtxpmtrueenddir,
    k0vtxpprecostartdir,
    k0vtxpmrecostartdir,
    k0vtxpprecoenddir,
    k0vtxpmrecoenddir,
    k0vtxpptrueproc,
    k0vtxpmtrueproc,
    k0vtxpptrueendproc,
    k0vtxpmtrueendproc,
    k0vtxtruepptruedaupdg,
    k0vtxtruepmtruedaupdg,

    // Brother variables (proton daughters of K0 parent)
    k0brothnproton,
    k0brothtruepdg,
    k0brothrecopdg,
    k0brothtrueproc,
    k0brothtrueendproc,
    k0brothtruendau,
    k0brothrecondau,
    k0brothlength,

  enumNeutralKaonMicroTreesLast
  };
}

#endif