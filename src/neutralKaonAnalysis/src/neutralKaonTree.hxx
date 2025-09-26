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
  void FillNeutralKaonVariables_SingleCandidate(OutputManager& output, AnaNeutralParticlePD* candidate);

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Candidate info
    nk0 = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    k0id,
    k0length,
    k0recoendpos,
    k0trueendpos,
    k0dautruepdg,

    enumNeutralKaonMicroTreesLast
  };
}

#endif