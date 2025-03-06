#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

    void AddNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart);

// Enum with unique indexes for output tree variables
enum enumKaonMicroTrees{

    //true kaon candidates info
    truenkaons = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    truekaon_truemom,
    truekaon_trueendmom,
    truekaon_trueparentpdg,
    truekaon_trueparentid,
    truekaon_truepdg,
    truekaon_trueproc,
    truekaon_trueendproc,
    truekaon_truedecay,
    truekaon_truechainmuon,
    truekaon_truendau,
    truekaon_truepos,
    truekaon_trueendpos,
    truekaon_truedir,
    truekaon_trueenddir,
    truekaon_truegeneration,
    truekaon_trueeff,
    truekaon_truepur,
    truekaon_branch,
    };
}

#endif