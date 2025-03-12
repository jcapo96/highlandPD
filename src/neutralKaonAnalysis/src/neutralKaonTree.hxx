#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

    void AddNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart);

    void AddNeutralKaonVariables_TrueDaughter1Candidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueDaughter1Candidates(OutputManager& output, const AnaTrueParticlePD* truePart);

    void AddNeutralKaonVariables_TrueDaughter2Candidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueDaughter2Candidates(OutputManager& output, const AnaTrueParticlePD* truePart);

    void AddNeutralKaonVariables_TrueParentCandidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueParentCandidates(OutputManager& output, const AnaTrueParticlePD* truePart);

    void AddNeutralKaonVariables_TrueGrandParentCandidates(OutputManager& output);
    void FillNeutralKaonVariables_TrueGrandParentCandidates(OutputManager& output, const AnaTrueParticlePD* truePart);

    AnaTrueParticleB::ProcessEnum NewFunction(const AnaTrueParticlePD* truePart);

// Enum with unique indexes for output tree neutralKaon variables
enum enumNeutralKaonMicroTrees{

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

    truedau1_truemom,
    truedau1_trueendmom,
    truedau1_trueparentpdg,
    truedau1_trueparentid,
    truedau1_truepdg,
    truedau1_trueproc,
    truedau1_trueendproc,
    truedau1_truedecay,
    truedau1_truechainmuon,
    truedau1_truendau,
    truedau1_truepos,
    truedau1_trueendpos,
    truedau1_truedir,
    truedau1_trueenddir,
    truedau1_truegeneration,
    truedau1_trueeff,
    truedau1_truepur,
    truedau1_branch,

    truedau2_truemom,
    truedau2_trueendmom,
    truedau2_trueparentpdg,
    truedau2_trueparentid,
    truedau2_truepdg,
    truedau2_trueproc,
    truedau2_trueendproc,
    truedau2_truedecay,
    truedau2_truechainmuon,
    truedau2_truendau,
    truedau2_truepos,
    truedau2_trueendpos,
    truedau2_truedir,
    truedau2_trueenddir,
    truedau2_truegeneration,
    truedau2_trueeff,
    truedau2_truepur,
    truedau2_branch,

    truepar_truemom,
    truepar_trueendmom,
    truepar_trueparentpdg,
    truepar_trueparentid,
    truepar_truepdg,
    truepar_trueproc,
    truepar_trueendproc,
    truepar_truedecay,
    truepar_truechainmuon,
    truepar_truendau,
    truepar_truepos,
    truepar_trueendpos,
    truepar_truedir,
    truepar_trueenddir,
    truepar_truegeneration,
    truepar_trueeff,
    truepar_truepur,
    truepar_branch,

    truegpar_truemom,
    truegpar_trueendmom,
    truegpar_trueparentpdg,
    truegpar_trueparentid,
    truegpar_truepdg,
    truegpar_trueproc,
    truegpar_trueendproc,
    truegpar_truedecay,
    truegpar_truechainmuon,
    truegpar_truendau,
    truegpar_truepos,
    truegpar_trueendpos,
    truegpar_truedir,
    truegpar_trueenddir,
    truegpar_truegeneration,
    truegpar_trueeff,
    truegpar_truepur,
    truegpar_branch,
    };

}

#endif