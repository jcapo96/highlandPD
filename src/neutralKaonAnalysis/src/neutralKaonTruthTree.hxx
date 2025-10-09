#ifndef neutralKaonTruthTree_h
#define neutralKaonTruthTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdBaseAnalysis.hxx"
#include "neutralKaonTree.hxx"

namespace neutralKaonTruthTree{

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonTruthVariables(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonParentTruthVariables(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonDaughter1TruthVariables(OutputManager& output, UInt_t nmax);
  void AddNeutralKaonDaughter2TruthVariables(OutputManager& output, UInt_t nmax);
  void FillNeutralKaonTruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject = false);
  void FillNeutralKaonParentTruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject = false);
  void FillNeutralKaonDaughter1TruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject = false);
  void FillNeutralKaonDaughter2TruthVariables(OutputManager& output, const AnaTrueParticlePD& part, bool hasRecoObject = false);
  void FillVertexReconstructionDebugVariables(OutputManager& output, const AnaTrueParticlePD& part,
                                              AnaParticlePD* daughter1Reco, AnaParticlePD* daughter2Reco,
                                              AnaParticlePD* parentReco,
                                              double maxDaughterDistance, double trackFitLength, double vertexRadius);

  // enum with unique indexes for output tree variables
  enum enumNeutralKaonTruthMicroTrees{
    ntruek0 = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,

    // K0 variables
    k0truepdg,
    k0truendau,
    k0trueproc,
    k0trueendproc,
    k0truestartpos,
    k0trueendpos,
    k0truestartdir,
    k0trueenddir,
    k0truestartmom,
    k0trueendmom,
    k0truelength,
    k0truestartenddir,
    k0truegeneration,
    k0hasrecoobject,

    // Vertex reconstruction debugging variables
    k0dau1hasreco,
    k0dau2hasreco,
    k0parenthasreco,
    k0dau1validstartpos,
    k0dau1validendpos,
    k0dau2validstartpos,
    k0dau2validendpos,
    k0daucloseenough,
    k0daudistance,
    k0dau1fitok,
    k0dau2fitok,
    k0vtxpositionfound,
    k0vtxmindistance,
    k0parentrecodist,
    k0parentwithinradius,

    //k0 parent variables
    k0partruepdg,
    k0partruendau,
    k0partrueproc,
    k0partrueendproc,
    k0partruestartpos,
    k0partrueendpos,
    k0partruestartdir,
    k0partrueenddir,
    k0partruestartmom,
    k0partrueendmom,
    k0partruelength,
    k0partruestartenddir,
    k0parisbeam,
    k0partruegeneration,
    k0parhasrecoobject,

    //k0 daughter1 variables
    k0dau1truepdg,
    k0dau1truendau,
    k0dau1trueproc,
    k0dau1trueendproc,
    k0dau1truestartpos,
    k0dau1trueendpos,
    k0dau1truestartdir,
    k0dau1trueenddir,
    k0dau1truestartmom,
    k0dau1trueendmom,
    k0dau1truelength,
    k0dau1truestartenddir,
    k0dau1truegeneration,
    k0dau1hasrecoobject,

    //k0 daughter2 variables
    k0dau2truepdg,
    k0dau2truendau,
    k0dau2trueproc,
    k0dau2trueendproc,
    k0dau2truestartpos,
    k0dau2trueendpos,
    k0dau2truestartdir,
    k0dau2trueenddir,
    k0dau2truestartmom,
    k0dau2trueendmom,
    k0dau2truelength,
    k0dau2truestartenddir,
    k0dau2truegeneration,
    k0dau2hasrecoobject,

    enumNeutralKaonTruthMicroTreesLast_neutralKaonTruthTree
  };
}

#endif