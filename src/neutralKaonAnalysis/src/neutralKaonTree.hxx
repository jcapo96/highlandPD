#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_VertexCandidates(OutputManager& output, UInt_t nmax);

  // Methods to fill the neutralKaonAnalysis sets of variables in the output tree
  void FillNeutralKaonVariables_VertexCandidates(OutputManager& output, const std::vector<AnaVertexPD*>& reconCandidates, const std::vector<AnaTrueVertexPD*>& trueCandidates, const AnaEventB& event); // this should go out

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Vertex candidates info
    nvcandidates = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    ntruevertexcandidates, //number of true vertex candidates
    vrecoparticles, //reconstructed vertex reconstructed number of particles
    vtrueparticles, //reconstructed vertex true number of particles
    truevertexnparticles, //true vertex true number of particles
    vrecoposition, //reconstructed vertex reconstructed position
    vtrueposition, //true vertex true position
    truevertexposition, //true vertex true position
    vrecoquality, //reconstructed vertex reconstructed quality
    vtruequality, //reconstructed vertex true quality
    truevertexquality, //true vertex true quality
    vparrecopdg, //reconstructed vertex parent reconstructed PDG
    vpartruepdg, //reconstructed vertex parent true PDG
    vparrecodir, //reconstructed vertex parent reconstructed direction
    vpartruedir, //reconstructed vertex parent true direction
    vpardaurecoang, //reconstructed vertex: reconstructed angles between parent end and daughter start positions
    vpardautrueang, //reconstructed vertex: true angles between parent end and daughter start positions
    vdaurecodistss, //reconstructed vertex: reconstructed distances between daughter daughters (start-start)
    vdaurecodistse, //reconstructed vertex: reconstructed distances between daughter daughters (start-end)
    vdaurecodistes, //reconstructed vertex: reconstructed distances between daughter daughters (end-start)
    vdaurecodistee, //reconstructed vertex: reconstructed distances between daughter daughters (end-end)
    vdaurecodist, //reconstructed vertex: reconstructed distances between daughter daughters (min of 4 combinations)
    vdautruedist, //reconstructed vertex: true distances between daughter daughters (min of 4 combinations)
    vpardaurecodist, //reconstructed vertex: reconstructed distances between parent end and daughter start positions
    vpardautruedist, //reconstructed vertex: true distances between parent end and daughter start positions
    vpardaurecosystang, //reconstructed vertex: reconstructed system angles between parent end and daughter start positions
    vpardautruesystang, //reconstructed vertex: true system angles between parent end and daughter start positions
    vdaurecolength,
    vdautruelength,
    vdaurecopdg,
    vdautruepdg,
    vdautrueproc,
    vdaurecoproc,
    vdaurecoendproc,
    vdautrueendproc,
    vdauparrecoang,
    vdaupartrueang,
    vdaupartruepdg,
    vdaurecomom,
    vdautruemom,
    vdaurecodedx,
    vdautruededx,
    vrecomass,
    vtruemass,
    vdaupirecodist,
    vdaupitruedist,

    // New variables for parent daughters
    vparrecondau,
    vpartruendau,
    vpardaurecopdgs,
    vpardautruepdgs,
    vpardaurecoproc,
    vpardautrueproc,
    vpardaurecoendproc,
    vpardautrueendproc,

    // New variables for daughter daughters
    vdaurecondau,
    vdautruendau,
    vdaudaurecopdgs,
    vdaudautruepdgs,
    vdaudaurecoproc,
    vdaudautrueproc,
    vdaudaurecoendproc,
    vdaudautrueendproc,

    enumNeutralKaonMicroTreesLast
  };
}

#endif