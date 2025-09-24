#ifndef neutralKaonTree_h
#define neutralKaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace neutralKaonTree {

  // Methods to add to the output tree the neutralKaonAnalysis sets of variables
  void AddNeutralKaonVariables_VertexCandidates(OutputManager& output, UInt_t nmax);

  // Methods to fill the neutralKaonAnalysis sets of variables in the output tree
  void FillNeutralKaonVariables_VertexCandidates(OutputManager& output, const std::vector<AnaVertexPD*>& reconCandidates, const std::vector<AnaTrueVertexPD*>& trueCandidates); // this should go out

  // Enum with unique indexes for output tree variables
  enum enumNeutralKaonMicroTrees{

    // Vertex candidates info
    n_recovtx_candidates = standardPDTree::enumStandardMicroTreesLast_standardPDTree+1,
    n_true_vertex_candidates,
    recovtx_recoparticles,
    recovtx_trueparticles,
    true_vertex_nparticles,
    recovtx_recoposition,
    recovtx_trueposition,
    true_vertex_position,
    recovtx_recoquality,
    recovtx_truequality,
    true_vertex_quality,
    recovtx_par_truepdg,
    recovtx_par_recopdg,
    recovtx_par_recodir,
    recovtx_par_truedir,
    recovtx_par_dau_recoangles,
    recovtx_par_dau_trueangles,
    recovtx_dau_recodistances_ss,
    recovtx_dau_recodistances_se,
    recovtx_dau_recodistances_es,
    recovtx_dau_recodistances_ee,
    recovtx_dau_recodistances,
    recovtx_dau_truedistances,
    recovtx_par_dau_recodistances,
    recovtx_par_dau_truedistances,
    recovtx_par_dau_recosyst_angle,
    recovtx_par_dau_truesyst_angle,
    recovtx_dau_recolengths,
    recovtx_dau_truelengths,
    recovtx_dau_truepdgs,
    recovtx_dau_recopdgs,
    recovtx_dau_trueproc,
    recovtx_dau_trueendproc,
    recovtx_dau_recoendproc,
    recovtx_dau_par_recoangles,
    recovtx_dau_par_trueangles,
    recovtx_dau_par_truepdg,
    recovtx_dau_recomom,
    recovtx_dau_truemom,
    recovtx_dau_recodedx,
    recovtx_dau_truededx,
    recovtx_recomass,
    recovtx_truemass,
    recovtx_daupi_recodistance,
    recovtx_daupi_truedistance,

    enumNeutralKaonMicroTreesLast
  };
}

#endif