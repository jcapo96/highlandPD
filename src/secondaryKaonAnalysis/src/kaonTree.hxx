#ifndef kaonTree_h
#define kaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "standardPDTree.hxx"

namespace kaonTree{

  // Methods to add to the output tree the kaonAnalysis sets of variables
  void AddKaonVariables_TrueKaonCandidates(OutputManager& output);

  void AddKaonVariables_KaonCandidatesReco(OutputManager& output, UInt_t nmax);
  void AddKaonVariables_KaonCandidatesHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void AddKaonVariables_KaonCandidatesTrue(OutputManager& output, UInt_t nmax);

  void AddKaonVariables_KaonBestCandidateReco(OutputManager& output);
  void AddKaonVariables_KaonBestCandidateHitsReco(OutputManager& output, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void AddKaonVariables_KaonBestCandidateTrue(OutputManager& output);

  // Methods to fill the kaonAnalysis sets of variables in the output tree
  void FillKaonVariables_TrueKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart);

  void FillKaonVariables_KaonCandidatesReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent = NULL);
  void FillKaonVariables_KaonCandidatesHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void FillKaonVariables_KaonCandidatesTrue(OutputManager& output, AnaParticlePD* part);

  void FillKaonVariables_KaonBestCandidateReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent = NULL);
  void FillKaonVariables_KaonBestCandidateHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void FillKaonVariables_KaonBestCandidateTrue(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent = NULL);

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

    truekaon_truemuon_truemom,
    truekaon_truemuon_trueendmom,
    truekaon_truemuon_truepdg,
    truekaon_truemuon_trueproc,
    truekaon_truemuon_trueendproc,
    truekaon_truemuon_truendau,
    truekaon_truemuon_truepos,
    truekaon_truemuon_trueendpos,
    truekaon_truemuon_truedir,
    truekaon_truemuon_trueenddir,
    truekaon_truemuon_truegeneration,
    truekaon_truemuon_trueeff,
    truekaon_truemuon_truepur,

    //kaon candidates info
    ncandidates,
    candidates_generation,
    candidates_parentID,
    candidates_ndau,
    candidates_nsisters,
    candidates_pos,
    candidates_dir,
    candidates_endpos,
    candidates_enddir,
    candidates_length,
    candidates_mom_muon,
    candidates_mom_prot,
    candidates_type,
    candidates_CNNscore,
    candidates_chi2_prot,
    candidates_chi2_kaon,
    candidates_chi2_muon,
    candidates_chi2_ndf,
    candidates_kaon_PID,
    candidates_kaon_PID_ndf,
    candidates_distance_mother,
    candidates_distance_dau,
    candidates_cos_dau,
    candidates_averagedEdx,
    candidates_truncated_dedx,
    candidates_vtx_michelscore,
    candidates_vtx_nhits,
    candidates_distance_tcp,
    candidates_lkl_proton,
    candidates_freelkl_proton,
    candidates_freelkl_proton_range,

    candidates_nhits,
    candidates_hit_x,
    candidates_hit_y,
    candidates_hit_z,
    candidates_hit_dedx,
    candidates_hit_dedx_cal,
    //candidates_hit_dqdx_raw,    
    candidates_hit_resrange,

    //kaon candidates true info
    candidates_truendau,
    candidates_truegeneration,
    candidates_truepdg,
    candidates_trueorigin,
    candidates_truepos,
    candidates_trueendpos,
    candidates_trueproc,
    candidates_trueendproc,
    candidates_truemom,
    candidates_trueendmom,

    //kaon candidates dau info
    candidates_dau_ndau,
    candidates_dau_pos,
    candidates_dau_dir,
    candidates_dau_endpos,
    candidates_dau_enddir,
    candidates_dau_length,
    candidates_dau_mom_muon,
    candidates_dau_mom_prot,
    candidates_dau_type,
    candidates_dau_CNNscore,
    candidates_dau_chi2_prot,
    candidates_dau_chi2_muon,
    candidates_dau_chi2_ndf,
    candidates_dau_kaon_PID,
    candidates_dau_kaon_PID_ndf,
    candidates_dau_calE,
    candidates_dau_averagedEdx,
    candidates_dau_vtx_michelscore,
    candidates_dau_vtx_nhits,
    candidates_dau_forced,
    candidates_dau_forced_matched,

    candidates_dau_nhits,
    candidates_dau_hit_x,
    candidates_dau_hit_y,
    candidates_dau_hit_z,
    candidates_dau_hit_dedx,
    candidates_dau_hit_dedx_cal,
    //candidates_dau_hit_dqdx_raw,    
    candidates_dau_hit_resrange,

    //kaon candidates dau true info
    candidates_dau_truendau,
    candidates_dau_truepdg,
    candidates_dau_truepos,
    candidates_dau_trueendpos,
    candidates_dau_trueproc,
    candidates_dau_trueendproc,
    candidates_dau_truemom,
    candidates_dau_trueendmom,

    //kaon candidates dau info
    candidates_gdau_ndau,
    candidates_gdau_pos,
    candidates_gdau_dir,
    candidates_gdau_endpos,
    candidates_gdau_enddir,
    candidates_gdau_length,
    candidates_gdau_mom_muon,
    candidates_gdau_mom_prot,
    candidates_gdau_type,
    candidates_gdau_CNNscore,
    candidates_gdau_chi2_prot,
    candidates_gdau_chi2_muon,
    candidates_gdau_chi2_ndf,
    candidates_gdau_calE,
    candidates_gdau_averagedEdx,
    candidates_gdau_vtx_michelscore,
    candidates_gdau_vtx_nhits,

    //kaon candidates gdau true info
    candidates_gdau_truendau,
    candidates_gdau_truepdg,
    candidates_gdau_truepos,
    candidates_gdau_trueendpos,
    candidates_gdau_trueproc,
    candidates_gdau_trueendproc,
    candidates_gdau_truemom,
    candidates_gdau_trueendmom,

    //temporary solution for systematics
    //variables for storing the most advanced candidate in the selection
    //kaon candidates info
    bestcandidate_generation,
    bestcandidate_ndau,
    bestcandidate_pos,
    bestcandidate_dir,
    bestcandidate_endpos,
    bestcandidate_enddir,
    bestcandidate_length,
    bestcandidate_mom_muon,
    bestcandidate_mom_prot,
    bestcandidate_type,
    bestcandidate_CNNscore,
    bestcandidate_chi2_prot,
    bestcandidate_chi2_muon,
    bestcandidate_chi2_kaon,
    bestcandidate_chi2_ndf,
    bestcandidate_distance_mother,
    bestcandidate_distance_dau,
    bestcandidate_cos_dau,
    bestcandidate_averagedEdx,
    bestcandidate_vtx_michelscore,
    bestcandidate_vtx_nhits,
    bestcandidate_calE,

    bestcandidate_chi2_prot_perndf_5,
    bestcandidate_chi2_muon_perndf_5,
    bestcandidate_chi2_kaon_perndf_5,
    bestcandidate_chi2_prot_perndf_10,
    bestcandidate_chi2_muon_perndf_10,
    bestcandidate_chi2_kaon_perndf_10,
    bestcandidate_chi2_prot_perndf_15,
    bestcandidate_chi2_muon_perndf_15,
    bestcandidate_chi2_kaon_perndf_15,
    bestcandidate_chi2_prot_perndf_20,
    bestcandidate_chi2_muon_perndf_20,
    bestcandidate_chi2_kaon_perndf_20,
    bestcandidate_chi2_prot_perndf_25,
    bestcandidate_chi2_muon_perndf_25,
    bestcandidate_chi2_kaon_perndf_25,

    bestcandidate_lkl_proton,
    bestcandidate_freelkl_proton,
    bestcandidate_freelkl_proton_range,
    bestcandidate_lkl_kaon,
    bestcandidate_freelkl_kaon,
    bestcandidate_freelkl_kaon_range,

    bestcandidate_nhits,
    bestcandidate_hit_x,
    bestcandidate_hit_y,
    bestcandidate_hit_z,
    bestcandidate_hit_dedx,
    bestcandidate_hit_dedx_cal,
    //bestcandidate_hit_dqdx_raw,    
    bestcandidate_hit_resrange,

    //kaon bestcandidate true info
    bestcandidate_truendau,
    bestcandidate_truegeneration,
    bestcandidate_truepdg,
    bestcandidate_truepos,
    bestcandidate_trueendpos,
    bestcandidate_trueproc,
    bestcandidate_trueendproc,
    bestcandidate_truemom,
    bestcandidate_trueendmom,
    bestcandidate_trueId,
    bestcandidate_trueendmom_atAPA,

    //kaon bestcandidate dau info
    bestcandidate_dau_ndau,
    bestcandidate_dau_pos,
    bestcandidate_dau_dir,
    bestcandidate_dau_endpos,
    bestcandidate_dau_enddir,
    bestcandidate_dau_length,
    bestcandidate_dau_mom_muon,
    bestcandidate_dau_mom_prot,
    bestcandidate_dau_type,
    bestcandidate_dau_CNNscore,
    bestcandidate_dau_chi2_prot,
    bestcandidate_dau_chi2_muon,
    bestcandidate_dau_chi2_ndf,
    bestcandidate_dau_calE,
    bestcandidate_dau_vtx_michelscore,
    bestcandidate_dau_vtx_nhits,

    bestcandidate_dau_nhits,
    bestcandidate_dau_hit_x,
    bestcandidate_dau_hit_y,
    bestcandidate_dau_hit_z,
    bestcandidate_dau_hit_dedx,
    bestcandidate_dau_hit_dedx_cal,
    //bestcandidate_dau_hit_dqdx_raw,    
    bestcandidate_dau_hit_resrange,

    //kaon bestcandidate dau true info
    bestcandidate_dau_truendau,
    bestcandidate_dau_truepdg,
    bestcandidate_dau_truepos,
    bestcandidate_dau_trueendpos,
    bestcandidate_dau_trueproc,
    bestcandidate_dau_trueendproc,
    bestcandidate_dau_truemom,
    bestcandidate_dau_trueendmom,


    //kaon bestcandidate parent info
    bestcandidate_parent_ndau,
    bestcandidate_parent_pos,
    bestcandidate_parent_dir,
    bestcandidate_parent_endpos,
    bestcandidate_parent_enddir,

    //kaon bestcandidate parent true info
    bestcandidate_parent_trueId,
    bestcandidate_parent_truepos,
    bestcandidate_parent_truedir,
    bestcandidate_parent_trueendpos,
    bestcandidate_parent_trueenddir,
    bestcandidate_parent_trueendmom_atAPA,
    
    enumKaonMicroTreesLast
  };
}



#endif
