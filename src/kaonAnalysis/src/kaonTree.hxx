#ifndef kaonTree_h
#define kaonTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdDataClasses.hxx"
#include "kaonClasses.hxx"
#include "standardPDTree.hxx"

namespace kaonTree{

  // Methods to add to the output tree the kaonAnalysis sets of variables
  //void AddKaonVariables_CountersTrue(OutputManager& output);
  
  void AddKaonVariables_CandidateDaughtersTrue(OutputManager& output, UInt_t nmax);

  void AddKaonVariables_CandidateGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);
  void AddKaonVariables_CandidateGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);

  void AddKaonVariables_CandidateGGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters);
  void AddKaonVariables_CandidateGGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters);

  void AddKaonVariables_TrueKaonCandidates(OutputManager& output, UInt_t nmax);

  void AddKaonVariables_KaonCandidatesReco(OutputManager& output, UInt_t nmax);
  void AddKaonVariables_KaonCandidatesHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void AddKaonVariables_KaonCandidatesTrue(OutputManager& output, UInt_t nmax);
  

  // Methods to fill the kaonAnalysis sets of variables in the output tree
  //void FillKaonVariables_CountersTrue(OutputManager& output, PDCounters& counters);
  void FillKaonVariables_CandidateDaughterTrue(OutputManager& output, AnaParticlePD* part, AnaParticlePD* dau);

  void FillKaonVariables_CandidateGDaughterReco(OutputManager& output, AnaParticlePD* part, Int_t index);
  void FillKaonVariables_CandidateGDaughterTrue(OutputManager& output, AnaParticlePD* part, Int_t index);  

  void FillKaonVariables_CandidateGGDaughterReco(OutputManager& output, AnaParticlePD* part, Int_t index1, Int_t index2);
  void FillKaonVariables_CandidateGGDaughterTrue(OutputManager& output, AnaParticlePD* part, Int_t index1, Int_t index2);  

  //void FillKaonVariables_TrueKaonCandidates(OutputManager& output, AnaTrueParticlePD* truePart);
  void FillKaonVariables_TrueKaonCandidates(OutputManager& output, const kaonAnaTrueVertex& kvtx);

  void FillKaonVariables_KaonCandidatesReco(OutputManager& output, AnaParticlePD* part);
  void FillKaonVariables_KaonCandidatesHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);
  void FillKaonVariables_KaonCandidatesTrue(OutputManager& output, AnaParticlePD* part);

  // Enum with unique indexes for output tree variables  
  enum enumKaonMicroTrees{

    seltrk_dau_truesecondary = standardPDTree::enumStandardMicroTreesLast_standardPDTree,

    // selected track gdaughters reco info
    seltrk_gdau_ndau,
    seltrk_gdau_pos,
    seltrk_gdau_dir,
    seltrk_gdau_endpos,
    seltrk_gdau_enddir,
    seltrk_gdau_length,
    seltrk_gdau_mom_muon,
    seltrk_gdau_mom_prot,
    seltrk_gdau_type,
    seltrk_gdau_CNNscore,
    seltrk_gdau_chi2_prot,
    seltrk_gdau_chi2_muon,
    seltrk_gdau_chi2_ndf,
    seltrk_gdau_nhits,
    seltrk_gdau_hit_dedx,
    seltrk_gdau_hit_resrange,
    //seltrk_gdau_nhits2,
    //seltrk_gdau_hit_dqdx_raw,    
    

    // selected track gdaughters true info
    seltrk_gdau_truendau,
    seltrk_gdau_truepdg,
    seltrk_gdau_truepos,
    seltrk_gdau_trueendpos,
    seltrk_gdau_trueproc,
    seltrk_gdau_trueendproc,
    seltrk_gdau_truemom,
    seltrk_gdau_trueendmom,

    // selected track ggdaughters reco info
    seltrk_ggdau_ndau,
    seltrk_ggdau_posX,
    seltrk_ggdau_posY,
    seltrk_ggdau_posZ,
    seltrk_ggdau_dirX,
    seltrk_ggdau_dirY,
    seltrk_ggdau_dirZ,
    seltrk_ggdau_endposX,
    seltrk_ggdau_endposY,
    seltrk_ggdau_endposZ,
    seltrk_ggdau_enddirX,
    seltrk_ggdau_enddirY,
    seltrk_ggdau_enddirZ,
    seltrk_ggdau_length,
    seltrk_ggdau_mom_muon,
    seltrk_ggdau_mom_prot,
    seltrk_ggdau_type,
    seltrk_ggdau_CNNscore0,
    seltrk_ggdau_CNNscore1,
    seltrk_ggdau_CNNscore2,
    seltrk_ggdau_chi2_prot,
    seltrk_ggdau_chi2_muon,
    seltrk_ggdau_chi2_ndf,
    seltrk_ggdau_nhits,

    // selected track ggdaughters info
    seltrk_ggdau_truepdg,    
    seltrk_ggdau_truendau,   
    seltrk_ggdau_trueposX,
    seltrk_ggdau_trueposY,    
    seltrk_ggdau_trueposZ,    
    seltrk_ggdau_trueendposX, 
    seltrk_ggdau_trueendposY, 
    seltrk_ggdau_trueendposZ, 
    seltrk_ggdau_trueproc,   
    seltrk_ggdau_trueendproc,
    seltrk_ggdau_truemom,    
    seltrk_ggdau_trueendmom, 

    //true kaon candidates info
    truenkaons,
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
    candidates,
    candidates_ndau,
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
    candidates_chi2_muon,
    candidates_chi2_ndf,

    candidates_nhits,
    candidates_hit_x,
    candidates_hit_y,
    candidates_hit_z,
    candidates_hit_dedx,
    //candidates_hit_dedx_cal,
    //candidates_hit_dqdx_raw,    
    candidates_hit_resrange,

    //kaon candidates true info
    candidates_truendau,
    candidates_truepdg,
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

    candidates_dau_nhits,
    candidates_dau_hit_x,
    candidates_dau_hit_y,
    candidates_dau_hit_z,
    candidates_dau_hit_dedx,
    //candidates_dau_hit_dedx_cal,
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

    enumKaonMicroTreesLast
  };
}



#endif
