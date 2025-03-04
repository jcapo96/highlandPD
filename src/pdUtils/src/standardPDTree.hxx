#ifndef standardPDTree_h
#define standardPDTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdDataClasses.hxx"

namespace standardPDTree{

  // Methods to add to the output tree the standard sets of variables
  void AddStandardVariables_EventInfo(OutputManager& output);

  void AddStandardVariables_BeamInstrumentationReco(OutputManager& output);
  void AddStandardVariables_BeamInstrumentationTrue(OutputManager& output);

  void AddStandardVariables_BeamParticleReco(OutputManager& output);
  void AddStandardVariables_BeamParticleTrue(OutputManager& output);
  void AddStandardVariables_BeamParticleHitsReco(OutputManager& output);

  void AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_AllParticlesTrue(OutputManager& output, UInt_t nmax);

  void AddStandardVariables_BeamParticleDaughtersReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_BeamParticleDaughtersTrue(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_BeamParticleDaughtersHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);

  void AddStandardVariables_BeamParticleGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);
  void AddStandardVariables_BeamParticleGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);

  void AddStandardVariables_BeamParticleGGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters);
  void AddStandardVariables_BeamParticleGGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters);

  // Methods to fill the standard sets of variables in the output tree
  void FillStandardVariables_CountersTrue(OutputManager& output, PDCounters& counters);
  
  void FillStandardVariables_EventInfo(OutputManager& output, AnaEventInfoPD* info);

  void FillStandardVariables_BeamInstrumentationTrue(OutputManager& output, AnaBeamB* beamB);
  void FillStandardVariables_BeamInstrumentationReco(OutputManager& output, AnaBeamB* beamB);

  void FillStandardVariables_AllParticlesReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_BeamParticleTrue(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beamPart = NULL);
  void FillStandardVariables_BeamParticleHitsReco(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_BeamParticleDaughtersReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_BeamParticleDaughtersTrue(OutputManager& output, AnaParticlePD* part);  
  void FillStandardVariables_BeamParticleDaughtersHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxsavedhits);

  void FillStandardVariables_BeamParticleGDaughtersReco(OutputManager& output, AnaParticlePD* part, Int_t index);
  void FillStandardVariables_BeamParticleGDaughtersTrue(OutputManager& output, AnaParticlePD* part, Int_t index);  

  void FillStandardVariables_BeamParticleGGDaughtersReco(OutputManager& output, AnaParticlePD* part, Int_t index1, Int_t index2);
  void FillStandardVariables_BeamParticleGGDaughtersTrue(OutputManager& output, AnaParticlePD* part, Int_t index1, Int_t index2);  

  // Enum with unique indexes for output tree variables  
  enum enumStandardMicroTrees_standardPDTree{

    // selected track (beam particle) true info
    seltrk_truemom = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    seltrk_trueendmom,
    seltrk_truepdg,
    seltrk_trueproc,
    seltrk_trueendproc,
    seltrk_truendau,
    seltrk_truepos,
    seltrk_trueendpos,
    seltrk_truedir,
    seltrk_trueeff,
    seltrk_truepur,
    seltrk_trueId,
    seltrk_true_matched,

    // selected track reco info
    seltrk_ndau,
    seltrk_pos,
    seltrk_endpos,
    seltrk_dir,
    seltrk_enddir,
    seltrk_length,
    seltrk_mom_muon,
    seltrk_mom_prot,
    seltrk_csdarange_muon,
    seltrk_csdarange_prot,
    seltrk_nhits,
    seltrk_CNNscore,
    seltrk_chi2_prot,
    seltrk_chi2_kaon,
    seltrk_chi2_muon,
    seltrk_chi2_ndf,
    seltrk_truncated_dedx,
    seltrk_calE,

    seltrk_hit_dedx,
    seltrk_hit_dqdx,
    seltrk_hit_x,
    seltrk_hit_y,
    seltrk_hit_z,
    seltrk_hit_resrange,

    // selected track daughters reco info
    seltrk_dau_ndau,
    seltrk_dau_pos,
    seltrk_dau_dir,
    seltrk_dau_endpos,
    seltrk_dau_enddir,
    seltrk_dau_length,
    seltrk_dau_mom_muon,
    seltrk_dau_mom_prot,
    seltrk_dau_CNNscore,
    seltrk_dau_chi2_prot,
    seltrk_dau_chi2_kaon,
    seltrk_dau_chi2_muon,
    seltrk_dau_chi2_ndf,
    seltrk_dau_calE,
    seltrk_dau_truncated_dedx,
    seltrk_dau_nhits,
    seltrk_dau_hit_x,
    seltrk_dau_hit_y,
    seltrk_dau_hit_z,
    seltrk_dau_hit_dedx,
    seltrk_dau_hit_resrange,

    // selected track daughters true info
    seltrk_dau_truendau,
    seltrk_dau_truepdg,
    seltrk_dau_truepos,
    seltrk_dau_trueendpos,
    seltrk_dau_trueproc,
    seltrk_dau_trueendproc,
    seltrk_dau_truemom,
    seltrk_dau_trueendmom,

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
    seltrk_gdau_chi2_kaon,
    seltrk_gdau_chi2_muon,
    seltrk_gdau_chi2_ndf,
    seltrk_gdau_nhits,
    seltrk_gdau_hit_dedx,
    seltrk_gdau_hit_resrange,
    
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

    // all particles in the event, reco info
    ntracks,
    trk_generation,
    trk_ndau,
    trk_pos,
    trk_dir,
    trk_endpos,
    trk_enddir,
    trk_length,
    trk_mom_muon,
    trk_mom_prot,
    trk_type,
    trk_CNNscore,
    trk_chi2_prot,
    trk_chi2_kaon,
    trk_chi2_muon,
    trk_chi2_ndf,
    trk_nhits,

    // all particles in the event, true info    
    trk_truendau,
    trk_truegeneration,
    trk_truepdg,
    trk_trueorigin,
    trk_truepos,
    trk_trueendpos,
    trk_trueproc,
    trk_trueendproc,
    trk_truemom,
    trk_trueendmom,

    // beam instrumentation true info
    beam_truepos,
    beam_trueendpos,
    beam_truemom,
    beam_truedir,
    beam_truepdg,
    beam_trueendproc,

    // beam reco info
    beam_endpos,          
    beam_enddir,          
    beam_mom,                                 
    beam_nominal_mom,                                 
    beam_tof,
    beam_pdg,
    beam_ntracks,
    
    enumStandardMicroTreesLast_standardPDTree
  };
}



#endif
