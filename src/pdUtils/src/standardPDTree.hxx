#ifndef standardPDTree_h
#define standardPDTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdDataClasses.hxx"

namespace standardPDTree{

  // Methods to add to the output tree the standard sets of variables
  void AddStandardVariables_CountersTrue(OutputManager& output);
  
  void AddStandardVariables_BeamReco(OutputManager& output);
  void AddStandardVariables_BeamTrue(OutputManager& output);

  void AddStandardVariables_CandidateReco(OutputManager& output);
  void AddStandardVariables_CandidateTrue(OutputManager& output);
  void AddStandardVariables_CandidateHitsReco(OutputManager& output);

  void AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_AllParticlesTrue(OutputManager& output, UInt_t nmax);

  void AddStandardVariables_CandidateDaughtersReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_CandidateDaughtersTrue(OutputManager& output, UInt_t nmax);

  void AddStandardVariables_CandidateGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);
  void AddStandardVariables_CandidateGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters);


  // Methods to fill the standard sets of variables in the output tree
  void FillStandardVariables_CountersTrue(OutputManager& output, PDCounters& counters);
  
  void FillStandardVariables_BeamTrue(OutputManager& output, AnaBeamB* beamB);
  void FillStandardVariables_BeamReco(OutputManager& output, AnaBeamB* beamB);

  void FillStandardVariables_AllParticlesReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_CandidateTrue(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateHitsReco(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_CandidateDaughterReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateDaughterTrue(OutputManager& output, AnaParticlePD* part);  

  void FillStandardVariables_CandidateGDaughterReco(OutputManager& output, AnaParticlePD* part, Int_t index);
  void FillStandardVariables_CandidateGDaughterTrue(OutputManager& output, AnaParticlePD* part, Int_t index);  


  // Enum with unique indexes for output tree variables  
  enum enumStandardMicroTrees_standardPDTree{

    // selected track true info
    seltrk_truemom = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    seltrk_trueendmom,
    seltrk_truepdg,
    seltrk_trueproc,
    seltrk_trueendproc,
    seltrk_truendau,  

    // selected track reco info
    seltrk_mom_muon,
    seltrk_mom_prot,
    seltrk_mom_muon_alt,
    seltrk_mom_prot_alt,
    seltrk_costheta,
    seltrk_nhits,
    seltrk_length,
    seltrk_length_alt,
    seltrk_length_raw,
    seltrk_ndau,
    seltrk_dqdx,
    seltrk_dedx,
    seltrk_dedx_raw,
    //    seltrk_dedx_exp,
    seltrk_truepos,
    seltrk_trueendpos,
    seltrk_truedir,
    seltrk_trueeff,
    seltrk_truepur,
    seltrk_dir,
    seltrk_enddir, 
    seltrk_pos,
    seltrk_endpos,
    seltrk_pid,
    seltrk_calo,

    seltrk_hit_dedx,
    seltrk_hit_dedx_raw,
    seltrk_hit_dedx_cal,
    seltrk_hit_dqdx,
    seltrk_hit_dqdx_raw,
    seltrk_hit_dqdx_noSCE,
    seltrk_hit_x,
    seltrk_hit_y,
    seltrk_hit_z,
    seltrk_hit_x_raw,
    seltrk_hit_y_raw,
    seltrk_hit_z_raw,
    seltrk_hit_resrange,
    seltrk_hit_cnn,
    seltrk_nhitsperplane,

    // selected track daughters reco info
    seltrk_dau_ndau,
    seltrk_dau_mom_muon,
    seltrk_dau_mom_prot,
    seltrk_dau_pos,
    seltrk_dau_dir,
    seltrk_dau_endpos,
    seltrk_dau_enddir,
    seltrk_dau_length,
    seltrk_dau_nhits,
    seltrk_dau_nhits2,
    seltrk_dau_hit_dedx,
    seltrk_dau_hit_dqdx_raw,    
    seltrk_dau_hit_resrange,
    seltrk_dau_hit_cnn,

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
    seltrk_gdau_nhits,
    seltrk_gdau_nhits2,
    seltrk_gdau_hit_dedx,
    seltrk_gdau_hit_dqdx_raw,    
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

    // all particles in the event, true info
    ntracks,
    trk_truemom,
    trk_trueendmom,
    trk_truepdg,
    trk_truepos,
    trk_trueendpos,
    trk_trueproc,
    trk_trueendproc,
    trk_truedir,

    // all particles in the event, reco info
    trk_length,
    trk_dedx,
    trk_mom_muon,
    trk_mom_prot,
    trk_pos,
    trk_endpos,
    trk_dir,
    trk_ndau,

    // beam true info
    beam_truepos,
    beam_trueendpos,
    beam_truemom,
    beam_truemom_tpc,
    beam_truedir,
    beam_truepdg,
    beam_trueendproc,

    // beam reco info
    beam_endpos,          
    beam_enddir,          
    beam_mom,                                 
    beam_mom_raw,                                 
    beam_mom_tpc,                                 
    beam_mom_tpc_raw,                                 
    beam_trigger,      
    beam_tof,
    beam_time,                                 
    beam_ckov_status,  
    beam_ckov_time,    
    beam_ckov_pressure,
    beam_npdgs,
    beam_pdgs,
    
    // global counters
    truebeamdau_npi0,
    truebeamdau_npiplus,
    truebeamdau_npiminus,
    truebeamdau_nproton,
    truebeamdau_nneutron,
    truebeamdau_nnucleus,
        
    enumStandardMicroTreesLast_standardPDTree
  };
}



#endif
