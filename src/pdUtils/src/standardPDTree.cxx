#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void standardPDTree::AddStandardVariables_CountersTrue(OutputManager& output){
//********************************************************************
  
  AddVarI(  output, truebeamdau_npi0,     "number of true beam daughter pi0");
  AddVarI(  output, truebeamdau_npiplus,  "number of true beam daughter pi+");
  AddVarI(  output, truebeamdau_npiminus, "number of true beam daughter pi-");
  AddVarI(  output, truebeamdau_nproton,  "number of true beam daughter proton");
  AddVarI(  output, truebeamdau_nneutron, "number of true beam daughter neutron");
  AddVarI(  output, truebeamdau_nnucleus, "number of true beam daughter nucleus");
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamTrue(OutputManager& output){
//********************************************************************

  AddVarI(  output, beam_truepdg,      "beam track true pdg");

  AddVar4VF(output, beam_truepos,      "beam track true start position");
  AddVar3VF(output, beam_truedir,      "beam track true direction");
  AddVarF(  output, beam_truemom,      "beam track true momentum");
  AddVarF(  output, beam_truemom_tpc,  "beam track true momentum in TPC");

  AddVar4VF(output, beam_trueendpos,   "beam track true end position");
  AddVarI(  output, beam_trueendproc,  "beam track true end process");
}
  
//********************************************************************
void standardPDTree::AddStandardVariables_BeamReco(OutputManager& output){
//********************************************************************
  
  AddVar3VF(output, beam_endpos,           "beam track end position");
  AddVar3VF(output, beam_enddir,           "beam track end direction");
  AddVarF(  output, beam_mom,              "beam track momentum");
  AddVarF(  output, beam_mom_tpc,          "beam track momentum after crossing cryowall");
  AddVarF(  output, beam_mom_raw,          "beam track momentum without BI corrections");
  //AddVarF(  output, beam_mom_tpc_raw,      "beam track momentum after crossing cryowall without BI corrections");

  AddVarI(  output, beam_trigger,          "beam trigger");
  AddVarD(  output, beam_tof,              "beam TOF");
  AddVarD(  output, beam_time,             "beam track time");

  AddVarFixVI(output,  beam_ckov_status,   "beam ckov status",  2);
  AddVarFixVD(output,  beam_ckov_time,     "beam ckov time",    2);
  AddVarFixVD(output,  beam_ckov_pressure, "beam ckov pressure",2);  

  AddVarVI(  output, beam_pdgs,             "beam pdgs", beam_npdgs);
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateReco(OutputManager& output){
//********************************************************************

  AddVar4VF(  output, seltrk_pos,            "candidate reconstructed position");
  AddVar3VF(  output, seltrk_dir,            "candidate reconstructed direction");
  AddVar4VF(  output, seltrk_endpos,         "candidate reconstructed end position");
  AddVar3VF(  output, seltrk_enddir,         "candidate reconstructed end direction");
  AddVarF(    output, seltrk_costheta,       "candidate reconstructed cos(theta)");
  AddVarI(    output, seltrk_nhits,          "candidate number of hits");
  AddVarF(    output, seltrk_length,         "candidate length");
  AddVarF(    output, seltrk_length_alt,     "candidate alternate length");
  AddVarF(    output, seltrk_length_raw,     "candidate length");
  AddVarFixVI(output, seltrk_nhitsperplane,  "candidate number of hits per plane",3);

  AddVarF(    output, seltrk_mom_muon,       "candidate momentum muon");
  AddVarF(    output, seltrk_mom_prot,       "candidate momentum proton");
  AddVarF(    output, seltrk_mom_muon_alt,   "candidate alternate momentum muon");
  AddVarF(    output, seltrk_mom_prot_alt,   "candidate alternate momentum proton");
  AddVarF(    output, seltrk_dqdx,           "candidate average dQdx");
  AddVarF(    output, seltrk_dedx,           "candidate average dEdx");
  AddVarF(    output, seltrk_dedx_raw,       "candidate average dEdx");

  AddVarFixMF(output, seltrk_pid,            "candidate PID variables", 3,8);
  AddVarFixMF(output, seltrk_calo,           "candidate CALO variables",3,5); 
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateHitsReco(OutputManager& output){
//********************************************************************

  AddVarFixMF(output, seltrk_hit_x,          "candidate x per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_y,          "candidate y per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_z,          "candidate z per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_x_raw,      "candidate x per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_y_raw,      "candidate y per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_z_raw,      "candidate z per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dedx_raw,   "candidate dEdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dedx,       "candidate calibrated dEdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dedx_cal,   "candidate LAr calibrated dEdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dqdx_raw,   "candidate dQdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dqdx,       "candidate calibrated dQdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dqdx_noSCE, "candidate LAr noSCE dQdx per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_resrange,   "candidate Residual Range per hit",3,NMAXHITSPERPLANE_SELTRK);
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateTrue(OutputManager& output){
//********************************************************************

  AddVar4VF(output, seltrk_truepos,     "candidate true position");         
  AddVar3VF(output, seltrk_truedir,     "candidate true direction");         
  AddVarF(  output, seltrk_truemom,     "candidate true momentum");
  AddVarI(  output, seltrk_trueproc,    "candidate true initial process");
  AddVarI(  output, seltrk_truendau,    "candidate true number of daughters");

  AddVar4VF(output, seltrk_trueendpos,  "candidate true end position");      
  AddVarF(  output, seltrk_trueendmom,  "candidate true end momentum");
  AddVarI(  output, seltrk_trueendproc, "candidate true final process");
  
  AddVarI(  output, seltrk_truepdg,     "candidate true pdg code");
  AddVarF(  output, seltrk_trueeff,     "candidate true-reco matching eff");
  AddVarF(  output, seltrk_truepur,     "candidate true-reco matching pur");
}

//********************************************************************
void standardPDTree::AddStandardVariables_AllParticlesTrue(OutputManager& output, UInt_t nmax){
//********************************************************************
  
  AddVarMaxSizeVI(  output, trk_truepdg,      "pdg code",               ntracks,nmax);

  AddVarMaxSize4MF( output, trk_truepos,      "true start position",    ntracks,nmax);
  AddVarMaxSize3MF( output, trk_truedir,      "true start direction",   ntracks,nmax);
  AddVarMaxSizeVF(  output, trk_truemom,      "true momentum",          ntracks,nmax);
  AddVarMaxSizeVI(  output, trk_trueproc,     "initial process",        ntracks,nmax);

  AddVarMaxSize4MF( output, trk_trueendpos,   "true end position",      ntracks,nmax);
  AddVarMaxSizeVF(  output, trk_trueendmom,   "true end momentum",      ntracks,nmax);
  AddVarMaxSizeVI(  output, trk_trueendproc,  "final process",          ntracks,nmax);

}

//********************************************************************
void standardPDTree::AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSize4MF( output, trk_pos,          "start position",         ntracks,nmax);
  AddVarMaxSize4MF( output, trk_endpos,       "end position",           ntracks,nmax);
  AddVarMaxSize3MF( output, trk_dir,          "start direction",        ntracks,nmax);  
  AddVarMaxSizeVF(  output, trk_length,       "length",                 ntracks,nmax);
  AddVarMaxSizeVF(  output, trk_dedx,         "average dedx",           ntracks,nmax);
  AddVarMaxSizeVF(  output, trk_mom_muon,     "range momentum (muon)",  ntracks,nmax);
  AddVarMaxSizeVF(  output, trk_mom_prot,     "range momentum (proton)",ntracks,nmax);
  AddVarMaxSizeVI(  output, trk_ndau,         "#daughters",             ntracks,nmax);

}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateDaughtersTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSize4MF(output, seltrk_dau_truepos,     "daugthers true position",    seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_trueendpos,  "daugthers true end position",seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_trueproc,    "daugthers true process",     seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_trueendproc, "daughters true end process", seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_truemom,     "daughters true momentum",    seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_trueendmom,  "daughters true end momentum",seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_truepdg,     "daughters true pdg",         seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_truendau,    "daughters' true ndaughters", seltrk_ndau,nmax);
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateDaughtersReco(OutputManager& output, UInt_t nmax){
//********************************************************************
  
  AddVarMaxSizeVF(output,  seltrk_dau_mom_muon,   "daughters momentum (muon)",   seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_mom_prot,   "daughters momentum (proton)", seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_pos,        "daughters position",          seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_dir,        "daughters direction",         seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_endpos,     "daughters position",          seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_enddir,     "daughters direction",         seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_length,     "daughters length",            seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_nhits,      "daughters #hits",             seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_nhits2,     "daughters hits in plane 2",   seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_ndau,       "daughters' daughters"     ,   seltrk_ndau,nmax);

  AddVarMF(output, seltrk_dau_hit_dedx,     "daughters dEdx per hit",       seltrk_ndau,-nmax,NMAXHITSPERPLANE);
  AddVarMF(output, seltrk_dau_hit_dqdx_raw, "daughters raw dQdx per hit",   seltrk_ndau,-nmax,NMAXHITSPERPLANE);
  AddVarMF(output, seltrk_dau_hit_resrange, "daughters residual range",     seltrk_ndau,-nmax,NMAXHITSPERPLANE);
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters){
//********************************************************************

  AddVarMI(   output, seltrk_gdau_truepdg,     "gdaughters true pdg",                 seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI(   output, seltrk_gdau_truendau,    "gdaughters true number of daughters", seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_truepos,     "gdaughters true position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_trueendpos,  "gdaughters true position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI(   output, seltrk_gdau_trueproc,    "gdaughters true process",             seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI(   output, seltrk_gdau_trueendproc, "gdaughters true end process",         seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF(   output, seltrk_gdau_truemom,     "gdaughters true mom",                 seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF(   output, seltrk_gdau_trueendmom,  "gdaughters true end mom",             seltrk_ndau, -nmax, nmaxgdaughters);

}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters){
//********************************************************************

  AddVarMI(   output, seltrk_gdau_ndau,         "gdaughters number of daughters", seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_pos,          "gdaughters position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D3MF(output, seltrk_gdau_dir,          "gdaughters direction",           seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_endpos,       "gdaughters end position",        seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D3MF(output, seltrk_gdau_enddir,       "gdaughters end direction",       seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF(   output, seltrk_gdau_length,       "gdaughters length",              seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI(   output, seltrk_gdau_nhits,        "gdaughters nhits",               seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3DMF( output, seltrk_gdau_hit_dedx,     "gdaughters hit dedx",            seltrk_ndau, -nmax, nmaxgdaughters, NMAXHITSPERPLANE);
  AddVar3DMF( output, seltrk_gdau_hit_resrange, "gdaughters hit residual range",  seltrk_ndau, -nmax, nmaxgdaughters, NMAXHITSPERPLANE);
}  



//********************************************************************
void standardPDTree::AddStandardVariables_CandidateGGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters){
//********************************************************************

  AddVar3DMI(output, seltrk_ggdau_truepdg,     "ggdaughters true pdg",                 seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_truendau,    "ggdaughters true number of daughters", seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueposX,    "ggdaughters true position X",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueposY,    "ggdaughters true position Y",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueposZ,    "ggdaughters true position Z",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueendposX, "ggdaughters true position X",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueendposY, "ggdaughters true position Y",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueendposZ, "ggdaughters true position Z",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_trueproc,    "ggdaughters true process",             seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_trueendproc, "ggdaughters true end process",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_truemom,     "ggdaughters true mom",                 seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_trueendmom,  "ggdaughters true end mom",             seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateGGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters){
//********************************************************************

  AddVar3DMI(output, seltrk_ggdau_ndau,    "ggdaughters number of daughters", seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posX,    "ggdaughters position X",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posY,    "ggdaughters position Y",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posZ,    "ggdaughters position Z",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirX,    "ggdaughters direction X",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirY,    "ggdaughters direction Y",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirZ,    "ggdaughters direction Z",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposX, "ggdaughters end position X",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposY, "ggdaughters end position Y",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposZ, "ggdaughters end position Z",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirX, "ggdaughters end direction X",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirY, "ggdaughters end direction Y",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirZ, "ggdaughters end direction Z",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_length,  "ggdaughters length",              seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_nhits,   "ggdaughters #hits",               seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
}

//********************************************************************
void standardPDTree::FillStandardVariables_CountersTrue(OutputManager& output, PDCounters& counters){
//********************************************************************

  output.FillVar(truebeamdau_npiplus,                  counters.ntrue_beamdaughter_piplus);
  output.FillVar(truebeamdau_npiminus,                 counters.ntrue_beamdaughter_piminus);
  output.FillVar(truebeamdau_npi0,                     counters.ntrue_beamdaughter_pi0);
  output.FillVar(truebeamdau_nproton,                  counters.ntrue_beamdaughter_proton);
  output.FillVar(truebeamdau_nneutron,                 counters.ntrue_beamdaughter_neutron);
  output.FillVar(truebeamdau_nnucleus,                 counters.ntrue_beamdaughter_nucleus);
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamTrue(OutputManager& output, AnaBeamB* beamB){
//********************************************************************

  if (!beamB) return;  
  AnaBeam* beam = static_cast<AnaBeam*>(beamB);
  if (!beam->BeamParticle) return;
  AnaTrueParticlePD* beamTruePart= static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);   
  if (!beamTruePart) return;

  //      output.FillVar(beam_truelengthInTPC,                 beamTruePart->Length);
  output.FillVar(beam_truemom,                    beamTruePart->Momentum);
  output.FillVar(beam_truemom_tpc,                beamTruePart->MomentumInTPC);
  output.FillVar(beam_truepdg,                    beamTruePart->PDG);
  output.FillVar(beam_trueendproc,                beamTruePart->ProcessEnd);
  output.FillVectorVarFromArray(beam_truedir,     beamTruePart->Direction,   3);
  output.FillVectorVarFromArray(beam_truepos,     beamTruePart->Position, 4);
  output.FillVectorVarFromArray(beam_trueendpos,  beamTruePart->PositionEnd, 4);
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamReco(OutputManager& output, AnaBeamB* beamB){
//********************************************************************

  if (!beamB) return;
  AnaBeamPD* beam         = static_cast<AnaBeamPD*>(beamB);

  output.FillVar(beam_tof,     beam->TOF);
  output.FillVar(beam_trigger, beam->BeamTrigger);
  output.FillVar(beam_time,    beam->BeamTrackTime);

  output.FillVectorVarFromArray(beam_ckov_pressure,       beam->CerenkovPressure,2);
  output.FillVectorVarFromArray(beam_ckov_status,         beam->CerenkovStatus,  2);
  output.FillVectorVarFromArray(beam_ckov_time,           beam->CerenkovTime,    2);

  AnaParticleMomB* beamPart = beam->BeamParticle;    
  if (!beamPart) return;
  
  output.FillVectorVarFromArray(beam_endpos,         beamPart->PositionEnd,3);
  output.FillVectorVarFromArray(beam_enddir,         beamPart->DirectionEnd,3);
  output.FillVar(               beam_mom,            beamPart->Momentum);   
  output.FillVar(               beam_mom_tpc,        (Float_t)beam->BeamMomentumInTPC);
  
  if (beamPart->TrueObject)output.FillVar(               beam_mom_raw,        beamPart->Momentum);   
  else                     output.FillVar(               beam_mom_raw,        static_cast<const AnaParticle*>(beamPart->Original->Original)->Momentum);   

  for (UInt_t i=0;i<beam->PDGs.size();i++){
    output.FillVectorVar(beam_pdgs,beam->PDGs[i]);
    output.IncrementCounter(beam_npdgs);
  }
  
}

//********************************************************************
void standardPDTree::FillStandardVariables_AllParticlesReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  
  output.FillMatrixVarFromArray(trk_pos,   part->PositionStart,4);
  output.FillMatrixVarFromArray(trk_endpos,part->PositionEnd,4);
  output.FillMatrixVarFromArray(trk_dir,   part->DirectionStart,3);
  output.FillVectorVar(trk_mom_prot,       part->RangeMomentum[0]);
  output.FillVectorVar(trk_mom_muon,       part->RangeMomentum[1]);
  output.FillVectorVar(trk_dedx,           part->AveragedEdx);
  output.FillVectorVar(trk_length,         part->Length);
  output.FillVectorVar(trk_ndau,    (Int_t)part->Daughters.size());
}
  
//********************************************************************
void standardPDTree::FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  if (!part->TrueObject) return;
  
  if (part->TrueObject){
    output.FillMatrixVarFromArray(trk_truepos,    part->GetTrueParticle()->Position,4);
    output.FillMatrixVarFromArray(trk_trueendpos, part->GetTrueParticle()->PositionEnd,4);
    output.FillMatrixVarFromArray(trk_truedir,    part->GetTrueParticle()->Direction,3);
    output.FillVectorVar(trk_truemom,             part->GetTrueParticle()->Momentum);
    output.FillVectorVar(trk_trueendmom,          part->GetTrueParticle()->MomentumEnd);
    output.FillVectorVar(trk_truepdg,             part->GetTrueParticle()->PDG);
    output.FillVectorVar(trk_trueproc,     (Int_t)part->GetTrueParticle()->ProcessStart);
    output.FillVectorVar(trk_trueendproc,  (Int_t)part->GetTrueParticle()->ProcessEnd);    
  }
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  
  AnaTrueParticle* truePart = static_cast<AnaTrueParticle*>(part->TrueObject);
  if (!truePart) return;

  output.FillVar(seltrk_truepdg,                  truePart->PDG);

  output.FillVectorVarFromArray(seltrk_truepos,   truePart->Position, 4);
  output.FillVectorVarFromArray(seltrk_trueendpos,truePart->PositionEnd, 4);

  output.FillVectorVarFromArray(seltrk_truedir,   truePart->Direction, 3);

  output.FillVar(seltrk_truemom,                  truePart->Momentum);
  output.FillVar(seltrk_trueendmom,               truePart->MomentumEnd);
    
  output.FillVar(seltrk_trueproc,          (Int_t)truePart->ProcessStart);
  output.FillVar(seltrk_trueendproc,       (Int_t)truePart->ProcessEnd);
    
  output.FillVar(seltrk_truepur,                  part->TruePur);
  output.FillVar(seltrk_trueeff,                  part->TrueEff);

  output.FillVar(seltrk_truendau,          (Int_t)truePart->Daughters.size());
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;

  output.FillVar(seltrk_mom_prot,          (Float_t)part->RangeMomentum[0]);
  output.FillVar(seltrk_mom_muon,          (Float_t)part->RangeMomentum[1]);
  output.FillVar(seltrk_mom_prot_alt,      (Float_t)part->RangeMomentum_alt[0]);
  output.FillVar(seltrk_mom_muon_alt,      (Float_t)part->RangeMomentum_alt[1]);

  output.FillVar(seltrk_dedx,              part->AveragedEdx);
  output.FillVar(seltrk_dqdx,              part->AveragedQdx);

  output.FillVar(seltrk_dedx_raw,              static_cast<const AnaParticlePD*>(part->Original)->AveragedEdx);
  output.FillVar(seltrk_nhits,                 part->NHits);
  output.FillVar(seltrk_length_alt,            part->Length_alt);
  output.FillVar(seltrk_length,                pdAnaUtils::ComputeTrackLengthFromHitPosition(part));
  output.FillVar(seltrk_length_raw,            part->Original->Length);
  output.FillVar(seltrk_costheta,              part->DirectionStart[2]);
  output.FillVectorVarFromArray(seltrk_pos,    part->PositionStart, 4);
  output.FillVectorVarFromArray(seltrk_endpos, part->PositionEnd, 4);
  output.FillVectorVarFromArray(seltrk_dir,    part->DirectionStart, 3);
  output.FillVectorVarFromArray(seltrk_enddir, part->DirectionEnd, 3);

  for (int i=0;i<3;i++){
    output.FillMatrixVarFromArray(seltrk_pid,    part->PID[i],  i, 8);
    output.FillMatrixVarFromArray(seltrk_calo,   part->CALO[i], i, 5);      
  }
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateHitsReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  for (int i = 0; i < 3; i++){
    if(part->Hits[i].empty()){
      for (int j = 0; j < (int)NMAXHITSPERPLANE_SELTRK; j++){
        output.FillMatrixVar(seltrk_hit_x       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_y       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_z       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_x_raw   ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_y_raw   ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_z_raw   ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dedx    ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dedx_cal,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dqdx    ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dqdx_noSCE,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_resrange,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dedx_raw,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dqdx_raw,(Float_t)-999., i, j);
      }
    }
    else{
      int jmin = std::max<int>(0,(int)part->Hits[i].size()-NMAXHITSPERPLANE_SELTRK);
      int jmax = (int)part->Hits[i].size();
      for (int j = jmin; j < jmax; j++){
        output.FillMatrixVar(seltrk_hit_x    ,(Float_t)part->Hits[i][j].Position.X(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_y    ,(Float_t)part->Hits[i][j].Position.Y(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_z    ,(Float_t)part->Hits[i][j].Position.Z(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_x_raw,(Float_t)static_cast<const AnaParticlePD*>(part->Original->Original->Original)->Hits[i][j].Position.X(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_y_raw,(Float_t)static_cast<const AnaParticlePD*>(part->Original->Original->Original)->Hits[i][j].Position.Y(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_z_raw,(Float_t)static_cast<const AnaParticlePD*>(part->Original->Original->Original)->Hits[i][j].Position.Z(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dedx,       part->Hits[i][j].dEdx,              i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dedx_cal,   part->Hits[i][j].dEdx_calib,        i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dqdx,       part->Hits[i][j].dQdx,              i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dqdx_noSCE, part->Hits[i][j].dQdx_NoSCE,        i, j-jmin);
        output.FillMatrixVar(seltrk_hit_resrange,   part->Hits[i][j].ResidualRange,     i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dedx_raw,  static_cast<const AnaParticlePD*>(part->Original->Original->Original)->Hits[i][j].dEdx, i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dqdx_raw,  static_cast<const AnaParticlePD*>(part->Original->Original->Original)->Hits[i][j].dQdx, i, j-jmin);
      }
    }
    output.FillVectorVarFromArray(seltrk_nhitsperplane, part->NHitsPerPlane,3);
  }
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateDaughterReco(OutputManager& output, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;
  output.FillVectorVar(seltrk_dau_mom_prot,        dau->RangeMomentum[0]);
  output.FillVectorVar(seltrk_dau_mom_muon,        dau->RangeMomentum[1]);
  output.FillMatrixVarFromArray(seltrk_dau_pos,    dau->PositionStart,4);
  output.FillMatrixVarFromArray(seltrk_dau_dir,    dau->DirectionStart,3); 
  output.FillMatrixVarFromArray(seltrk_dau_endpos, dau->PositionEnd,4);
  output.FillMatrixVarFromArray(seltrk_dau_enddir, dau->DirectionEnd,3); 
  output.FillVectorVar(seltrk_dau_length,          dau->Length);
  output.FillVectorVar(seltrk_dau_nhits2,          dau->NHitsPerPlane[2] );
  output.FillVectorVar(seltrk_dau_nhits,           dau->NHits);
  output.FillVectorVar(seltrk_dau_ndau,     (Int_t)dau->Daughters.size());

  if(dau->Hits[2].empty()){
    for (int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(seltrk_dau_hit_dedx,     (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_dqdx_raw, (Float_t)-999., -1, j);  
      output.FillMatrixVar(seltrk_dau_hit_resrange, (Float_t)-999., -1, j);
    }
  }
  else{
    int jmin = std::max<int>(0,(int)dau->Hits[2].size()-NMAXHITSPERPLANE);
    int jmax = (int)dau->Hits[2].size();
    for (int j = jmin; j < jmax; j++){
      output.FillMatrixVar(seltrk_dau_hit_dedx,      dau->Hits[2][j].dEdx, -1, j-jmin);
      //output.FillMatrixVar(seltrk_dau_hit_dqdx_raw,  static_cast<const AnaParticlePD*>(dau->Original->Original->Original)->Hits[2][j].dQdx, -1, j);  
      output.FillMatrixVar(seltrk_dau_hit_resrange,  dau->Hits[2][j].ResidualRange, -1, j-jmin);
    }
  }
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateDaughterTrue(OutputManager& output, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;  
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  if(!dauTruePart) return;
  
  output.FillVectorVar(seltrk_dau_truepdg,             dauTruePart->PDG);
  output.FillVectorVar(seltrk_dau_trueproc,     (Int_t)dauTruePart->ProcessStart);
  output.FillVectorVar(seltrk_dau_trueendproc,  (Int_t)dauTruePart->ProcessEnd);
  output.FillMatrixVarFromArray(seltrk_dau_truepos,    dauTruePart->Position,4);
  output.FillMatrixVarFromArray(seltrk_dau_trueendpos, dauTruePart->PositionEnd,4);
  output.FillVectorVar(seltrk_dau_truemom,             dauTruePart->Momentum);
  output.FillVectorVar(seltrk_dau_trueendmom,          dauTruePart->MomentumEnd);
  output.FillVectorVar(seltrk_dau_truendau,     (Int_t)dauTruePart->Daughters.size());
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateGDaughterReco(OutputManager& output, AnaParticlePD* gdau, Int_t gdau_index){
//********************************************************************

  if(!gdau) return;
  
  output.FillMatrixVar           (seltrk_gdau_ndau,   (Int_t)gdau->Daughters.size(), -1, gdau_index   );
  output.Fill3DMatrixVarFromArray(seltrk_gdau_pos,           gdau->PositionStart,    -1, gdau_index, 4);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_dir,           gdau->DirectionStart,   -1, gdau_index, 3);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_endpos,        gdau->PositionEnd,      -1, gdau_index, 4);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_enddir,        gdau->DirectionEnd,     -1, gdau_index, 3);
  output.FillMatrixVar           (seltrk_gdau_length,        gdau->Length,           -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_nhits,  (Int_t)gdau->NHits,            -1, gdau_index   );

  if(gdau->Hits[2].empty()){
    for (int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.Fill3DMatrixVar(seltrk_gdau_hit_dedx,     (Float_t)-999., -1, gdau_index, j);
      output.Fill3DMatrixVar(seltrk_gdau_hit_resrange, (Float_t)-999., -1, gdau_index, j);
    }
  }
  else{
    int jmin = std::max<int>(0,(int)gdau->Hits[2].size()-NMAXHITSPERPLANE);
    int jmax = (int)gdau->Hits[2].size();
    for (int j = jmin; j < jmax; j++){
      output.Fill3DMatrixVar(seltrk_gdau_hit_dedx,      gdau->Hits[2][j].dEdx,          -1, gdau_index, j-jmin);
      output.Fill3DMatrixVar(seltrk_gdau_hit_resrange,  gdau->Hits[2][j].ResidualRange, -1, gdau_index, j-jmin);
    }
  }
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateGDaughterTrue(OutputManager& output, AnaParticlePD* gdau, Int_t gdau_index){
//********************************************************************

  if(!gdau) return;  
  AnaTrueParticle* gdauTruePart = static_cast<AnaTrueParticle*>(gdau->TrueObject);
  if(!gdauTruePart) return;
  
  output.FillMatrixVar           (seltrk_gdau_truepdg,  (Int_t)gdauTruePart->PDG,              -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_truendau, (Int_t)gdauTruePart->Daughters.size(), -1, gdau_index   );
  output.Fill3DMatrixVarFromArray(seltrk_gdau_truepos,         gdauTruePart->Position,         -1, gdau_index, 4);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_trueendpos,      gdauTruePart->PositionEnd,      -1, gdau_index, 4);
  output.FillMatrixVar           (seltrk_gdau_trueproc,        gdauTruePart->ProcessStart,     -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_trueendproc,     gdauTruePart->ProcessEnd,       -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_truemom,         gdauTruePart->Momentum,         -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_trueendmom,      gdauTruePart->MomentumEnd,      -1, gdau_index   );
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateGGDaughterReco(OutputManager& output, AnaParticlePD* ggdau, Int_t gdau_index, Int_t ggdau_index){
//********************************************************************

  if(!ggdau) return;
  
  output.Fill3DMatrixVar(seltrk_ggdau_ndau,  (Int_t)ggdau->Daughters.size(),  -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_posX,         ggdau->PositionStart[0],  -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_posY,         ggdau->PositionStart[1],  -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_posZ,         ggdau->PositionStart[2],  -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_dirX,         ggdau->DirectionStart[0], -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_dirY,         ggdau->DirectionStart[1], -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_dirZ,         ggdau->DirectionStart[2], -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_endposX,      ggdau->PositionEnd[0],    -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_endposY,      ggdau->PositionEnd[1],    -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_endposZ,      ggdau->PositionEnd[2],    -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_enddirX,      ggdau->DirectionEnd[0],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_enddirY,      ggdau->DirectionEnd[1],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_enddirZ,      ggdau->DirectionEnd[2],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_length,       ggdau->Length,            -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_nhits,        ggdau->NHits,             -1, gdau_index, ggdau_index);
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateGGDaughterTrue(OutputManager& output, AnaParticlePD* ggdau, Int_t gdau_index, Int_t ggdau_index){
//********************************************************************

  if(!ggdau) return;  
  AnaTrueParticle* ggdauTruePart = static_cast<AnaTrueParticle*>(ggdau->TrueObject);
  if(!ggdauTruePart) return;
  
  output.Fill3DMatrixVar(seltrk_ggdau_truepdg,  (Int_t)ggdauTruePart->PDG,              -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_truendau, (Int_t)ggdauTruePart->Daughters.size(), -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueposX,        ggdauTruePart->Position[0],      -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueposY,        ggdauTruePart->Position[1],      -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueposZ,        ggdauTruePart->Position[2],      -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueendposX,     ggdauTruePart->PositionEnd[0],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueendposY,     ggdauTruePart->PositionEnd[1],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueendposZ,     ggdauTruePart->PositionEnd[2],   -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueproc,        ggdauTruePart->ProcessStart,     -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueendproc,     ggdauTruePart->ProcessEnd,       -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_truemom,         ggdauTruePart->Momentum,         -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_trueendmom,      ggdauTruePart->MomentumEnd,      -1, gdau_index, ggdau_index);
}
