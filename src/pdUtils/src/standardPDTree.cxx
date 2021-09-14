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
  //AddVarF(    output, seltrk_length_alt,     "candidate alternate length");
  //AddVarF(    output, seltrk_length_raw,     "candidate length");
  //AddVarFixVI(output, seltrk_nhitsperplane,  "candidate number of hits per plane",3);

  AddVarF(    output, seltrk_mom_muon,       "candidate momentum muon");
  AddVarF(    output, seltrk_mom_prot,       "candidate momentum proton");
  //AddVarF(    output, seltrk_mom_muon_alt,   "candidate alternate momentum muon");
  //AddVarF(    output, seltrk_mom_prot_alt,   "candidate alternate momentum proton");
  AddVarF(    output, seltrk_dqdx,           "candidate average dQdx");
  AddVarF(    output, seltrk_dedx,           "candidate average dEdx");
  AddVarF(    output, seltrk_dedx_raw,       "candidate average dEdx");
  AddVar3VF(  output, seltrk_CNNscore,       "candidate reconstructed CNN score");
  AddVarF(    output, seltrk_chi2_prot,      "candidate chi2 prot");
  AddVarF(    output, seltrk_chi2_muon,      "candidate chi2 muon");
  AddVarF(    output, seltrk_chi2_ndf,       "candidate chi2 ndf");
  //AddVarFixMF(output, seltrk_pid,            "candidate PID variables", 3,8);
  //AddVarFixMF(output, seltrk_calo,           "candidate CALO variables",3,5); 
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
  AddVarFixMI(output, seltrk_hit_ch,         "candidate channel per hit",3,NMAXHITSPERPLANE_SELTRK);
  AddVarFixMI(output, seltrk_hit_t0,         "candidate t0 per hit",3,NMAXHITSPERPLANE_SELTRK);

  AddVarFixMF(output, seltrk_hit_cnn,        "candidate cnn for hit",NMAXHITSPERPLANE_SELTRK,3);

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

  AddVarMaxSizeVI (output, trk_truendau,       "#true ndaughters",   ntracks, nmax);
  AddVarMaxSizeVI (output, trk_truegeneration, "true generation",    ntracks, nmax);
  AddVarMaxSizeVI (output, trk_truepdg,        "true pdg",           ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueorigin,     "true origin",        ntracks, nmax);
  AddVarMaxSize4MF(output, trk_truepos,        "true position",      ntracks, nmax);
  AddVarMaxSize4MF(output, trk_trueendpos,     "true end position",  ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueproc,       "true process",ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueendproc,    "true end process",   ntracks, nmax);
  AddVarMaxSizeVF (output, trk_truemom,        "true momentum",      ntracks, nmax);
  AddVarMaxSizeVF (output, trk_trueendmom,     "true end momentum",  ntracks, nmax);

}

//********************************************************************
void standardPDTree::AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, trk_ndau,         "#daughters",        ntracks, nmax);
  AddVarMaxSize4MF(output, trk_pos,          "position",          ntracks, nmax); 
  AddVarMaxSize3MF(output, trk_dir,          "direction",         ntracks, nmax);
  AddVarMaxSize4MF(output, trk_endpos,       "position",          ntracks, nmax); 
  AddVarMaxSize3MF(output, trk_enddir,       "direction",         ntracks, nmax);
  AddVarMaxSizeVF (output, trk_length,       "length",            ntracks, nmax);
  AddVarMaxSizeVF (output, trk_mom_muon,     "momentum (muon)",   ntracks, nmax);
  AddVarMaxSizeVF (output, trk_mom_prot,     "momentum (proton)", ntracks, nmax);
  AddVarMaxSizeVI (output, trk_type,         "object type",       ntracks, nmax);
  AddVarMaxSize3MF(output, trk_CNNscore,     "CNN score",         ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_prot,    "chi2 proton",       ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_muon,    "chi2 proton",       ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_ndf,     "chi2 ndf",          ntracks, nmax);
  AddVarMaxSizeVF (output, trk_averagedEdx,  "average dEdx/hit",  ntracks, nmax);
  AddVarMaxSizeVI (output, trk_nhits,        "#hits",             ntracks, nmax);

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
  
  AddVarMaxSizeVF(output,  seltrk_dau_mom_muon,   "daughters momentum (muon)",       seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_mom_prot,   "daughters momentum (proton)",     seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_pos,        "daughters position",              seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_dir,        "daughters direction",             seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_endpos,     "daughters position",              seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_enddir,     "daughters direction",             seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_length,     "daughters length",                seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_nhits,      "daughters #hits",                 seltrk_ndau,nmax);
  //AddVarMaxSizeVI( output, seltrk_dau_nhits2,     "daughters hits in plane 2",       seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_ndau,       "daughters' daughters"     ,       seltrk_ndau,nmax);
  AddVarMaxSize3MF(output, seltrk_dau_CNNscore,   "candidate daughters CNN score",   seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_prot,  "candidate daughters chi2 proton", seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_muon,  "candidate daughters chi2 proton", seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_ndf,   "candidate daughters chi2 ndf",    seltrk_ndau,nmax);

}

//********************************************************************
void standardPDTree::AddStandardVariables_CandidateDaughtersHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane){
//********************************************************************
  
  AddVarMF(output, seltrk_dau_hit_x,        "daughters x per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_y,        "daughters y per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_z,        "daughters z per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);

  AddVarMF(output, seltrk_dau_hit_dedx,     "daughters dEdx per hit",       seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_dedx_cal, "daughters dEdx per hit calib", seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_dqdx_raw, "daughters raw dQdx per hit",   seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_resrange, "daughters hit residual range", seltrk_ndau,-nmax,nmaxhitsperplane);

  AddVarMF(output, seltrk_dau_hit_cnn,      "daughters hit cnn",            seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMI(output, seltrk_dau_hit_ch,       "daughters hit channel",        seltrk_ndau,-nmax,nmaxhitsperplane);

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
  output.FillVectorVar         (trk_generation,       part->Generation       );
  output.FillVectorVar         (trk_ndau,      (Int_t)part->Daughters.size() );
  output.FillMatrixVarFromArray(trk_pos,              part->PositionStart,  4);
  output.FillMatrixVarFromArray(trk_dir,              part->DirectionStart, 3); 
  output.FillMatrixVarFromArray(trk_endpos,           part->PositionEnd,    4);
  output.FillMatrixVarFromArray(trk_enddir,           part->DirectionEnd,   3); 
  output.FillVectorVar         (trk_length,           part->Length           );
  output.FillVectorVar         (trk_mom_prot,         pdAnaUtils::ComputeRangeMomentum(part->Length,2212));
  output.FillVectorVar         (trk_mom_muon,         pdAnaUtils::ComputeRangeMomentum(part->Length,13)  );
  output.FillVectorVar         (trk_type,             part->Type             );
  output.FillMatrixVarFromArray(trk_CNNscore,         part->CNNscore,       3); 
  output.FillVectorVar         (trk_chi2_prot,        part->Chi2Proton       );
  output.FillVectorVar         (trk_chi2_muon,        part->Chi2Muon         );
  output.FillVectorVar         (trk_chi2_ndf,         part->Chi2ndf          );
  output.FillVectorVar         (trk_nhits,            part->NHits            );
  output.FillVectorVar         (trk_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );

}
  
//********************************************************************
void standardPDTree::FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  if (!part->TrueObject) return;
 
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart) return;

  output.FillVectorVar         (trk_truendau,    (Int_t)truePart->Daughters.size());
  output.FillVectorVar         (trk_truegeneration,     truePart->Generation      );
  output.FillVectorVar         (trk_truepdg,            truePart->PDG             );
  output.FillVectorVar         (trk_trueorigin,         truePart->Origin          );
  output.FillMatrixVarFromArray(trk_truepos,            truePart->Position,      4);
  output.FillMatrixVarFromArray(trk_trueendpos,         truePart->PositionEnd,   4);
  output.FillVectorVar         (trk_trueproc,    (Int_t)truePart->ProcessStart    );
  output.FillVectorVar         (trk_trueendproc, (Int_t)truePart->ProcessEnd      );
  output.FillVectorVar         (trk_truemom,            truePart->Momentum        );
  output.FillVectorVar         (trk_trueendmom,         truePart->MomentumEnd     );
  
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

  output.FillVar(seltrk_mom_prot,                (Float_t)part->RangeMomentum[0]);
  output.FillVar(seltrk_mom_muon,                (Float_t)part->RangeMomentum[1]);
  //output.FillVar(seltrk_mom_prot_alt,          (Float_t)part->RangeMomentum_alt[0]);
  //output.FillVar(seltrk_mom_muon_alt,          (Float_t)part->RangeMomentum_alt[1]);

  output.FillVar(seltrk_dedx,                    part->AveragedEdx);
  output.FillVar(seltrk_dqdx,                    part->AveragedQdx);

  output.FillVar(seltrk_dedx_raw,                static_cast<const AnaParticlePD*>(part->Original)->AveragedEdx);
  output.FillVar(seltrk_nhits,                   part->NHits);
  //output.FillVar(seltrk_length_alt,            part->Length_alt);
  output.FillVar(seltrk_length,                  pdAnaUtils::ComputeTrackLengthFromHitPosition(part));
  //output.FillVar(seltrk_length_raw,            part->Original->Length);
  output.FillVar(seltrk_costheta,                part->DirectionStart[2]);
  output.FillVectorVarFromArray(seltrk_pos,      part->PositionStart, 4);
  output.FillVectorVarFromArray(seltrk_endpos,   part->PositionEnd, 4);
  output.FillVectorVarFromArray(seltrk_dir,      part->DirectionStart, 3);
  output.FillVectorVarFromArray(seltrk_enddir,   part->DirectionEnd, 3);
  output.FillVectorVarFromArray(seltrk_CNNscore, part->CNNscore,3);
  output.FillVar(seltrk_chi2_prot,               part->Chi2Proton);
  output.FillVar(seltrk_chi2_muon,               part->Chi2Muon);
  output.FillVar(seltrk_chi2_ndf,                part->Chi2ndf);

  //for (int i=0;i<3;i++){
  //  output.FillMatrixVarFromArray(seltrk_pid,    part->PID[i],  i, 8);
  //  output.FillMatrixVarFromArray(seltrk_calo,   part->CALO[i], i, 5);      
  //}
}

//********************************************************************
void standardPDTree::FillStandardVariables_CandidateHitsReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  for (int i = 0; i < 3; i++){
    if(part->Hits[i].empty()){
      for (int j = 0; j < (int)NMAXHITSPERPLANE_SELTRK; j++){
        output.FillMatrixVar(seltrk_hit_ch      ,  (Int_t)-999, i, j);
        output.FillMatrixVar(seltrk_hit_t0      ,  (Int_t)-999, i, j);
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
        output.FillMatrixVar(seltrk_hit_ch   ,  (Int_t)part->Hits[i][j].Channel,      i, j-jmin);
        output.FillMatrixVar(seltrk_hit_t0   ,  (Int_t)part->Hits[i][j].StartTick,    i, j-jmin);
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
    //output.FillVectorVarFromArray(seltrk_nhitsperplane, part->NHitsPerPlane,3);
  }


  int jmin = std::max<int>(0,(int)part->Hits[2].size()-NMAXHITSPERPLANE_SELTRK);
  int jmax = (int)part->Hits[2].size();

  if(!part->Hits[2].empty()){
    for (int j = jmin; j < jmax; j++){
      for (int i = jmin; i < 3; i++){
	output.FillMatrixVar(seltrk_hit_cnn,       part->Hits[2][j].CNN[i], j-jmin, i);
      }
    }
  }
}
  
//********************************************************************
void standardPDTree::FillStandardVariables_CandidateDaughterReco(OutputManager& output, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;
  output.FillVectorVar         (seltrk_dau_mom_prot,         dau->RangeMomentum[0]);
  output.FillVectorVar         (seltrk_dau_mom_muon,         dau->RangeMomentum[1]);
  output.FillMatrixVarFromArray(seltrk_dau_pos,              dau->PositionStart, 4);
  output.FillMatrixVarFromArray(seltrk_dau_dir,              dau->DirectionStart,3); 
  output.FillMatrixVarFromArray(seltrk_dau_endpos,           dau->PositionEnd,   4);
  output.FillMatrixVarFromArray(seltrk_dau_enddir,           dau->DirectionEnd,  3); 
  output.FillVectorVar         (seltrk_dau_length,           dau->Length);
  //output.FillVectorVar         (seltrk_dau_nhits2,           dau->NHitsPerPlane[2]);
  output.FillVectorVar         (seltrk_dau_nhits,            dau->NHits           );
  output.FillVectorVar         (seltrk_dau_ndau,      (Int_t)dau->Daughters.size());
  output.FillMatrixVarFromArray(seltrk_dau_CNNscore,         dau->CNNscore,      3); 
  output.FillVectorVar         (seltrk_dau_chi2_prot,        dau->Chi2Proton      );
  output.FillVectorVar         (seltrk_dau_chi2_muon,        dau->Chi2Muon        );
  output.FillVectorVar         (seltrk_dau_chi2_ndf,         dau->Chi2ndf         );

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
void standardPDTree::FillStandardVariables_CandidateDaughterHitsReco(OutputManager& output, AnaParticlePD* dau, UInt_t nmaxhitsperplane){
//********************************************************************

  if (!dau) return;
  if(dau->Hits[2].empty()){
    for (int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(seltrk_dau_hit_x,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_y,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_z,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_dedx,     (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_dqdx_raw, (Float_t)-999., -1, j);  
      output.FillMatrixVar(seltrk_dau_hit_resrange, (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_cnn,      (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_ch,       (Int_t)-999   , -1, j);
    }
  }
  else{
    int jmin = std::max<int>(0,(int)dau->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)dau->Hits[2].size(); 

    for (int j = jmin; j < jmax; j++){
      output.FillMatrixVar(seltrk_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,           -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_dedx_cal,  (Float_t)dau->Hits[2][j].dEdx_calib,     -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,  -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_cnn,       (Float_t)dau->Hits[2][j].CNN[1],         -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_ch,        (Int_t)  dau->Hits[2][j].Channel,        -1, j-jmin);
    }
  }
}
