#include "kaonTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void kaonTree::AddKaonVariables_TrueKaonCandidates(OutputManager& output){
//********************************************************************

  AddVarF  (output, truekaon_truemom,        "true kaon candidate true momentum"     );
  AddVarF  (output, truekaon_trueendmom,     "true kaon candidate true end momentum" );
  AddVarI  (output, truekaon_truepdg,        "true kaon candidate true pdg"          );
  AddVarI  (output, truekaon_trueparentpdg,  "true kaon candidate parent true pdg"   );
  AddVarI  (output, truekaon_trueparentid,   "true kaon candidate parent true ID"    );
  AddVarI  (output, truekaon_trueproc,       "true kaon candidate true process"      );
  AddVarI  (output, truekaon_trueendproc,    "true kaon candidate true end process"  );
  AddVarI  (output, truekaon_truedecay,      "true kaon candidate true decay"        );
  AddVarI  (output, truekaon_truechainmuon,  "true kaon candidate chain to muon"     );
  AddVarI  (output, truekaon_truendau,       "true kaon candidate true ndaughters"   );
  AddVar4VF(output, truekaon_truepos,        "true kaon candidate true position"     );
  AddVar4VF(output, truekaon_trueendpos,     "true kaon candidate true end position" );
  AddVar3VF(output, truekaon_truedir,        "true kaon candidate true direction"    );
  AddVar3VF(output, truekaon_trueenddir,     "true kaon candidate true end direction");
  AddVarI  (output, truekaon_truegeneration, "true kaon candidate true generation"   );
  AddVarF  (output, truekaon_trueeff,        "true kaon candidate true efficciency"  );
  AddVarF  (output, truekaon_truepur,        "true kaon candidate true purity"       );
  AddVarI  (output, truekaon_branch,         "selection branch associated to this true kaon");

  AddVarF  (output, truekaon_truemuon_truemom,        "true kaon candidate muon daughter true momentum"     );
  AddVarF  (output, truekaon_truemuon_trueendmom,     "true kaon candidate muon daughter true end momentum" );
  AddVarI  (output, truekaon_truemuon_truepdg,        "true kaon candidate muon daughter true pdg"          );
  AddVarI  (output, truekaon_truemuon_trueproc,       "true kaon candidate muon daughter true process"      );
  AddVarI  (output, truekaon_truemuon_trueendproc,    "true kaon candidate muon daughter true end process"  );
  AddVarI  (output, truekaon_truemuon_truendau,       "true kaon candidate muon daughter true ndaughters"   );
  AddVar4VF(output, truekaon_truemuon_truepos,        "true kaon candidate muon daughter true position"     );
  AddVar4VF(output, truekaon_truemuon_trueendpos,     "true kaon candidate muon daughter true end position" );
  AddVar3VF(output, truekaon_truemuon_truedir,        "true kaon candidate muon daughter true direction"    );
  AddVar3VF(output, truekaon_truemuon_trueenddir,     "true kaon candidate muon daughter true end direction");
  AddVarI  (output, truekaon_truemuon_truegeneration, "true kaon candidate muon daughter true generation"   );
  AddVarF  (output, truekaon_truemuon_trueeff,        "true kaon candidate muon daughter true efficciency"  );
  AddVarF  (output, truekaon_truemuon_truepur,        "true kaon candidate muon daughter true purity"       );
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesReco(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, candidates_generation,   "candidates generation",        ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_parentID,     "candidates parent ID",         ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_ndau,         "candidates' daughters",        ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_nsisters,     "candidates' sistets",          ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_pos,          "candidates position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dir,          "candidates direction",         ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_endpos,       "candidates position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_enddir,       "candidates direction",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_length,       "candidates length",            ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_muon,     "candidates momentum (muon)",   ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_prot,     "candidates momentum (proton)", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_type,         "candidates object type",       ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_CNNscore,     "candidates CNN score",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_prot,    "candidates chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_kaon,    "candidates chi2 kaon",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_muon,    "candidates chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_ndf,     "candidates chi2 ndf",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_kaon_PID,     "candidates kaon chi2",         ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_kaon_PID_ndf, "candidates kaon ndf",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_mother, "candidates-mother distance",   ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_dau,    "candidates-daughter distance", ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_cos_dau,         "candidates-daughter cos",      ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_nhits,           "candidates #hits",             ncandidates, nmax);
 
  AddVarMaxSizeVF (output, candidates_averagedEdx,     "candidates average dEdx/hit",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_vtx_michelscore, "candidates michelscore in the vertex", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_vtx_nhits,       "candidates points in the vertex",      ncandidates, nmax);

  AddVarMaxSizeVI (output, candidates_dau_ndau,       "candidates daughter' daughters",        ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_pos,        "candidates daughter position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_dir,        "candidates daughter direction",         ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_endpos,     "candidates daughter position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_enddir,     "candidates daughter direction",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_length,     "candidates daughter length",            ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_muon,   "candidates daughter momentum (muon)",   ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_prot,   "candidates daughter momentum (proton)", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_type,       "candidates daughter object type",       ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_dau_CNNscore,   "candidates daughter CNN score",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_prot,  "candidates daughter chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_muon,  "candidates daughter chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_ndf,   "candidates daughter chi2 ndf",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_kaon_PID,   "candidates daughter chi2 ndf",          ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_kaon_PID_ndf, "candidates daughter chi2 ndf",          ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_nhits,      "candidates daughter #hits",             ncandidates, nmax);

  AddVarMaxSizeVF (output, candidates_dau_averagedEdx,     "candidates dau average dEdx",                  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_calE,            "candidates dau calorimetry energy deposition", ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_vtx_michelscore, "candidates dau michelscore in the vertex",     ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_vtx_nhits,       "candidates dau points in the vertex",          ncandidates, nmax);
  
  AddVarMaxSizeVI (output, candidates_gdau_ndau,       "candidates gdaughter' gdaughters",        ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_pos,        "candidates gdaughter position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_gdau_dir,        "candidates gdaughter direction",         ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_endpos,     "candidates gdaughter position",          ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_gdau_enddir,     "candidates gdaughter direction",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_length,     "candidates gdaughter length",            ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_mom_muon,   "candidates gdaughter momentum (muon)",   ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_mom_prot,   "candidates gdaughter momentum (proton)", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_type,       "candidates gdaughter object type",       ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_gdau_CNNscore,   "candidates gdaughter CNN score",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_prot,  "candidates gdaughter chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_muon,  "candidates gdaughter chi2 proton",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_ndf,   "candidates gdaughter chi2 ndf",          ncandidates, nmax);

  AddVarMaxSizeVF (output, candidates_gdau_averagedEdx,     "candidates gdau average dEdx",                  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_calE,            "candidates gdau calorimetry energy deposition", ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_vtx_michelscore, "candidates gdau michelscore in the vertex",     ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_vtx_nhits,       "candidates gdau points in the vertex",          ncandidates, nmax);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane){
//********************************************************************

  AddVarMF(output, candidates_hit_x,        "candidates x per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_y,        "candidates y per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_z,        "candidates z per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_dedx,     "candidates dEdx per hit",       ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_dedx_cal, "candidates dEdx cal per hit",   ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_resrange, "candidates hit residual range", ncandidates, -nmax, nmaxhitsperplane);

  AddVarMF(output, candidates_dau_hit_x,        "candidates daughter x per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_y,        "candidates daughter y per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_z,        "candidates daughter z per hit",          ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_dedx,     "candidates daughter dEdx per hit",       ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_dedx_cal, "candidates daughter dEdx cal per hit",   ncandidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_resrange, "candidates daughter hit residual range", ncandidates, -nmax, nmaxhitsperplane);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, candidates_truendau,       "candidates' true ndaughters",  ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_truegeneration, "candidates true generation",   ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_truepdg,        "candidates true pdg",          ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_trueorigin,     "candidates true origin",       ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_truepos,        "daugthers true position",      ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_trueendpos,     "daugthers true end position",  ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_trueproc,       "daugthers true process",       ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_trueendproc,    "candidates true end process",  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_truemom,        "candidates true momentum",     ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_trueendmom,     "candidates true end momentum", ncandidates, nmax);

  AddVarMaxSizeVI (output, candidates_dau_truendau,    "candidates daughter' true ndaughters",  ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_truepdg,     "candidates daughter true pdg",          ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_truepos,     "candidates daughter true position",     ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_trueendpos,  "candidates daughter true end position", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_trueproc,    "candidates daughter true process",      ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_trueendproc, "candidates daughter true end process",  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_truemom,     "candidates daughter true momentum",     ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_trueendmom,  "candidates daughter true end momentum", ncandidates, nmax);

  AddVarMaxSizeVI (output, candidates_gdau_truendau,    "candidates gdaughter' true ndaughters",  ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_truepdg,     "candidates gdaughter true pdg",          ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_truepos,     "candidates gdaughter true position",     ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_trueendpos,  "candidates gdaughter true end position", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_trueproc,    "candidates gdaughter true process",      ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_trueendproc, "candidates gdaughter true end process",  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_truemom,     "candidates gdaughter true momentum",     ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_trueendmom,  "candidates gdaughter true end momentum", ncandidates, nmax);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonBestCandidateReco(OutputManager& output){
//********************************************************************

  AddVarI  (output, bestcandidate_generation,   "bestcandidate generation"       );
  AddVarI  (output, bestcandidate_ndau,         "bestcandidate' daughters"       );
  AddVar4VF(output, bestcandidate_pos,          "bestcandidate position"         ); 
  AddVar3VF(output, bestcandidate_dir,          "bestcandidate direction"        );
  AddVar4VF(output, bestcandidate_endpos,       "bestcandidate position"         ); 
  AddVar3VF(output, bestcandidate_enddir,       "bestcandidate direction"        );
  AddVarF  (output, bestcandidate_length,       "bestcandidate length"           );
  AddVarF  (output, bestcandidate_mom_muon,     "bestcandidate momentum (muon)"  );
  AddVarF  (output, bestcandidate_mom_prot,     "bestcandidate momentum (proton)");
  AddVarI  (output, bestcandidate_type,         "bestcandidate object type"      );
  AddVar3VF(output, bestcandidate_CNNscore,     "bestcandidate CNN score"        );
  AddVarF  (output, bestcandidate_chi2_prot,    "bestcandidate chi2 proton"      );
  AddVarF  (output, bestcandidate_chi2_muon,    "bestcandidate chi2 proton"      );
  AddVarF  (output, bestcandidate_chi2_ndf,     "bestcandidate chi2 ndf"         );
  AddVarF  (output, bestcandidate_distance_mother, "bestcandidate-mother distance");
  AddVarF  (output, bestcandidate_distance_dau, "bestcandidate-daughter distance");
  AddVarF  (output, bestcandidate_cos_dau,      "bestcandidate-daughter cos"     );
  AddVarI  (output, bestcandidate_nhits,        "bestcandidate #hits"            );
 
  AddVarF  (output, bestcandidate_averagedEdx,     "bestcandidate average dEdx/hit"         );
  AddVarF  (output, bestcandidate_vtx_michelscore, "bestcandidate michelscore in the vertex");
  AddVarI  (output, bestcandidate_vtx_nhits,       "bestcandidate points in the vertex"     );

  AddVarI  (output, bestcandidate_dau_ndau,       "bestcandidate daughter' daughters"       );
  AddVar4VF(output, bestcandidate_dau_pos,        "bestcandidate daughter position"         ); 
  AddVar3VF(output, bestcandidate_dau_dir,        "bestcandidate daughter direction"        );
  AddVar4VF(output, bestcandidate_dau_endpos,     "bestcandidate daughter position"         ); 
  AddVar3VF(output, bestcandidate_dau_enddir,     "bestcandidate daughter direction"        );
  AddVarF  (output, bestcandidate_dau_length,     "bestcandidate daughter length"           );
  AddVarF  (output, bestcandidate_dau_mom_muon,   "bestcandidate daughter momentum (muon)"  );
  AddVarF  (output, bestcandidate_dau_mom_prot,   "bestcandidate daughter momentum (proton)");
  AddVarI  (output, bestcandidate_dau_type,       "bestcandidate daughter object type"      );
  AddVar3VF(output, bestcandidate_dau_CNNscore,   "bestcandidate daughter CNN score"        );
  AddVarF  (output, bestcandidate_dau_chi2_prot,  "bestcandidate daughter chi2 proton"      );
  AddVarF  (output, bestcandidate_dau_chi2_muon,  "bestcandidate daughter chi2 proton"      );
  AddVarF  (output, bestcandidate_dau_chi2_ndf,   "bestcandidate daughter chi2 ndf"         );
  AddVarI  (output, bestcandidate_dau_nhits,      "bestcandidate daughter #hits"            );

  AddVarF  (output, bestcandidate_dau_calE,            "bestcandidate dau calorimetry energy deposition");
  AddVarF  (output, bestcandidate_dau_vtx_michelscore, "bestcandidate dau michelscore in the vertex"    );
  AddVarI  (output, bestcandidate_dau_vtx_nhits,       "bestcandidate dau points in the vertex"         );
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonBestCandidateHitsReco(OutputManager& output, UInt_t nmaxhitsperplane){
//********************************************************************

  AddVarFixVF(output, bestcandidate_hit_x,        "bestcandidate x per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_hit_y,        "bestcandidate y per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_hit_z,        "bestcandidate z per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_hit_dedx,     "bestcandidate dEdx per hit",       nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_hit_dedx_cal, "bestcandidate dEdx cal per hit",   nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_hit_resrange, "bestcandidate hit residual range", nmaxhitsperplane);

  AddVarFixVF(output, bestcandidate_dau_hit_x,        "bestcandidate daughter x per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_dau_hit_y,        "bestcandidate daughter y per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_dau_hit_z,        "bestcandidate daughter z per hit",          nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_dau_hit_dedx,     "bestcandidate daughter dEdx per hit",       nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_dau_hit_dedx_cal, "bestcandidate daughter dEdx cal per hit",   nmaxhitsperplane);
  AddVarFixVF(output, bestcandidate_dau_hit_resrange, "bestcandidate daughter hit residual range", nmaxhitsperplane);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonBestCandidateTrue(OutputManager& output){
//********************************************************************

  AddVarI (output, bestcandidate_truendau,       "bestcandidate' true ndaughters" );
  AddVarI (output, bestcandidate_truegeneration, "bestcandidate true generation"  );
  AddVarI (output, bestcandidate_truepdg,        "bestcandidate true pdg"         );
  AddVar4VF(output, bestcandidate_truepos,        "daugthers true position"       );
  AddVar4VF(output, bestcandidate_trueendpos,     "daugthers true end position"   );
  AddVarI (output, bestcandidate_trueproc,       "daugthers true process"         );
  AddVarI (output, bestcandidate_trueendproc,    "bestcandidate true end process" );
  AddVarF (output, bestcandidate_truemom,        "bestcandidate true momentum"    );
  AddVarF (output, bestcandidate_trueendmom,     "bestcandidate true end momentum");

  AddVarI (output, bestcandidate_dau_truendau,    "bestcandidate daughter' true ndaughters"  );
  AddVarI (output, bestcandidate_dau_truepdg,     "bestcandidate daughter true pdg"          );
  AddVar4VF(output, bestcandidate_dau_truepos,     "bestcandidate daughter true position"    );
  AddVar4VF(output, bestcandidate_dau_trueendpos,  "bestcandidate daughter true end position");
  AddVarI (output, bestcandidate_dau_trueproc,    "bestcandidate daughter true process"      );
  AddVarI (output, bestcandidate_dau_trueendproc, "bestcandidate daughter true end process"  );
  AddVarF (output, bestcandidate_dau_truemom,     "bestcandidate daughter true momentum"     );
  AddVarF (output, bestcandidate_dau_trueendmom,  "bestcandidate daughter true end momentum" );
}

//********************************************************************
void kaonTree::FillKaonVariables_TrueKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart){
//********************************************************************

  if(!truePart)return;

  output.FillVar               (truekaon_truemom,            truePart->Momentum        );
  output.FillVar               (truekaon_trueendmom,         truePart->MomentumEnd     );
  output.FillVar               (truekaon_truepdg,            truePart->PDG             );
  output.FillVar               (truekaon_trueparentpdg,      truePart->ParentPDG       );
  output.FillVar               (truekaon_trueparentid,       truePart->ParentID        );
  output.FillVar               (truekaon_trueproc,           truePart->ProcessStart    ); 
  output.FillVar               (truekaon_trueendproc,        truePart->ProcessEnd      );
  //output.FillVar               (truekaon_truedecay,          kvtx.DecayMode            );
  //output.FillVar               (truekaon_truechainmuon,      kvtx.ChainMuon            );
  output.FillVar               (truekaon_truendau,    (Int_t)truePart->Daughters.size()); 
  output.FillVectorVarFromArray(truekaon_truepos,            truePart->Position,      4);
  output.FillVectorVarFromArray(truekaon_trueendpos,         truePart->PositionEnd,   4);
  output.FillVectorVarFromArray(truekaon_truedir,            truePart->Direction,     3);
  output.FillVectorVarFromArray(truekaon_trueenddir,         truePart->DirectionEnd,  3);
  //output.FillVar               (truekaon_branch,      (Int_t)kvtx.Branch               ); 
  
  /*if(kvtx.TrueParticlesVect.size()>1){
    AnaTrueParticlePD* trueDau = static_cast<AnaTrueParticlePD*>(kvtx.TrueParticlesVect[1]);
    if(trueDau){
      output.FillVar               (truekaon_truemuon_truemom,            trueDau->Momentum        );
      output.FillVar               (truekaon_truemuon_trueendmom,         trueDau->MomentumEnd     );
      output.FillVar               (truekaon_truemuon_truepdg,            trueDau->PDG             );
      output.FillVar               (truekaon_truemuon_truepdg,            trueDau->PDG             );
      output.FillVar               (truekaon_truemuon_trueproc,           trueDau->ProcessStart    ); 
      output.FillVar               (truekaon_truemuon_trueendproc,        trueDau->ProcessEnd      );
      output.FillVar               (truekaon_truemuon_truendau,    (Int_t)trueDau->Daughters.size()); 
      output.FillVectorVarFromArray(truekaon_truemuon_truepos,            trueDau->Position,      4);
      output.FillVectorVarFromArray(truekaon_truemuon_trueendpos,         trueDau->PositionEnd,   4);
      output.FillVectorVarFromArray(truekaon_truemuon_truedir,            trueDau->Direction,     3);
      output.FillVectorVarFromArray(truekaon_truemuon_trueenddir,         trueDau->DirectionEnd,  3);
    }
    }*/
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonCandidatesReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent){
//********************************************************************

  if (!part) return;
  output.FillVectorVar         (candidates_generation,       part->Generation       );
  output.FillVectorVar         (candidates_parentID,         part->ParentID      );
  output.FillVectorVar         (candidates_ndau,      (Int_t)part->Daughters.size() );
  if(parent)
    output.FillVectorVar         (candidates_nsisters,  (Int_t)parent->DaughtersIDs.size() );
  output.FillMatrixVarFromArray(candidates_pos,              part->PositionStart,  4);
  output.FillMatrixVarFromArray(candidates_dir,              part->DirectionStart, 3); 
  output.FillMatrixVarFromArray(candidates_endpos,           part->PositionEnd,    4);
  output.FillMatrixVarFromArray(candidates_enddir,           part->DirectionEnd,   3); 
  output.FillVectorVar         (candidates_length,           part->Length           );
  output.FillVectorVar         (candidates_mom_prot,         pdAnaUtils::ComputeRangeMomentum(part->Length,2212));
  output.FillVectorVar         (candidates_mom_muon,         pdAnaUtils::ComputeRangeMomentum(part->Length,13)  );
  output.FillVectorVar         (candidates_type,             part->Type             );
  output.FillMatrixVarFromArray(candidates_CNNscore,         part->CNNscore,       3); 
  output.FillVectorVar         (candidates_chi2_prot,        part->Chi2Proton       );
  output.FillVectorVar         (candidates_chi2_kaon,        (Float_t)pdAnaUtils::Chi2PID(*part,321).first);
  output.FillVectorVar         (candidates_chi2_muon,        part->Chi2Muon         );
  output.FillVectorVar         (candidates_chi2_ndf,         part->Chi2ndf          );
  std::pair<double,int>kaon_PID = pdAnaUtils::kaonPID(*part);
  output.FillVectorVar         (candidates_kaon_PID,         (Float_t)kaon_PID.first);
  output.FillVectorVar         (candidates_kaon_PID_ndf,     kaon_PID.second        );
  output.FillVectorVar         (candidates_nhits,            part->NHits            );

  output.FillVectorVar         (candidates_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );
  output.FillVectorVar         (candidates_vtx_michelscore,  part->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_vtx_nhits,        part->vtx_CNN_NHits);

  if(parent)output.FillVectorVar(candidates_distance_mother, pdAnaUtils::ComputeDistanceMotherDaughter(parent,part));

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if (!dau) return;
  output.FillVectorVar         (candidates_distance_dau,         pdAnaUtils::ComputeDistanceMotherDaughter(part,dau));
  output.FillVectorVar         (candidates_cos_dau,              pdAnaUtils::ComputeCosMotherDaughter(part,dau)     );
     
  output.FillVectorVar         (candidates_dau_ndau,      (Int_t)dau->Daughters.size() );
  output.FillMatrixVarFromArray(candidates_dau_pos,              dau->PositionStart,  4);
  output.FillMatrixVarFromArray(candidates_dau_dir,              dau->DirectionStart, 3); 
  output.FillMatrixVarFromArray(candidates_dau_endpos,           dau->PositionEnd,    4);
  output.FillMatrixVarFromArray(candidates_dau_enddir,           dau->DirectionEnd,   3); 
  output.FillVectorVar         (candidates_dau_length,           dau->Length           );
  output.FillVectorVar         (candidates_dau_mom_prot,         pdAnaUtils::ComputeRangeMomentum(dau->Length,2212));
  output.FillVectorVar         (candidates_dau_mom_muon,         pdAnaUtils::ComputeRangeMomentum(dau->Length,13)  );
  output.FillVectorVar         (candidates_dau_type,             dau->Type             );
  output.FillMatrixVarFromArray(candidates_dau_CNNscore,         dau->CNNscore,       3); 
  output.FillVectorVar         (candidates_dau_chi2_prot,        dau->Chi2Proton       );
  output.FillVectorVar         (candidates_dau_chi2_muon,        dau->Chi2Muon         );
  output.FillVectorVar         (candidates_dau_chi2_ndf,         dau->Chi2ndf          );
  std::pair<double,int>kaon_PID_dau = pdAnaUtils::kaonPID(*dau);
  output.FillVectorVar         (candidates_dau_kaon_PID,         (Float_t)kaon_PID_dau.first);
  output.FillVectorVar         (candidates_dau_kaon_PID_ndf,     kaon_PID_dau.second        );
  output.FillVectorVar         (candidates_dau_nhits,            dau->NHits            );

  output.FillVectorVar         (candidates_dau_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(dau,5));
  output.FillVectorVar         (candidates_dau_calE,             pdAnaUtils::ComputeKineticEnergy(*dau));
  output.FillVectorVar         (candidates_dau_vtx_michelscore,  dau->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_dau_vtx_nhits,        dau->vtx_CNN_NHits);

  if(dau->Daughters.empty())return;
  AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[0]);
  if (!gdau) return;
  output.FillVectorVar         (candidates_gdau_ndau,      (Int_t)gdau->Daughters.size() );
  output.FillMatrixVarFromArray(candidates_gdau_pos,              gdau->PositionStart,  4);
  output.FillMatrixVarFromArray(candidates_gdau_dir,              gdau->DirectionStart, 3); 
  output.FillMatrixVarFromArray(candidates_gdau_endpos,           gdau->PositionEnd,    4);
  output.FillMatrixVarFromArray(candidates_gdau_enddir,           gdau->DirectionEnd,   3); 
  output.FillVectorVar         (candidates_gdau_length,           gdau->Length           );
  output.FillVectorVar         (candidates_gdau_mom_prot,         pdAnaUtils::ComputeRangeMomentum(gdau->Length,2212));
  output.FillVectorVar         (candidates_gdau_mom_muon,         pdAnaUtils::ComputeRangeMomentum(gdau->Length,13)  );
  output.FillVectorVar         (candidates_gdau_type,             gdau->Type             );
  output.FillMatrixVarFromArray(candidates_gdau_CNNscore,         gdau->CNNscore,       3); 
  output.FillVectorVar         (candidates_gdau_chi2_prot,        gdau->Chi2Proton       );
  output.FillVectorVar         (candidates_gdau_chi2_muon,        gdau->Chi2Muon         );
  output.FillVectorVar         (candidates_gdau_chi2_ndf,         gdau->Chi2ndf          );

  output.FillVectorVar         (candidates_gdau_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(gdau,5));
  output.FillVectorVar         (candidates_gdau_calE,             pdAnaUtils::ComputeKineticEnergy(*gdau));
  output.FillVectorVar         (candidates_gdau_vtx_michelscore,  gdau->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_gdau_vtx_nhits,        gdau->vtx_CNN_NHits);

}

//********************************************************************
void kaonTree::FillKaonVariables_KaonCandidatesHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxhitsperplane){
//********************************************************************

  if (!part) return;
  if(part->Hits[2].empty()){
    for(int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(candidates_hit_x,            (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_hit_y,            (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_hit_z,            (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_hit_dedx,         (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_hit_dedx_cal,     (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_hit_resrange,     (Float_t)-999., -1, j);
    }
  }
  else{
    for (int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(candidates_hit_x,         (Float_t)part->Hits[2][j].Position.X(),   -1, j);
      output.FillMatrixVar(candidates_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),   -1, j);
      output.FillMatrixVar(candidates_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),   -1, j);
      output.FillMatrixVar(candidates_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,           -1, j);
      output.FillMatrixVar(candidates_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,     -1, j);
      output.FillMatrixVar(candidates_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,  -1, j);
    }
  }
  
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return;
  if(dau->Hits[2].empty()){
    for(int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(candidates_dau_hit_x,        (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_dau_hit_y,        (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_dau_hit_z,        (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_dau_hit_dedx,     (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_dau_hit_dedx_cal, (Float_t)-999., -1, j);
      output.FillMatrixVar(candidates_dau_hit_resrange, (Float_t)-999., -1, j);
    }
  }
  else{
    for(int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(candidates_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),   -1, j);
      output.FillMatrixVar(candidates_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),   -1, j);
      output.FillMatrixVar(candidates_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),   -1, j);
      output.FillMatrixVar(candidates_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,           -1, j);
      output.FillMatrixVar(candidates_dau_hit_dedx_cal,  (Float_t)dau->Hits[2][j].dEdx_calib,     -1, j);
      output.FillMatrixVar(candidates_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,  -1, j);
    }
  }
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonCandidatesTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;  
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart) return;

  output.FillVectorVar         (candidates_truendau,    (Int_t)truePart->Daughters.size());
  output.FillVectorVar         (candidates_truegeneration,     truePart->Generation      );
  output.FillVectorVar         (candidates_truepdg,            truePart->PDG             );
  output.FillVectorVar         (candidates_trueorigin,         truePart->Origin          );
  output.FillMatrixVarFromArray(candidates_truepos,            truePart->Position,      4);
  output.FillMatrixVarFromArray(candidates_trueendpos,         truePart->PositionEnd,   4);
  output.FillVectorVar         (candidates_trueproc,    (Int_t)truePart->ProcessStart    );
  output.FillVectorVar         (candidates_trueendproc, (Int_t)truePart->ProcessEnd      );
  output.FillVectorVar         (candidates_truemom,            truePart->Momentum        );
  output.FillVectorVar         (candidates_trueendmom,         truePart->MomentumEnd     );

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return;
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  if(!dauTruePart) return;

  output.FillVectorVar         (candidates_dau_truendau,    (Int_t)dauTruePart->Daughters.size());
  output.FillVectorVar         (candidates_dau_truepdg,            dauTruePart->PDG             );
  output.FillMatrixVarFromArray(candidates_dau_truepos,            dauTruePart->Position,      4);
  output.FillMatrixVarFromArray(candidates_dau_trueendpos,         dauTruePart->PositionEnd,   4);
  output.FillVectorVar         (candidates_dau_trueproc,    (Int_t)dauTruePart->ProcessStart    );
  output.FillVectorVar         (candidates_dau_trueendproc, (Int_t)dauTruePart->ProcessEnd      );
  output.FillVectorVar         (candidates_dau_truemom,            dauTruePart->Momentum        );
  output.FillVectorVar         (candidates_dau_trueendmom,         dauTruePart->MomentumEnd     );

  if(dau->Daughters.empty())return;
  AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[0]);
  if(!gdau)return;
  AnaTrueParticle* gdauTruePart = static_cast<AnaTrueParticle*>(gdau->TrueObject);
  if(!gdauTruePart) return;

  output.FillVectorVar         (candidates_gdau_truendau,    (Int_t)gdauTruePart->Daughters.size());
  output.FillVectorVar         (candidates_gdau_truepdg,            gdauTruePart->PDG             );
  output.FillMatrixVarFromArray(candidates_gdau_truepos,            gdauTruePart->Position,      4);
  output.FillMatrixVarFromArray(candidates_gdau_trueendpos,         gdauTruePart->PositionEnd,   4);
  output.FillVectorVar         (candidates_gdau_trueproc,    (Int_t)gdauTruePart->ProcessStart    );
  output.FillVectorVar         (candidates_gdau_trueendproc, (Int_t)gdauTruePart->ProcessEnd      );
  output.FillVectorVar         (candidates_gdau_truemom,            gdauTruePart->Momentum        );
  output.FillVectorVar         (candidates_gdau_trueendmom,         gdauTruePart->MomentumEnd     );
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonBestCandidateReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent){
//********************************************************************

  if (!part) return;
  output.FillVar               (bestcandidate_generation,       part->Generation       );
  output.FillVar               (bestcandidate_ndau,      (Int_t)part->Daughters.size() );
  output.FillVectorVarFromArray(bestcandidate_pos,              part->PositionStart,  4);
  output.FillVectorVarFromArray(bestcandidate_dir,              part->DirectionStart, 3); 
  output.FillVectorVarFromArray(bestcandidate_endpos,           part->PositionEnd,    4);
  output.FillVectorVarFromArray(bestcandidate_enddir,           part->DirectionEnd,   3); 
  output.FillVar               (bestcandidate_length,           part->Length           );
  output.FillVar               (bestcandidate_mom_prot,         pdAnaUtils::ComputeRangeMomentum(part->Length,2212));
  output.FillVar               (bestcandidate_mom_muon,         pdAnaUtils::ComputeRangeMomentum(part->Length,13)  );
  output.FillVar               (bestcandidate_type,             part->Type             );
  output.FillVectorVarFromArray(bestcandidate_CNNscore,         part->CNNscore,       3); 
  output.FillVar               (bestcandidate_chi2_prot,        part->Chi2Proton       );
  output.FillVar               (bestcandidate_chi2_muon,        part->Chi2Muon         );
  output.FillVar               (bestcandidate_chi2_ndf,         part->Chi2ndf          );
  output.FillVar               (bestcandidate_nhits,            part->NHits            );

  output.FillVar               (bestcandidate_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );
  output.FillVar               (bestcandidate_vtx_michelscore,  part->vtx_CNN_michelscore);
  output.FillVar               (bestcandidate_vtx_nhits,        part->vtx_CNN_NHits);

  if(parent)output.FillVar(bestcandidate_distance_mother,         pdAnaUtils::ComputeDistanceMotherDaughter(parent,part));

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if (!dau) return;
  output.FillVar               (bestcandidate_distance_dau,         pdAnaUtils::ComputeDistanceMotherDaughter(part,dau));
  output.FillVar               (bestcandidate_cos_dau,              pdAnaUtils::ComputeCosMotherDaughter(part,dau)     );

  output.FillVar               (bestcandidate_dau_ndau,      (Int_t)dau->Daughters.size() );
  output.FillVectorVarFromArray(bestcandidate_dau_pos,              dau->PositionStart,  4);
  output.FillVectorVarFromArray(bestcandidate_dau_dir,              dau->DirectionStart, 3); 
  output.FillVectorVarFromArray(bestcandidate_dau_endpos,           dau->PositionEnd,    4);
  output.FillVectorVarFromArray(bestcandidate_dau_enddir,           dau->DirectionEnd,   3); 
  output.FillVar               (bestcandidate_dau_length,           dau->Length           );
  output.FillVar               (bestcandidate_dau_mom_prot,         pdAnaUtils::ComputeRangeMomentum(dau->Length,2212));
  output.FillVar               (bestcandidate_dau_mom_muon,         pdAnaUtils::ComputeRangeMomentum(dau->Length,13)  );
  output.FillVar               (bestcandidate_dau_type,             dau->Type             );
  output.FillVectorVarFromArray(bestcandidate_dau_CNNscore,         dau->CNNscore,       3); 
  output.FillVar               (bestcandidate_dau_chi2_prot,        dau->Chi2Proton       );
  output.FillVar               (bestcandidate_dau_chi2_muon,        dau->Chi2Muon         );
  output.FillVar               (bestcandidate_dau_chi2_ndf,         dau->Chi2ndf          );
  output.FillVar               (bestcandidate_dau_nhits,            dau->NHits            );

  output.FillVar               (bestcandidate_dau_calE,             pdAnaUtils::ComputeKineticEnergy(*part));
  output.FillVar               (bestcandidate_dau_vtx_michelscore,  dau->vtx_CNN_michelscore);
  output.FillVar               (bestcandidate_dau_vtx_nhits,        dau->vtx_CNN_NHits);

}

//********************************************************************
void kaonTree::FillKaonVariables_KaonBestCandidateHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxhitsperplane){
//********************************************************************

  if (!part) return;
  if(part->Hits[2].empty()){
    for(int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillVectorVar(bestcandidate_hit_x,            (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_hit_y,            (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_hit_z,            (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_hit_dedx,         (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_hit_dedx_cal,     (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_hit_resrange,     (Float_t)-999.,  j);
    }
  }
  else{
    for (int j = 0; j < (int)nmaxhitsperplane; j++){
      output.FillVectorVar(bestcandidate_hit_x,         (Float_t)part->Hits[2][j].Position.X(),    j);
      output.FillVectorVar(bestcandidate_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),    j);
      output.FillVectorVar(bestcandidate_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),    j);
      output.FillVectorVar(bestcandidate_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,            j);
      output.FillVectorVar(bestcandidate_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,      j);
      output.FillVectorVar(bestcandidate_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,   j);
    }
  }
  
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return;
  if(dau->Hits[2].empty()){
    for(int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillVectorVar(bestcandidate_dau_hit_x,        (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_dau_hit_y,        (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_dau_hit_z,        (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_dau_hit_dedx,     (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_dau_hit_dedx_cal, (Float_t)-999.,  j);
      output.FillVectorVar(bestcandidate_dau_hit_resrange, (Float_t)-999.,  j);
    }
  }
  else{
    for (int j = 0; j < (int)nmaxhitsperplane; j++){
      output.FillVectorVar(bestcandidate_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),    j);
      output.FillVectorVar(bestcandidate_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),    j);
      output.FillVectorVar(bestcandidate_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),    j);
      output.FillVectorVar(bestcandidate_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,            j);
      output.FillVectorVar(bestcandidate_dau_hit_dedx_cal,  (Float_t)dau->Hits[2][j].dEdx_calib,      j);
      output.FillVectorVar(bestcandidate_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,   j);
    }
  }
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonBestCandidateTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;  
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart) return;

  output.FillVar               (bestcandidate_truendau,    (Int_t)truePart->Daughters.size());
  output.FillVar               (bestcandidate_truegeneration,     truePart->Generation      );
  output.FillVar               (bestcandidate_truepdg,            truePart->PDG             );
  output.FillVectorVarFromArray(bestcandidate_truepos,            truePart->Position,      4);
  output.FillVectorVarFromArray(bestcandidate_trueendpos,         truePart->PositionEnd,   4);
  output.FillVar               (bestcandidate_trueproc,    (Int_t)truePart->ProcessStart    );
  output.FillVar               (bestcandidate_trueendproc, (Int_t)truePart->ProcessEnd      );
  output.FillVar               (bestcandidate_truemom,            truePart->Momentum        );
  output.FillVar               (bestcandidate_trueendmom,         truePart->MomentumEnd     );

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if(!dau)return;
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  if(!dauTruePart) return;

  output.FillVar               (bestcandidate_dau_truendau,    (Int_t)dauTruePart->Daughters.size());
  output.FillVar               (bestcandidate_dau_truepdg,            dauTruePart->PDG             );
  output.FillVectorVarFromArray(bestcandidate_dau_truepos,            dauTruePart->Position,      4);
  output.FillVectorVarFromArray(bestcandidate_dau_trueendpos,         dauTruePart->PositionEnd,   4);
  output.FillVar               (bestcandidate_dau_trueproc,    (Int_t)dauTruePart->ProcessStart    );
  output.FillVar               (bestcandidate_dau_trueendproc, (Int_t)dauTruePart->ProcessEnd      );
  output.FillVar               (bestcandidate_dau_truemom,            dauTruePart->Momentum        );
  output.FillVar               (bestcandidate_dau_trueendmom,         dauTruePart->MomentumEnd     );
}
