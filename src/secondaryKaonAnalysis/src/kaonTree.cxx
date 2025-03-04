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

  AddVarMaxSizeVI (output, candidates_generation,      "candidates generation",                ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_parentID,        "candidates parent ID",                 ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_ndau,            "candidates' daughters",                ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_nsisters,        "candidates' sistets",                  ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_pos,             "candidates position",                  ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dir,             "candidates direction",                 ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_endpos,          "candidates position",                  ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_enddir,          "candidates direction",                 ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_length,          "candidates length",                    ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_muon,        "candidates momentum (muon)",           ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_prot,        "candidates momentum (proton)",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_calE,            "candidates calorimetric energy",       ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_type,            "candidates object type",               ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_CNNscore,        "candidates CNN score",                 ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_prot,       "candidates chi2 proton",               ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_kaon,       "candidates chi2 kaon",                 ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_muon,       "candidates chi2 proton",               ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_ndf,        "candidates chi2 ndf",                  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_mother, "candidates-mother distance",           ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_dau,    "candidates-daughter distance",         ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_tcp,    "candidates dist to closest part",      ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_cos_dau,         "candidates-daughter cos",              ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_nhits,           "candidates #hits",                     ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_average_dedx,    "candidates average dEdx/hit",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_truncated_dedx,  "candidates truncated dEdx/hit",        ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_vtx_michelscore, "candidates michelscore in the vertex", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_vtx_nhits,       "candidates points in the vertex",      ncandidates, nmax);

  AddVarMaxSizeVI (output, candidates_dau_ndau,           "candidates daughter' daughters",            ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_pos,            "candidates daughter position",              ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_dir,            "candidates daughter direction",             ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_endpos,         "candidates daughter position",              ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_enddir,         "candidates daughter direction",             ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_length,         "candidates daughter length",                ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_muon,       "candidates daughter momentum (muon)",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_prot,       "candidates daughter momentum (proton)",     ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_type,           "candidates daughter object type",           ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_dau_CNNscore,       "candidates daughter CNN score",             ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_prot,      "candidates daughter chi2 proton",           ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_muon,      "candidates daughter chi2 proton",           ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_ndf,       "candidates daughter chi2 ndf",              ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_nhits,          "candidates daughter #hits",                 ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_forced,         "candidates forced daughter",                ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_forced_matched, "candidates forced daughter matched",        ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_average_dedx,    "candidates dau average dEdx",              ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_truncated_dedx,  "candidates dau truncated dEdx",            ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_calE,            "candidates dau calorimetric energy",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_vtx_michelscore, "candidates dau michelscore in the vertex", ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_vtx_nhits,       "candidates dau points in the vertex",      ncandidates, nmax);
  
  AddVarMaxSizeVI (output, candidates_gdau_ndau,            "candidates gdaughter' gdaughters",              ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_pos,             "candidates gdaughter position",                 ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_gdau_dir,             "candidates gdaughter direction",                ncandidates, nmax);
  AddVarMaxSize4MF(output, candidates_gdau_endpos,          "candidates gdaughter position",                 ncandidates, nmax); 
  AddVarMaxSize3MF(output, candidates_gdau_enddir,          "candidates gdaughter direction",                ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_length,          "candidates gdaughter length",                   ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_mom_muon,        "candidates gdaughter momentum (muon)",          ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_mom_prot,        "candidates gdaughter momentum (proton)",        ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_type,            "candidates gdaughter object type",              ncandidates, nmax);
  AddVarMaxSize3MF(output, candidates_gdau_CNNscore,        "candidates gdaughter CNN score",                ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_prot,       "candidates gdaughter chi2 proton",              ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_muon,       "candidates gdaughter chi2 proton",              ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_chi2_ndf,        "candidates gdaughter chi2 ndf",                 ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_average_dedx,    "candidates gdau average dEdx",                  ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_calE,            "candidates gdau calorimetry energy deposition", ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_gdau_vtx_michelscore, "candidates gdau michelscore in the vertex",     ncandidates, nmax);
  AddVarMaxSizeVI (output, candidates_gdau_vtx_nhits,       "candidates gdau points in the vertex",          ncandidates, nmax);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesRecoPID(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVF (output, candidates_lkl_prot,     "candidates proton likelihood",     ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_lkl_kaon,     "candidates kaon likelihood",       ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_lkl_prot, "candidates dau proton likelihood", ncandidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_lkl_muon, "candidates dau muon likelihood",   ncandidates, nmax);
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

  AddVarI  (output, bestcandidate_generation,      "bestcandidate generation"       );
  AddVarI  (output, bestcandidate_ndau,            "bestcandidate' daughters"       );
  AddVar4VF(output, bestcandidate_pos,             "bestcandidate position"         ); 
  AddVar3VF(output, bestcandidate_dir,             "bestcandidate direction"        );
  AddVar4VF(output, bestcandidate_endpos,          "bestcandidate position"         ); 
  AddVar3VF(output, bestcandidate_enddir,          "bestcandidate direction"        );
  AddVarF  (output, bestcandidate_length,          "bestcandidate length"           );
  AddVarF  (output, bestcandidate_mom_muon,        "bestcandidate momentum (muon)"  );
  AddVarF  (output, bestcandidate_mom_prot,        "bestcandidate momentum (proton)");
  AddVarI  (output, bestcandidate_type,            "bestcandidate object type"      );
  AddVar3VF(output, bestcandidate_CNNscore,        "bestcandidate CNN score"        );
  AddVarF  (output, bestcandidate_chi2_prot,       "bestcandidate chi2 proton"      );
  AddVarF  (output, bestcandidate_chi2_muon,       "bestcandidate chi2 proton"      );
  AddVarF  (output, bestcandidate_chi2_kaon,       "bestcandidate chi2 kaon"      );
  AddVarF  (output, bestcandidate_chi2_ndf,        "bestcandidate chi2 ndf"         );
  AddVarF  (output, bestcandidate_distance_mother, "bestcandidate-mother distance");
  AddVarF  (output, bestcandidate_distance_dau,    "bestcandidate-daughter distance");
  AddVarF  (output, bestcandidate_cos_dau,         "bestcandidate-daughter cos"     );
  AddVarI  (output, bestcandidate_nhits,           "bestcandidate #hits"            );
  AddVarF  (output, bestcandidate_average_dedx,    "bestcandidate average dEdx/hit"         );
  AddVarF  (output, bestcandidate_vtx_michelscore, "bestcandidate michelscore in the vertex");
  AddVarI  (output, bestcandidate_vtx_nhits,       "bestcandidate points in the vertex"     );
  AddVarF  (output, bestcandidate_calE,            "bestcandidate calorimetry energy deposition"     );
  AddVarF  (output, bestcandidate_truncated_dedx,            "bestcandidate calorimetry energy deposition"     );

  AddVarI  (output, bestcandidate_dau_ndau,            "bestcandidate daughter' daughters"       );
  AddVar4VF(output, bestcandidate_dau_pos,             "bestcandidate daughter position"         ); 
  AddVar3VF(output, bestcandidate_dau_dir,             "bestcandidate daughter direction"        );
  AddVar4VF(output, bestcandidate_dau_endpos,          "bestcandidate daughter position"         ); 
  AddVar3VF(output, bestcandidate_dau_enddir,          "bestcandidate daughter direction"        );
  AddVarF  (output, bestcandidate_dau_length,          "bestcandidate daughter length"           );
  AddVarF  (output, bestcandidate_dau_mom_muon,        "bestcandidate daughter momentum (muon)"  );
  AddVarF  (output, bestcandidate_dau_mom_prot,        "bestcandidate daughter momentum (proton)");
  AddVarI  (output, bestcandidate_dau_type,            "bestcandidate daughter object type"      );
  AddVar3VF(output, bestcandidate_dau_CNNscore,        "bestcandidate daughter CNN score"        );
  AddVarF  (output, bestcandidate_dau_chi2_prot,       "bestcandidate daughter chi2 proton"      );
  AddVarF  (output, bestcandidate_dau_chi2_muon,       "bestcandidate daughter chi2 proton"      );
  AddVarF  (output, bestcandidate_dau_chi2_ndf,        "bestcandidate daughter chi2 ndf"         );
  AddVarI  (output, bestcandidate_dau_nhits,           "bestcandidate daughter #hits"            );
  AddVarF  (output, bestcandidate_dau_calE,            "bestcandidate dau calorimetry energy deposition");
  AddVarF  (output, bestcandidate_dau_vtx_michelscore, "bestcandidate dau michelscore in the vertex"    );
  AddVarI  (output, bestcandidate_dau_vtx_nhits,       "bestcandidate dau points in the vertex"         );

  AddVarI  (output, bestcandidate_parent_ndau,       "bestcandidate parent' daughters"       );
  AddVar4VF(output, bestcandidate_parent_pos,        "bestcandidate parent position"         ); 
  AddVar3VF(output, bestcandidate_parent_dir,        "bestcandidate parent direction"        );
  AddVar4VF(output, bestcandidate_parent_endpos,     "bestcandidate parent position"         ); 
  AddVar3VF(output, bestcandidate_parent_enddir,     "bestcandidate parent direction"        );
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonBestCandidateRecoPID(OutputManager& output){
//********************************************************************

  AddVarF  (output, bestcandidate_chi2_prot_perndf_5,  "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_muon_perndf_5,  "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_kaon_perndf_5,  "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_prot_perndf_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_muon_perndf_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_kaon_perndf_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_prot_perndf_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_muon_perndf_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_kaon_perndf_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_prot_perndf_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_muon_perndf_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_kaon_perndf_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_prot_perndf_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_muon_perndf_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_chi2_kaon_perndf_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_range, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_range, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_muon, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_muon_range, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot_5, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon_5, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon_5, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_prot_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_kaon_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_lkl_muon_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_5, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_5, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_10, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_15, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_20, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_prot_25, "bestcandidate chi2 proton");
  AddVarF  (output, bestcandidate_freelkl_kaon_25, "bestcandidate chi2 proton");
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
  AddVarI (output, bestcandidate_trueId,         "bestcandidate true Id");
  AddVarF (output, bestcandidate_trueendmom_atAPA,     "bestcandidate true end momentum at APA border");

  AddVarI (output, bestcandidate_dau_truendau,    "bestcandidate daughter' true ndaughters"  );
  AddVarI (output, bestcandidate_dau_truepdg,     "bestcandidate daughter true pdg"          );
  AddVar4VF(output, bestcandidate_dau_truepos,    "bestcandidate daughter true position"     );
  AddVar4VF(output, bestcandidate_dau_trueendpos, "bestcandidate daughter true end position" );
  AddVarI (output, bestcandidate_dau_trueproc,    "bestcandidate daughter true process"      );
  AddVarI (output, bestcandidate_dau_trueendproc, "bestcandidate daughter true end process"  );
  AddVarF (output, bestcandidate_dau_truemom,     "bestcandidate daughter true momentum"     );
  AddVarF (output, bestcandidate_dau_trueendmom,  "bestcandidate daughter true end momentum" );

  AddVarI (output, bestcandidate_parent_trueId,   "bestcandidate parent true Id");
  AddVar4VF (output, bestcandidate_parent_truepos,  "bestcandidate parent true pos");
  AddVar4VF (output, bestcandidate_parent_trueendpos,  "bestcandidate parent true end pos");
  AddVar3VF (output, bestcandidate_parent_truedir,  "bestcandidate parent true dir");
  AddVar3VF (output, bestcandidate_parent_trueenddir,  "bestcandidate parent true end dir");
  AddVarF (output, bestcandidate_parent_trueendmom_atAPA,   "bestcandidate parent true end mom at apa");
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
  Float_t calE = (Float_t)pdAnaUtils::ComputeDepositedEnergy(part);
  Float_t mom  = sqrt(pow(calE,2)+2*calE*938.27);
  Float_t csda = pdAnaUtils::ComputeCSDARange(mom, 2212);
  output.FillVectorVar         (candidates_calE,             calE           );
  output.FillVectorVar         (candidates_type,             part->Type             );
  output.FillMatrixVarFromArray(candidates_CNNscore,         part->CNNscore,       3); 
  output.FillVectorVar         (candidates_chi2_prot,        part->Chi2Proton       );
  output.FillVectorVar         (candidates_chi2_kaon,        (Float_t)pdAnaUtils::Chi2PID(*part,321).first);
  output.FillVectorVar         (candidates_chi2_muon,        part->Chi2Muon         );
  output.FillVectorVar         (candidates_chi2_ndf,         part->Chi2ndf          );
  output.FillVectorVar         (candidates_nhits,            part->NHits            );

  output.FillVectorVar         (candidates_average_dedx,      pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );
  output.FillVectorVar         (candidates_truncated_dedx,   (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,part->Hits[2]));
  output.FillVectorVar         (candidates_vtx_michelscore,  part->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_vtx_nhits,        part->vtx_CNN_NHits);
  output.FillVectorVar         (candidates_distance_tcp,        (Float_t)part->Distance_to_closest_particle);

  if(parent)output.FillVectorVar(candidates_distance_mother, pdAnaUtils::ComputeDistanceMotherDaughter(parent,part));

  if(part->Daughters.empty())return;
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
  output.FillVectorVar         (candidates_dau_nhits,            dau->NHits            );
  output.FillVectorVar         (candidates_dau_forced,            (int)part->forced_daughter            );
  output.FillVectorVar         (candidates_dau_forced_matched,            (int)part->forced_daughter_matched            );

  output.FillVectorVar         (candidates_dau_average_dedx,      pdAnaUtils::ComputeAveragedEdxOverResRange(dau,5));
  output.FillVectorVar         (candidates_dau_truncated_dedx,   (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dau->Hits[2]));
  Float_t dau_calE = (Float_t)pdAnaUtils::ComputeDepositedEnergy(dau);
  Float_t dau_mom  = sqrt(pow(calE,2)+2*calE*105.07);
  Float_t dau_csda = pdAnaUtils::ComputeCSDARange(dau_mom, 13);
  output.FillVectorVar         (candidates_dau_calE,            dau_calE           );
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

  output.FillVectorVar         (candidates_gdau_average_dedx,      pdAnaUtils::ComputeAveragedEdxOverResRange(gdau,5));
  output.FillVectorVar         (candidates_gdau_calE,             (Float_t)pdAnaUtils::ComputeDepositedEnergy(gdau));
  output.FillVectorVar         (candidates_gdau_vtx_michelscore,  gdau->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_gdau_vtx_nhits,        gdau->vtx_CNN_NHits);
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonCandidatesRecoPID(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  output.FillVectorVar(candidates_lkl_prot, (Float_t)pdAnaUtils::GetdEdxLikelihood(part,2212));
  output.FillVectorVar(candidates_lkl_kaon, (Float_t)pdAnaUtils::GetdEdxLikelihood(part,321));

  if(part->Daughters.empty())return;
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[0]);
  if (!dau) return;
  output.FillVectorVar(candidates_dau_lkl_prot, (Float_t)pdAnaUtils::GetdEdxLikelihood(dau,2212));
  output.FillVectorVar(candidates_dau_lkl_muon, (Float_t)pdAnaUtils::GetdEdxLikelihood(dau,13));
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
    for (int j = 0; j < (int)std::min((int)NMAXHITSPERPLANE,(int)part->Hits[2].size()); j++){
      output.FillMatrixVar(candidates_hit_x,         (Float_t)part->Hits[2][j].Position.X(),   -1, j);
      output.FillMatrixVar(candidates_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),   -1, j);
      output.FillMatrixVar(candidates_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),   -1, j);
      output.FillMatrixVar(candidates_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,           -1, j);
      output.FillMatrixVar(candidates_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,     -1, j);
      output.FillMatrixVar(candidates_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,  -1, j);
    }
  }

  if(part->Daughters.empty())return;
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
    for(int j = 0; j < (int)std::min((int)NMAXHITSPERPLANE,(int)dau->Hits[2].size()); j++){
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

  if(part->Daughters.empty())return;
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
  output.FillVar               (bestcandidate_chi2_kaon,        (Float_t)pdAnaUtils::Chi2PID(*part,321).first);
  output.FillVar               (bestcandidate_chi2_ndf,         part->Chi2ndf          );
  output.FillVar               (bestcandidate_nhits,            part->NHits            );
  output.FillVar               (bestcandidate_average_dedx,     pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );
  output.FillVar               (bestcandidate_vtx_michelscore,  part->vtx_CNN_michelscore);
  output.FillVar               (bestcandidate_vtx_nhits,        part->vtx_CNN_NHits);
  output.FillVar               (bestcandidate_calE,             (Float_t)pdAnaUtils::ComputeDepositedEnergy(part));
  output.FillVar               (bestcandidate_truncated_dedx,   (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,part->Hits[2]));

  if(parent){
    output.FillVar(bestcandidate_distance_mother,         pdAnaUtils::ComputeDistanceMotherDaughter(parent,part));
    output.FillVar               (bestcandidate_parent_ndau,      (Int_t)parent->Daughters.size() );
    output.FillVectorVarFromArray(bestcandidate_parent_pos,              parent->PositionStart,  4);
    output.FillVectorVarFromArray(bestcandidate_parent_dir,              parent->DirectionStart, 3); 
    output.FillVectorVarFromArray(bestcandidate_parent_endpos,           parent->PositionEnd,    4);
    output.FillVectorVarFromArray(bestcandidate_parent_enddir,           parent->DirectionEnd,   3); 
  }

  if(part->Daughters.empty())return;
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

  output.FillVar               (bestcandidate_dau_calE,             (Float_t)pdAnaUtils::ComputeDepositedEnergy(dau));
  output.FillVar               (bestcandidate_dau_vtx_michelscore,  dau->vtx_CNN_michelscore);
  output.FillVar               (bestcandidate_dau_vtx_nhits,        dau->vtx_CNN_NHits);
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonBestCandidateRecoPID(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  std::pair<double,int> result = pdAnaUtils::Chi2PID_UpToRR(*part,2212,5);
  output.FillVar(bestcandidate_chi2_prot_perndf_5,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,13,5);
  output.FillVar(bestcandidate_chi2_muon_perndf_5,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,321,5);
  output.FillVar(bestcandidate_chi2_kaon_perndf_5,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,2212,10);
  output.FillVar(bestcandidate_chi2_prot_perndf_10,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,13,10);
  output.FillVar(bestcandidate_chi2_muon_perndf_10,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,321,10);
  output.FillVar(bestcandidate_chi2_kaon_perndf_10,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,2212,15);
  output.FillVar(bestcandidate_chi2_prot_perndf_15,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,13,15);
  output.FillVar(bestcandidate_chi2_muon_perndf_15,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,321,15);
  output.FillVar(bestcandidate_chi2_kaon_perndf_15,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,2212,20);
  output.FillVar(bestcandidate_chi2_prot_perndf_20,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,13,20);
  output.FillVar(bestcandidate_chi2_muon_perndf_20,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,321,20);
  output.FillVar(bestcandidate_chi2_kaon_perndf_20,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,2212,25);
  output.FillVar(bestcandidate_chi2_prot_perndf_25,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,13,25);
  output.FillVar(bestcandidate_chi2_muon_perndf_25,(Float_t)result.first/result.second);
  result = pdAnaUtils::Chi2PID_UpToRR(*part,321,25);
  output.FillVar(bestcandidate_chi2_kaon_perndf_25,(Float_t)result.first/result.second);

  output.FillVar(bestcandidate_lkl_prot,pdAnaUtils::GetdEdxLikelihood(part,2212));
  output.FillVar(bestcandidate_lkl_kaon,pdAnaUtils::GetdEdxLikelihood(part,321));
  output.FillVar(bestcandidate_lkl_muon,pdAnaUtils::GetdEdxLikelihood(part,13));
  std::pair<Float_t,Float_t>lkl_result = pdAnaUtils::GetdEdxLikelihoodFreeRange(part,2212);
  output.FillVar(bestcandidate_freelkl_prot,lkl_result.first);
  output.FillVar(bestcandidate_freelkl_prot_range,lkl_result.second);
  lkl_result = pdAnaUtils::GetdEdxLikelihoodFreeRange(part,321);
  output.FillVar(bestcandidate_freelkl_kaon,lkl_result.first);
  output.FillVar(bestcandidate_freelkl_kaon_range,lkl_result.second);
  lkl_result = pdAnaUtils::GetdEdxLikelihoodFreeRange(part,13);
  output.FillVar(bestcandidate_freelkl_muon,lkl_result.first);
  output.FillVar(bestcandidate_freelkl_muon_range,lkl_result.second);
  output.FillVar(bestcandidate_lkl_prot_5,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,2212,5));
  output.FillVar(bestcandidate_lkl_kaon_5,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,321,5));
  output.FillVar(bestcandidate_lkl_muon_5,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,13,5));
  output.FillVar(bestcandidate_lkl_prot_10,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,2212,10));
  output.FillVar(bestcandidate_lkl_kaon_10,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,321,10));
  output.FillVar(bestcandidate_lkl_muon_10,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,13,10));
  output.FillVar(bestcandidate_lkl_prot_15,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,2212,15));
  output.FillVar(bestcandidate_lkl_kaon_15,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,321,15));
  output.FillVar(bestcandidate_lkl_muon_15,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,13,15));
  output.FillVar(bestcandidate_lkl_prot_20,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,2212,20));
  output.FillVar(bestcandidate_lkl_kaon_20,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,321,20));
  output.FillVar(bestcandidate_lkl_muon_20,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,13,20));
  output.FillVar(bestcandidate_lkl_prot_25,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,2212,25));
  output.FillVar(bestcandidate_lkl_kaon_25,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,321,25));
  output.FillVar(bestcandidate_lkl_muon_25,pdAnaUtils::GetdEdxLikelihood_UpToRR(part,13,25));
  output.FillVar(bestcandidate_freelkl_prot_5,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,2212,5).first);
  output.FillVar(bestcandidate_freelkl_kaon_5,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,321,5).first);
  output.FillVar(bestcandidate_freelkl_prot_10,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,2212,10).first);
  output.FillVar(bestcandidate_freelkl_kaon_10,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,321,10).first);
  output.FillVar(bestcandidate_freelkl_prot_15,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,2212,15).first);
  output.FillVar(bestcandidate_freelkl_kaon_15,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,321,15).first);
  output.FillVar(bestcandidate_freelkl_prot_20,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,2212,20).first);
  output.FillVar(bestcandidate_freelkl_kaon_20,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,321,20).first);
  output.FillVar(bestcandidate_freelkl_prot_25,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,2212,25).first);
  output.FillVar(bestcandidate_freelkl_kaon_25,pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(part,321,25).first);
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
    for (int j = 0; j < (int)std::min((int)nmaxhitsperplane,(int)part->Hits[2].size()); j++){
      output.FillVectorVar(bestcandidate_hit_x,         (Float_t)part->Hits[2][j].Position.X(),    j);
      output.FillVectorVar(bestcandidate_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),    j);
      output.FillVectorVar(bestcandidate_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),    j);
      output.FillVectorVar(bestcandidate_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,            j);
      output.FillVectorVar(bestcandidate_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,      j);
      output.FillVectorVar(bestcandidate_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,   j);
    }
  }
  
  if(part->Daughters.empty())return;
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
    for (int j = 0; j < (int)std::min((int)nmaxhitsperplane,(int)dau->Hits[2].size()); j++){
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
void kaonTree::FillKaonVariables_KaonBestCandidateTrue(OutputManager& output, AnaParticlePD* part, AnaParticlePD* parent){
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
  output.FillVar               (bestcandidate_trueendmom_atAPA,   (Float_t)pdAnaUtils::EstimateTrueMomAtAPABorder(part));
  output.FillVar               (bestcandidate_trueId,             truePart->ID     );

  if(parent){
    AnaTrueParticle* parentTruePart = static_cast<AnaTrueParticle*>(parent->TrueObject);
    if(parentTruePart){
      output.FillVar(bestcandidate_parent_trueId,             parentTruePart->ID);
      output.FillVectorVarFromArray(bestcandidate_parent_truepos,            parentTruePart->Position,      4);
      output.FillVectorVarFromArray(bestcandidate_parent_trueendpos,         parentTruePart->PositionEnd,   4);
      output.FillVectorVarFromArray(bestcandidate_parent_truedir,            parentTruePart->Direction,      3);
      output.FillVectorVarFromArray(bestcandidate_parent_trueenddir,         parentTruePart->DirectionEnd,   3);
      output.FillVar(bestcandidate_parent_trueendmom_atAPA,             (Float_t)pdAnaUtils::EstimateTrueMomAtAPABorder(parent));
    }
  }

  if(part->Daughters.empty())return;
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
