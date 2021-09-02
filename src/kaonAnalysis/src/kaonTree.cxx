#include "kaonTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void kaonTree::AddKaonVariables_CandidateDaughtersTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  Int_t seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMaxSizeVI(output,  seltrk_dau_truesecondary, "is this daughter a true secondary?", seltrk_ndau, nmax);
}

//********************************************************************
void kaonTree::AddKaonVariables_CandidateGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters){
//********************************************************************

  Int_t seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMI   (output, seltrk_gdau_truepdg,     "gdaughters true pdg",                 seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI   (output, seltrk_gdau_truendau,    "gdaughters true number of daughters", seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_truepos,     "gdaughters true position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_trueendpos,  "gdaughters true position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI   (output, seltrk_gdau_trueproc,    "gdaughters true process",             seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI   (output, seltrk_gdau_trueendproc, "gdaughters true end process",         seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_truemom,     "gdaughters true mom",                 seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_trueendmom,  "gdaughters true end mom",             seltrk_ndau, -nmax, nmaxgdaughters);

}


//********************************************************************
void kaonTree::AddKaonVariables_CandidateGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters){
//********************************************************************
  
  Int_t seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMI   (output, seltrk_gdau_ndau,         "gdaughters number of daughters", seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_pos,          "gdaughters position",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D3MF(output, seltrk_gdau_dir,          "gdaughters direction",           seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D4MF(output, seltrk_gdau_endpos,       "gdaughters end position",        seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D3MF(output, seltrk_gdau_enddir,       "gdaughters end direction",       seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_length,       "gdaughters length",              seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_mom_prot,     "gdaughters mom my range, muon",  seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_mom_muon,     "gdaughters mom my range, muon",  seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI   (output, seltrk_gdau_type,         "gdaughters track or shower",     seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3D3MF(output, seltrk_gdau_CNNscore,     "gdaughters CNN score",           seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_chi2_prot,    "gdaughters chi2 proton",         seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_chi2_muon,    "gdaughters chi2 ndf",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMF   (output, seltrk_gdau_chi2_ndf,     "gdaughters chi2 ndf",            seltrk_ndau, -nmax, nmaxgdaughters);
  AddVarMI   (output, seltrk_gdau_nhits,        "gdaughters nhits",               seltrk_ndau, -nmax, nmaxgdaughters);
  AddVar3DMF (output, seltrk_gdau_hit_dedx,     "gdaughters hit dedx",            seltrk_ndau, -nmax, nmaxgdaughters, NMAXHITSPERPLANE);
  AddVar3DMF (output, seltrk_gdau_hit_resrange, "gdaughters hit residual range",  seltrk_ndau, -nmax, nmaxgdaughters, NMAXHITSPERPLANE);
}  



//********************************************************************
void kaonTree::AddKaonVariables_CandidateGGDaughtersTrue(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters){
//********************************************************************

  Int_t seltrk_ndau = standardPDTree::seltrk_ndau;
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
void kaonTree::AddKaonVariables_CandidateGGDaughtersReco(OutputManager& output, UInt_t nmax, UInt_t nmaxgdaughters, UInt_t nmaxggdaughters){
//********************************************************************

  Int_t seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVar3DMI(output, seltrk_ggdau_ndau,      "ggdaughters number of daughters", seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posX,      "ggdaughters position X",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posY,      "ggdaughters position Y",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_posZ,      "ggdaughters position Z",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirX,      "ggdaughters direction X",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirY,      "ggdaughters direction Y",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_dirZ,      "ggdaughters direction Z",         seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposX,   "ggdaughters end position X",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposY,   "ggdaughters end position Y",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_endposZ,   "ggdaughters end position Z",      seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirX,   "ggdaughters end direction X",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirY,   "ggdaughters end direction Y",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_enddirZ,   "ggdaughters end direction Z",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_length,    "ggdaughters length",              seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_mom_muon,  "gdaughters mom range (muon)",     seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_mom_prot,  "gdaughters mom range (proton)",   seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_type,      "gdaughters object type",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_CNNscore0, "gdaughters CNN score 0",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_CNNscore1, "gdaughters CNN score 1",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_CNNscore2, "gdaughters CNN score 2",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_chi2_muon, "gdaughters chi2 muon",            seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_chi2_prot, "gdaughters chi2 proton",          seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMF(output, seltrk_ggdau_chi2_ndf,  "gdaughters chi2 ndf",             seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
  AddVar3DMI(output, seltrk_ggdau_nhits,     "ggdaughters #hits",               seltrk_ndau, -nmax, nmaxgdaughters, nmaxggdaughters);
}

//********************************************************************
void kaonTree::AddKaonVariables_TrueKaonCandidates(OutputManager& output, UInt_t nmax){
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

  AddVarMaxSizeVI (output, candidates_generation,   "candidates generation",        candidates, nmax);
  AddVarMaxSizeVI (output, candidates_ndau,         "candidates' daughters",        candidates, nmax);
  AddVarMaxSize4MF(output, candidates_pos,          "candidates position",          candidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dir,          "candidates direction",         candidates, nmax);
  AddVarMaxSize4MF(output, candidates_endpos,       "candidates position",          candidates, nmax); 
  AddVarMaxSize3MF(output, candidates_enddir,       "candidates direction",         candidates, nmax);
  AddVarMaxSizeVF (output, candidates_length,       "candidates length",            candidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_muon,     "candidates momentum (muon)",   candidates, nmax);
  AddVarMaxSizeVF (output, candidates_mom_prot,     "candidates momentum (proton)", candidates, nmax);
  AddVarMaxSizeVI (output, candidates_type,         "candidates object type",       candidates, nmax);
  AddVarMaxSize3MF(output, candidates_CNNscore,     "candidates CNN score",         candidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_prot,    "candidates chi2 proton",       candidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_muon,    "candidates chi2 proton",       candidates, nmax);
  AddVarMaxSizeVF (output, candidates_chi2_ndf,     "candidates chi2 ndf",          candidates, nmax);
  AddVarMaxSizeVF (output, candidates_distance_dau, "candidates-daughter distance", candidates, nmax);
  AddVarMaxSizeVF (output, candidates_cos_dau,      "candidates-daughter cos",      candidates, nmax);
  AddVarMaxSizeVI (output, candidates_nhits,        "candidates #hits",             candidates, nmax);
 
  AddVarMaxSizeVF (output, candidates_averagedEdx,     "candidates average dEdx/hit",          candidates, nmax);
  AddVarMaxSizeVF (output, candidates_vtx_michelscore, "candidates michelscore in the vertex", candidates, nmax);
  AddVarMaxSizeVI (output, candidates_vtx_nhits,       "candidates points in the vertex",      candidates, nmax);

  AddVarMaxSizeVI (output, candidates_dau_ndau,       "candidates daughter' daughters",        candidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_pos,        "candidates daughter position",          candidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_dir,        "candidates daughter direction",         candidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_endpos,     "candidates daughter position",          candidates, nmax); 
  AddVarMaxSize3MF(output, candidates_dau_enddir,     "candidates daughter direction",         candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_length,     "candidates daughter length",            candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_muon,   "candidates daughter momentum (muon)",   candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_mom_prot,   "candidates daughter momentum (proton)", candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_type,       "candidates daughter object type",       candidates, nmax);
  AddVarMaxSize3MF(output, candidates_dau_CNNscore,   "candidates daughter CNN score",         candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_prot,  "candidates daughter chi2 proton",       candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_muon,  "candidates daughter chi2 proton",       candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_chi2_ndf,   "candidates daughter chi2 ndf",          candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_nhits,      "candidates daughter #hits",             candidates, nmax);

  AddVarMaxSizeVF (output, candidates_dau_calE,            "candidates dau calorimetry energy deposition", candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_vtx_michelscore, "candidates dau michelscore in the vertex",     candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_vtx_nhits,       "candidates dau points in the vertex",          candidates, nmax);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane){
//********************************************************************

  AddVarMF(output, candidates_hit_x,        "candidates x per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_y,        "candidates y per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_z,        "candidates z per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_dedx,     "candidates dEdx per hit",       candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_dedx_cal, "candidates dEdx cal per hit",   candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_hit_resrange, "candidates hit residual range", candidates, -nmax, nmaxhitsperplane);

  AddVarMF(output, candidates_dau_hit_x,        "candidates daughter x per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_y,        "candidates daughter y per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_z,        "candidates daughter z per hit",          candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_dedx,     "candidates daughter dEdx per hit",       candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_dedx_cal, "candidates daughter dEdx cal per hit",   candidates, -nmax, nmaxhitsperplane);
  AddVarMF(output, candidates_dau_hit_resrange, "candidates daughter hit residual range", candidates, -nmax, nmaxhitsperplane);
}

//********************************************************************
void kaonTree::AddKaonVariables_KaonCandidatesTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, candidates_truendau,       "candidates' true ndaughters",  candidates, nmax);
  AddVarMaxSizeVI (output, candidates_truegeneration, "candidates true generation",   candidates, nmax);
  AddVarMaxSizeVI (output, candidates_truepdg,        "candidates true pdg",          candidates, nmax);
  AddVarMaxSize4MF(output, candidates_truepos,        "daugthers true position",      candidates, nmax);
  AddVarMaxSize4MF(output, candidates_trueendpos,     "daugthers true end position",  candidates, nmax);
  AddVarMaxSizeVI (output, candidates_trueproc,       "daugthers true process",       candidates, nmax);
  AddVarMaxSizeVI (output, candidates_trueendproc,    "candidates true end process",  candidates, nmax);
  AddVarMaxSizeVF (output, candidates_truemom,        "candidates true momentum",     candidates, nmax);
  AddVarMaxSizeVF (output, candidates_trueendmom,     "candidates true end momentum", candidates, nmax);

  AddVarMaxSizeVI (output, candidates_dau_truendau,    "candidates daughter' true ndaughters",  candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_truepdg,     "candidates daughter true pdg",          candidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_truepos,     "candidates daughter true position",     candidates, nmax);
  AddVarMaxSize4MF(output, candidates_dau_trueendpos,  "candidates daughter true end position", candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_trueproc,    "candidates daughter true process",      candidates, nmax);
  AddVarMaxSizeVI (output, candidates_dau_trueendproc, "candidates daughter true end process",  candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_truemom,     "candidates daughter true momentum",     candidates, nmax);
  AddVarMaxSizeVF (output, candidates_dau_trueendmom,  "candidates daughter true end momentum", candidates, nmax);
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
void kaonTree::FillKaonVariables_CandidateDaughterTrue(OutputManager& output, AnaParticlePD* part, AnaParticlePD* dau){
//********************************************************************

  int parentID  = -999;
  int gparentID = -999;
  int dparentID = -999;

  if(part->TrueObject){
    parentID  = static_cast<AnaTrueParticlePD*>(part->TrueObject)->ID;
    gparentID = static_cast<AnaTrueParticlePD*>(part->TrueObject)->ParentID;
  }

  int isSecondary = 0;
  dparentID = -999;
  if(dau->TrueObject)dparentID = static_cast<AnaTrueParticlePD*>(dau->TrueObject)->ParentID;
  if(dparentID == parentID && dparentID!=-999 && gparentID==0)isSecondary = 1;

  output.FillVectorVar(seltrk_dau_truesecondary,isSecondary);
}

//********************************************************************
void kaonTree::FillKaonVariables_CandidateGDaughterReco(OutputManager& output, AnaParticlePD* gdau, Int_t gdau_index){
//********************************************************************

  if(!gdau) return;
  
  output.FillMatrixVar           (seltrk_gdau_ndau,     (Int_t)gdau->Daughters.size(), -1, gdau_index   );
  output.Fill3DMatrixVarFromArray(seltrk_gdau_pos,             gdau->PositionStart,    -1, gdau_index, 4);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_dir,             gdau->DirectionStart,   -1, gdau_index, 3);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_endpos,          gdau->PositionEnd,      -1, gdau_index, 4);
  output.Fill3DMatrixVarFromArray(seltrk_gdau_enddir,          gdau->DirectionEnd,     -1, gdau_index, 3);
  output.FillMatrixVar           (seltrk_gdau_length,          gdau->Length,           -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_type,     (Int_t)gdau->Type,             -1, gdau_index   );
  output.Fill3DMatrixVarFromArray(seltrk_gdau_CNNscore,        gdau->CNNscore,         -1, gdau_index, 3);
  output.FillMatrixVar           (seltrk_gdau_chi2_prot,       gdau->Chi2Proton,       -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_chi2_muon,       gdau->Chi2Muon,         -1, gdau_index   );
  output.FillMatrixVar           (seltrk_gdau_chi2_ndf,        gdau->Chi2ndf,          -1, gdau_index   );
  
  output.FillMatrixVar           (seltrk_gdau_mom_prot,      pdAnaUtils::ComputeRangeMomentum(gdau->Length,2212), -1, gdau_index);
  output.FillMatrixVar           (seltrk_gdau_mom_muon,      pdAnaUtils::ComputeRangeMomentum(gdau->Length,13  ), -1, gdau_index);

  output.FillMatrixVar           (seltrk_gdau_nhits,    (Int_t)gdau->NHits,            -1, gdau_index   );
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
void kaonTree::FillKaonVariables_CandidateGDaughterTrue(OutputManager& output, AnaParticlePD* gdau, Int_t gdau_index){
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
void kaonTree::FillKaonVariables_CandidateGGDaughterReco(OutputManager& output, AnaParticlePD* ggdau, Int_t gdau_index, Int_t ggdau_index){
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
  output.Fill3DMatrixVar(seltrk_ggdau_type,         ggdau->Type,              -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_CNNscore0,    ggdau->CNNscore[0],       -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_CNNscore1,    ggdau->CNNscore[1],       -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_CNNscore2,    ggdau->CNNscore[2],       -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_chi2_prot,    ggdau->Chi2Proton,        -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_chi2_muon,    ggdau->Chi2Muon,          -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_chi2_ndf,     ggdau->Chi2ndf,           -1, gdau_index, ggdau_index);

  output.Fill3DMatrixVar(seltrk_ggdau_mom_prot, pdAnaUtils::ComputeRangeMomentum(ggdau->Length,2212), -1, gdau_index, ggdau_index);
  output.Fill3DMatrixVar(seltrk_ggdau_mom_muon, pdAnaUtils::ComputeRangeMomentum(ggdau->Length,13  ), -1, gdau_index, ggdau_index);

  output.Fill3DMatrixVar(seltrk_ggdau_nhits,        ggdau->NHits,             -1, gdau_index, ggdau_index);

}

//********************************************************************
void kaonTree::FillKaonVariables_CandidateGGDaughterTrue(OutputManager& output, AnaParticlePD* ggdau, Int_t gdau_index, Int_t ggdau_index){
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

//********************************************************************
void kaonTree::FillKaonVariables_TrueKaonCandidates(OutputManager& output, const kaonAnaTrueVertex& kvtx){
//********************************************************************

  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(kvtx.TrueParticlesVect[0]);

  output.FillVar               (truekaon_truemom,            truePart->Momentum        );
  output.FillVar               (truekaon_trueendmom,         truePart->MomentumEnd     );
  output.FillVar               (truekaon_truepdg,            truePart->PDG             );
  output.FillVar               (truekaon_trueparentpdg,      truePart->ParentPDG       );
  output.FillVar               (truekaon_trueparentid,       truePart->ParentID        );
  output.FillVar               (truekaon_trueproc,           truePart->ProcessStart    ); 
  output.FillVar               (truekaon_trueendproc,        truePart->ProcessEnd      );
  output.FillVar               (truekaon_truedecay,          kvtx.DecayMode            );
  output.FillVar               (truekaon_truechainmuon,      kvtx.ChainMuon            );
  output.FillVar               (truekaon_truendau,    (Int_t)truePart->Daughters.size()); 
  output.FillVectorVarFromArray(truekaon_truepos,            truePart->Position,      4);
  output.FillVectorVarFromArray(truekaon_trueendpos,         truePart->PositionEnd,   4);
  output.FillVectorVarFromArray(truekaon_truedir,            truePart->Direction,     3);
  output.FillVectorVarFromArray(truekaon_trueenddir,         truePart->DirectionEnd,  3);
  output.FillVar               (truekaon_branch,      (Int_t)kvtx.Branch               ); 
  
  if(kvtx.TrueParticlesVect.size()>1){
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
  }
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonCandidatesReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  output.FillVectorVar         (candidates_generation,       part->Generation       );
  output.FillVectorVar         (candidates_ndau,      (Int_t)part->Daughters.size() );
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
  output.FillVectorVar         (candidates_chi2_muon,        part->Chi2Muon         );
  output.FillVectorVar         (candidates_chi2_ndf,         part->Chi2ndf          );
  output.FillVectorVar         (candidates_nhits,            part->NHits            );

  output.FillVectorVar         (candidates_averagedEdx,      pdAnaUtils::ComputeAveragedEdxOverResRange(part)   );
  output.FillVectorVar         (candidates_vtx_michelscore,  part->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_vtx_nhits,        part->vtx_CNN_NHits);

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

  output.FillVectorVar         (candidates_dau_calE,             pdAnaUtils::ComputeKineticEnergy(*part));
  output.FillVectorVar         (candidates_dau_vtx_michelscore,  dau->vtx_CNN_michelscore);
  output.FillVectorVar         (candidates_dau_vtx_nhits,        dau->vtx_CNN_NHits);

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
    int jmin = std::max<int>(0,(int)part->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)part->Hits[2].size(); 
    
    for (int j = jmin; j < jmax; j++){
      output.FillMatrixVar(candidates_hit_x,         (Float_t)part->Hits[2][j].Position.X(),   -1, j-jmin);
      output.FillMatrixVar(candidates_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),   -1, j-jmin);
      output.FillMatrixVar(candidates_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),   -1, j-jmin);
      output.FillMatrixVar(candidates_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,           -1, j-jmin);
      output.FillMatrixVar(candidates_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,     -1, j-jmin);
      output.FillMatrixVar(candidates_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,  -1, j-jmin);
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
    int jmin = std::max<int>(0,(int)dau->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)dau->Hits[2].size(); 
    
    for (int j = jmin; j < jmax; j++){
      output.FillMatrixVar(candidates_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),   -1, j-jmin);
      output.FillMatrixVar(candidates_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),   -1, j-jmin);
      output.FillMatrixVar(candidates_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),   -1, j-jmin);
      output.FillMatrixVar(candidates_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,           -1, j-jmin);
      output.FillMatrixVar(candidates_dau_hit_dedx_cal,  (Float_t)dau->Hits[2][j].dEdx_calib,     -1, j-jmin);
      output.FillMatrixVar(candidates_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,  -1, j-jmin);
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
}

//********************************************************************
void kaonTree::FillKaonVariables_KaonBestCandidateReco(OutputManager& output, AnaParticlePD* part){
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
    int jmin = std::max<int>(0,(int)part->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)part->Hits[2].size(); 
    
    for (int j = jmin; j < jmax; j++){
      output.FillVectorVar(bestcandidate_hit_x,         (Float_t)part->Hits[2][j].Position.X(),    j-jmin);
      output.FillVectorVar(bestcandidate_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),    j-jmin);
      output.FillVectorVar(bestcandidate_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),    j-jmin);
      output.FillVectorVar(bestcandidate_hit_dedx,      (Float_t)part->Hits[2][j].dEdx,            j-jmin);
      output.FillVectorVar(bestcandidate_hit_dedx_cal,  (Float_t)part->Hits[2][j].dEdx_calib,      j-jmin);
      output.FillVectorVar(bestcandidate_hit_resrange,  (Float_t)part->Hits[2][j].ResidualRange,   j-jmin);
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
    int jmin = std::max<int>(0,(int)dau->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)dau->Hits[2].size(); 
    
    for (int j = jmin; j < jmax; j++){
      output.FillVectorVar(bestcandidate_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),    j-jmin);
      output.FillVectorVar(bestcandidate_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),    j-jmin);
      output.FillVectorVar(bestcandidate_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),    j-jmin);
      output.FillVectorVar(bestcandidate_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,            j-jmin);
      output.FillVectorVar(bestcandidate_dau_hit_dedx_cal,  (Float_t)dau->Hits[2][j].dEdx_calib,      j-jmin);
      output.FillVectorVar(bestcandidate_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,   j-jmin);
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
