#include "neutralKaonTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output){
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
}

// void neutralKaonTree::AddVariablesDynamic(OutputManager& output, const std::string& par1) {
//   floatVars[par1 + "_truemom"] = 0.0f;
//   floatVars[par1 + "_trueendmom"] = 0.0f;
//   // intVars[par1 + "_truepdg"] = 0;
//   // intVars[par1 + "_trueparentpdg"] = 0;
//   // intVars[par1 + "_trueparentid"] = 0;
//   // intVars[par1 + "_trueproc"] = 0;
//   // intVars[par1 + "_trueendproc"] = 0;
//   // intVars[par1 + "_truedecay"] = 0;
//   // intVars[par1 + "_truechainmuon"] = 0;
//   // intVars[par1 + "_truendau"] = 0;
//   // vec4Vars[par1 + "_truepos"] = {0.0f, 0.0f, 0.0f, 0.0f};
//   // vec4Vars[par1 + "_trueendpos"] = {0.0f, 0.0f, 0.0f, 0.0f};
//   // vec3Vars[par1 + "_truedir"] = {0.0f, 0.0f, 0.0f};
//   // vec3Vars[par1 + "_trueenddir"] = {0.0f, 0.0f, 0.0f};
//   // intVars[par1 + "_truegeneration"] = 0;
//   floatVars[par1 + "_trueeff"] = 0.0f;
//   floatVars[par1 + "_truepur"] = 0.0f;
//   // intVars[par1 + "_branch"] = 0;

//   AddVarF(output, floatVars[par1 + "_truemom"], par1 + " true momentum");
//   AddVarF(output, floatVars[par1 + "_trueendmom"], par1 + " true end momentum");
//   // AddVarI(output, intVars[par1 + "_truepdg"], par1 + " true pdg");
//   // AddVarI(output, intVars[par1 + "_trueparentpdg"], par1 + " parent true pdg");
//   // AddVarI(output, intVars[par1 + "_trueparentid"], par1 + " parent true ID");
//   // AddVarI(output, intVars[par1 + "_trueproc"], par1 + " true process");
//   // AddVarI(output, intVars[par1 + "_trueendproc"], par1 + " true end process");
//   // AddVarI(output, intVars[par1 + "_truedecay"], par1 + " true decay");
//   // AddVarI(output, intVars[par1 + "_truechainmuon"], par1 + " chain to muon");
//   // AddVarI(output, intVars[par1 + "_truendau"], par1 + " true ndaughters");
//   // AddVar4VF(output, vec4Vars[par1 + "_truepos"], par1 + " true position");
//   // AddVar4VF(output, vec4Vars[par1 + "_trueendpos"], par1 + " true end position");
//   // AddVar3VF(output, vec3Vars[par1 + "_truedir"], par1 + " true direction");
//   // AddVar3VF(output, vec3Vars[par1 + "_trueenddir"], par1 + " true end direction");
//   // AddVarI(output, intVars[par1 + "_truegeneration"], par1 + " true generation");
//   AddVarF(output, floatVars[par1 + "_trueeff"], par1 + " true efficiency");
//   AddVarF(output, floatVars[par1 + "_truepur"], par1 + " true purity");
//   // AddVarI(output, intVars[par1 + "_branch"], "selection branch associated to this " + par1);
// }


//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart){
    //********************************************************************

      if(!truePart)return;

      output.FillVar               (truekaon_truemom,            truePart->Momentum        );
      output.FillVar               (truekaon_trueendmom,         truePart->MomentumEnd     );
      output.FillVar               (truekaon_truepdg,            truePart->PDG             );
      output.FillVar               (truekaon_trueparentpdg,      truePart->ParentPDG       );
      output.FillVar               (truekaon_trueparentid,       truePart->ParentID        );
      output.FillVar               (truekaon_trueproc,           NewFunction(truePart));
      output.FillVar               (truekaon_trueendproc,        truePart->ProcessEnd      );
      output.FillVar               (truekaon_truegeneration,     truePart->Generation      );
      output.FillVar               (truekaon_truendau,    (Int_t)truePart->Daughters.size());
      output.FillVectorVarFromArray(truekaon_truepos,            truePart->Position,      4);
      output.FillVectorVarFromArray(truekaon_trueendpos,         truePart->PositionEnd,   4);
      output.FillVectorVarFromArray(truekaon_truedir,            truePart->Direction,     3);
      output.FillVectorVarFromArray(truekaon_trueenddir,         truePart->DirectionEnd,  3);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueDaughter1Candidates(OutputManager& output){
  //********************************************************************

    AddVarF  (output, truedau1_truemom,        "true dau1 candidate true momentum"     );
    AddVarF  (output, truedau1_trueendmom,     "true dau1 candidate true end momentum" );
    AddVarI  (output, truedau1_truepdg,        "true dau1 candidate true pdg"          );
    AddVarI  (output, truedau1_trueparentpdg,  "true dau1 candidate parent true pdg"   );
    AddVarI  (output, truedau1_trueparentid,   "true dau1 candidate parent true ID"    );
    AddVarI  (output, truedau1_trueproc,       "true dau1 candidate true process"      );
    AddVarI  (output, truedau1_trueendproc,    "true dau1 candidate true end process"  );
    AddVarI  (output, truedau1_truedecay,      "true dau1 candidate true decay"        );
    AddVarI  (output, truedau1_truechainmuon,  "true dau1 candidate chain to muon"     );
    AddVarI  (output, truedau1_truendau,       "true dau1 candidate true ndaughters"   );
    AddVar4VF(output, truedau1_truepos,        "true dau1 candidate true position"     );
    AddVar4VF(output, truedau1_trueendpos,     "true dau1 candidate true end position" );
    AddVar3VF(output, truedau1_truedir,        "true dau1 candidate true direction"    );
    AddVar3VF(output, truedau1_trueenddir,     "true dau1 candidate true end direction");
    AddVarI  (output, truedau1_truegeneration, "true dau1 candidate true generation"   );
    AddVarF  (output, truedau1_trueeff,        "true dau1 candidate true efficciency"  );
    AddVarF  (output, truedau1_truepur,        "true dau1 candidate true purity"       );
    AddVarI  (output, truedau1_branch,         "selection branch associated to this true dau1");
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueDaughter1Candidates(OutputManager& output, const AnaTrueParticlePD* truePart){
  //********************************************************************

    if(!truePart)return;

    output.FillVar               (truedau1_truemom,            truePart->Momentum        );
    output.FillVar               (truedau1_trueendmom,         truePart->MomentumEnd     );
    output.FillVar               (truedau1_truepdg,            truePart->PDG             );
    output.FillVar               (truedau1_trueparentpdg,      truePart->ParentPDG       );
    output.FillVar               (truedau1_trueparentid,       truePart->ParentID        );
    output.FillVar               (truedau1_trueproc,           NewFunction(truePart));
    output.FillVar               (truedau1_trueendproc,        truePart->ProcessEnd      );
    output.FillVar               (truedau1_truegeneration,     truePart->Generation      );
    output.FillVar               (truedau1_truendau,    (Int_t)truePart->Daughters.size());
    output.FillVectorVarFromArray(truedau1_truepos,            truePart->Position,      4);
    output.FillVectorVarFromArray(truedau1_trueendpos,         truePart->PositionEnd,   4);
    output.FillVectorVarFromArray(truedau1_truedir,            truePart->Direction,     3);
    output.FillVectorVarFromArray(truedau1_trueenddir,         truePart->DirectionEnd,  3);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueDaughter2Candidates(OutputManager& output){
  //********************************************************************

    AddVarF  (output, truedau2_truemom,        "true dau2 candidate true momentum"     );
    AddVarF  (output, truedau2_trueendmom,     "true dau2 candidate true end momentum" );
    AddVarI  (output, truedau2_truepdg,        "true dau2 candidate true pdg"          );
    AddVarI  (output, truedau2_trueparentpdg,  "true dau2 candidate parent true pdg"   );
    AddVarI  (output, truedau2_trueparentid,   "true dau2 candidate parent true ID"    );
    AddVarI  (output, truedau2_trueproc,       "true dau2 candidate true process"      );
    AddVarI  (output, truedau2_trueendproc,    "true dau2 candidate true end process"  );
    AddVarI  (output, truedau2_truedecay,      "true dau2 candidate true decay"        );
    AddVarI  (output, truedau2_truechainmuon,  "true dau2 candidate chain to muon"     );
    AddVarI  (output, truedau2_truendau,       "true dau2 candidate true ndaughters"   );
    AddVar4VF(output, truedau2_truepos,        "true dau2 candidate true position"     );
    AddVar4VF(output, truedau2_trueendpos,     "true dau2 candidate true end position" );
    AddVar3VF(output, truedau2_truedir,        "true dau2 candidate true direction"    );
    AddVar3VF(output, truedau2_trueenddir,     "true dau2 candidate true end direction");
    AddVarI  (output, truedau2_truegeneration, "true dau2 candidate true generation"   );
    AddVarF  (output, truedau2_trueeff,        "true dau2 candidate true efficciency"  );
    AddVarF  (output, truedau2_truepur,        "true dau2 candidate true purity"       );
    AddVarI  (output, truedau2_branch,         "selection branch associated to this true dau2");
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueDaughter2Candidates(OutputManager& output, const AnaTrueParticlePD* truePart){
  //********************************************************************

    if(!truePart)return;

    output.FillVar               (truedau2_truemom,            truePart->Momentum        );
    output.FillVar               (truedau2_trueendmom,         truePart->MomentumEnd     );
    output.FillVar               (truedau2_truepdg,            truePart->PDG             );
    output.FillVar               (truedau2_trueparentpdg,      truePart->ParentPDG       );
    output.FillVar               (truedau2_trueparentid,       truePart->ParentID        );
    output.FillVar               (truedau2_trueproc,           NewFunction(truePart));
    output.FillVar               (truedau2_trueendproc,        truePart->ProcessEnd      );
    output.FillVar               (truedau2_truegeneration,     truePart->Generation      );
    output.FillVar               (truedau2_truendau,    (Int_t)truePart->Daughters.size());
    output.FillVectorVarFromArray(truedau2_truepos,            truePart->Position,      4);
    output.FillVectorVarFromArray(truedau2_trueendpos,         truePart->PositionEnd,   4);
    output.FillVectorVarFromArray(truedau2_truedir,            truePart->Direction,     3);
    output.FillVectorVarFromArray(truedau2_trueenddir,         truePart->DirectionEnd,  3);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueDaughtersDynamic(OutputManager& output, int maxDaughters){
  //********************************************************************

    // Use truekaon_truendau as dynamic length (rows) and maxDaughters as columns
    AddVarMaxSizeVD(output, truedau_dynamic_truemom,        "true dynamic daughters true momentum",      truekaon_truendau, maxDaughters);
    AddVarMaxSizeVD(output, truedau_dynamic_trueendmom,     "true dynamic daughters true end momentum",  truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_truepdg,        "true dynamic daughters true pdg",           truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_trueparentpdg,  "true dynamic daughters parent true pdg",    truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_trueparentid,   "true dynamic daughters parent true ID",     truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_trueproc,       "true dynamic daughters true process",       truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_trueendproc,    "true dynamic daughters true end process",   truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_truedecay,      "true dynamic daughters true decay",         truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_truechainmuon,  "true dynamic daughters chain to muon",      truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_truendau,       "true dynamic daughters true ndaughters",    truekaon_truendau, maxDaughters);
    AddVarMaxSize4MF(output, truedau_dynamic_truepos,       "true dynamic daughters true position",      truekaon_truendau, maxDaughters);
    AddVarMaxSize4MF(output, truedau_dynamic_trueendpos,    "true dynamic daughters true end position",  truekaon_truendau, maxDaughters);
    AddVarMaxSize3MF(output, truedau_dynamic_truedir,       "true dynamic daughters true direction",     truekaon_truendau, maxDaughters);
    AddVarMaxSize3MF(output, truedau_dynamic_trueenddir,    "true dynamic daughters true end direction", truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_truegeneration, "true dynamic daughters true generation",    truekaon_truendau, maxDaughters);
    AddVarMaxSizeVD(output, truedau_dynamic_trueeff,        "true dynamic daughters true efficiency",     truekaon_truendau, maxDaughters);
    AddVarMaxSizeVD(output, truedau_dynamic_truepur,        "true dynamic daughters true purity",         truekaon_truendau, maxDaughters);
    AddVarMaxSizeVI(output, truedau_dynamic_branch,         "true dynamic daughters branch",              truekaon_truendau, maxDaughters);
  }

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueDaughtersDynamic(OutputManager& output, const AnaTrueParticlePD* truePart, int maxDaughters){
  //********************************************************************

    if(!truePart)return;

    // Reset the dynamic counter at the start of this particle
    output.InitializeCounter(truekaon_truendau);

    // Get the daughters of this true particle
    const std::vector<int>& daughters = truePart->Daughters;
    int nDaughters = std::min((int)daughters.size(), maxDaughters);

    // Fill variables for each daughter up to maxDaughters
    for(int i = 0; i < nDaughters; i++) {
        // Get the actual daughter particle from the TrueParticles collection
        // Note: This requires access to the TrueParticles collection, which should be passed as a parameter
        // For now, we'll need to modify the function signature or use a different approach

        // Since we don't have direct access to the TrueParticles collection here,
        // we'll need to modify the calling function to pass the daughter particles directly
        // This is a limitation of the current architecture

        // For now, we'll fill with placeholder values and add a comment explaining the issue
        // TODO: Modify function signature to accept daughter particles directly or pass TrueParticles collection

        // Fill momentum (placeholder - should be daughter's momentum)
        output.FillVectorVar(truedau_dynamic_truemom,        (Double_t)0.0);
        output.FillVectorVar(truedau_dynamic_trueendmom,     (Double_t)0.0);

        // Fill PDG and parent info (placeholder - should be daughter's info)
        output.FillVectorVar(truedau_dynamic_truepdg,        (Int_t)0);
        output.FillVectorVar(truedau_dynamic_trueparentpdg,  (Int_t)truePart->PDG); // Parent is the current particle
        output.FillVectorVar(truedau_dynamic_trueparentid,   (Int_t)truePart->ID);   // Parent ID is the current particle's ID

        // Fill process info (placeholder)
        output.FillVectorVar(truedau_dynamic_trueproc,       (Int_t)0);
        output.FillVectorVar(truedau_dynamic_trueendproc,    (Int_t)0);

        // Fill decay and chain info (placeholder)
        output.FillVectorVar(truedau_dynamic_truedecay,      (Int_t)0);
        output.FillVectorVar(truedau_dynamic_truechainmuon,  (Int_t)0);

        // Fill daughter count and generation (placeholder)
        output.FillVectorVar(truedau_dynamic_truendau,       (Int_t)0);
        output.FillVectorVar(truedau_dynamic_truegeneration, (Int_t)0);

        // Fill position and direction arrays (placeholder)
        float pos[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        float dir[3] = {0.0f, 0.0f, 0.0f};
        output.FillMatrixVarFromArray(truedau_dynamic_truepos,     pos, 4);
        output.FillMatrixVarFromArray(truedau_dynamic_trueendpos,  pos, 4);
        output.FillMatrixVarFromArray(truedau_dynamic_truedir,     dir, 3);
        output.FillMatrixVarFromArray(truedau_dynamic_trueenddir,  dir, 3);

        // Fill efficiency and purity (default values)
        output.FillVectorVar(truedau_dynamic_trueeff,        (Double_t)1.0);
        output.FillVectorVar(truedau_dynamic_truepur,        (Double_t)1.0);

        // Fill branch (default value)
        output.FillVectorVar(truedau_dynamic_branch,         (Int_t)0);
    }
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueDaughtersWithCollection(OutputManager& output, const AnaTrueParticlePD* truePart,
                                                                         const std::vector<AnaTrueParticleB*>& trueParticles, int maxDaughters){
  //********************************************************************

    if(!truePart)return;

    // Get the daughters of this true particle
    const std::vector<int>& daughters = truePart->Daughters;
    int nDaughters = std::min((int)daughters.size(), maxDaughters);

    // Fill variables for each daughter up to maxDaughters
    for(int i = 0; i < nDaughters; i++) {
        // Resolve daughter using helper (matches secondary/standard pattern)
        AnaTrueParticlePD* daughter = pdAnaUtils::GetTrueParticle(trueParticles, daughters[i]);
        if(!daughter) continue; // skip if not found; do not advance counter

        // Fill with actual daughter information
        output.FillVectorVar(truedau_dynamic_truemom,        (Double_t)daughter->Momentum);
        output.FillVectorVar(truedau_dynamic_trueendmom,     (Double_t)daughter->MomentumEnd);

        output.FillVectorVar(truedau_dynamic_truepdg,        (Int_t)daughter->PDG);
        output.FillVectorVar(truedau_dynamic_trueparentpdg,  (Int_t)truePart->PDG); // Parent is the current particle
        output.FillVectorVar(truedau_dynamic_trueparentid,   (Int_t)truePart->ID);  // Parent ID is the current particle's ID

        output.FillVectorVar(truedau_dynamic_trueproc,       (Int_t)NewFunction(daughter));
        output.FillVectorVar(truedau_dynamic_trueendproc,    (Int_t)daughter->ProcessEnd);

        // Fields not available: keep default zeros
        output.FillVectorVar(truedau_dynamic_truedecay,      (Int_t)0);
        output.FillVectorVar(truedau_dynamic_truechainmuon,  (Int_t)0);

        output.FillVectorVar(truedau_dynamic_truendau,       (Int_t)daughter->Daughters.size());
        output.FillVectorVar(truedau_dynamic_truegeneration, (Int_t)daughter->Generation);

        output.FillMatrixVarFromArray(truedau_dynamic_truepos,     daughter->Position,     4);
        output.FillMatrixVarFromArray(truedau_dynamic_trueendpos,  daughter->PositionEnd,  4);
        output.FillMatrixVarFromArray(truedau_dynamic_truedir,     daughter->Direction,    3);
        output.FillMatrixVarFromArray(truedau_dynamic_trueenddir,  daughter->DirectionEnd, 3);

        // Efficiency and purity might need to be calculated or set to default values
        output.FillVectorVar(truedau_dynamic_trueeff,        (Double_t)1.0);
        output.FillVectorVar(truedau_dynamic_truepur,        (Double_t)1.0);

        output.FillVectorVar(truedau_dynamic_branch,         (Int_t)0);
        // advance dynamic counter for filled daughter
        output.IncrementCounter(truekaon_truendau);
    }
    // scalar truekaon_truendau is already filled elsewhere; do not overwrite here
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueParentCandidates(OutputManager& output){
  //********************************************************************

    AddVarF  (output, truepar_truemom,        "true par candidate true momentum"     );
    AddVarF  (output, truepar_trueendmom,     "true par candidate true end momentum" );
    AddVarI  (output, truepar_truepdg,        "true par candidate true pdg"          );
    AddVarI  (output, truepar_trueparentpdg,  "true par candidate parent true pdg"   );
    AddVarI  (output, truepar_trueparentid,   "true par candidate parent true ID"    );
    AddVarI  (output, truepar_trueproc,       "true par candidate true process"      );
    AddVarI  (output, truepar_trueendproc,    "true par candidate true end process"  );
    AddVarI  (output, truepar_truedecay,      "true par candidate true decay"        );
    AddVarI  (output, truepar_truechainmuon,  "true par candidate chain to muon"     );
    AddVarI  (output, truepar_truendau,       "true par candidate true ndaughters"   );
    AddVar4VF(output, truepar_truepos,        "true par candidate true position"     );
    AddVar4VF(output, truepar_trueendpos,     "true par candidate true end position" );
    AddVar3VF(output, truepar_truedir,        "true par candidate true direction"    );
    AddVar3VF(output, truepar_trueenddir,     "true par candidate true end direction");
    AddVarI  (output, truepar_truegeneration, "true par candidate true generation"   );
    AddVarF  (output, truepar_trueeff,        "true par candidate true efficciency"  );
    AddVarF  (output, truepar_truepur,        "true par candidate true purity"       );
    AddVarI  (output, truepar_branch,         "selection branch associated to this true par");
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueParentCandidates(OutputManager& output, const AnaTrueParticlePD* truePart){
  //********************************************************************

    if(!truePart)return;

    output.FillVar               (truepar_truemom,            truePart->Momentum        );
    output.FillVar               (truepar_trueendmom,         truePart->MomentumEnd     );
    output.FillVar               (truepar_truepdg,            truePart->PDG             );
    output.FillVar               (truepar_trueparentpdg,      truePart->ParentPDG       );
    output.FillVar               (truepar_trueparentid,       truePart->ParentID        );
    output.FillVar               (truepar_trueproc,           NewFunction(truePart));
    output.FillVar               (truepar_trueendproc,        truePart->ProcessEnd      );
    output.FillVar               (truepar_truegeneration,     truePart->Generation      );
    output.FillVar               (truepar_truendau,    (Int_t)truePart->Daughters.size());
    output.FillVectorVarFromArray(truepar_truepos,            truePart->Position,      4);
    output.FillVectorVarFromArray(truepar_trueendpos,         truePart->PositionEnd,   4);
    output.FillVectorVarFromArray(truepar_truedir,            truePart->Direction,     3);
    output.FillVectorVarFromArray(truepar_trueenddir,         truePart->DirectionEnd,  3);
}

//********************************************************************
void neutralKaonTree::AddNeutralKaonVariables_TrueGrandParentCandidates(OutputManager& output){
  //********************************************************************

    AddVarF  (output, truegpar_truemom,        "true gpar candidate true momentum"     );
    AddVarF  (output, truegpar_trueendmom,     "true gpar candidate true end momentum" );
    AddVarI  (output, truegpar_truepdg,        "true gpar candidate true pdg"          );
    AddVarI  (output, truegpar_trueparentpdg,  "true gpar candidate gparent true pdg"   );
    AddVarI  (output, truegpar_trueparentid,   "true gpar candidate gparent true ID"    );
    AddVarI  (output, truegpar_trueproc,       "true gpar candidate true process"      );
    AddVarI  (output, truegpar_trueendproc,    "true gpar candidate true end process"  );
    AddVarI  (output, truegpar_truedecay,      "true gpar candidate true decay"        );
    AddVarI  (output, truegpar_truechainmuon,  "true gpar candidate chain to muon"     );
    AddVarI  (output, truegpar_truendau,       "true gpar candidate true ndaughters"   );
    AddVar4VF(output, truegpar_truepos,        "true gpar candidate true position"     );
    AddVar4VF(output, truegpar_trueendpos,     "true gpar candidate true end position" );
    AddVar3VF(output, truegpar_truedir,        "true gpar candidate true direction"    );
    AddVar3VF(output, truegpar_trueenddir,     "true gpar candidate true end direction");
    AddVarI  (output, truegpar_truegeneration, "true gpar candidate true generation"   );
    AddVarF  (output, truegpar_trueeff,        "true gpar candidate true efficciency"  );
    AddVarF  (output, truegpar_truepur,        "true gpar candidate true purity"       );
    AddVarI  (output, truegpar_branch,         "selection branch associated to this true gpar");
}

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueGrandParentCandidates(OutputManager& output, const AnaTrueParticlePD* truePart){
  //********************************************************************

    if(!truePart)return;

    output.FillVar               (truegpar_truemom,            truePart->Momentum        );
    output.FillVar               (truegpar_trueendmom,         truePart->MomentumEnd     );
    output.FillVar               (truegpar_truepdg,            truePart->PDG             );
    output.FillVar               (truegpar_trueparentpdg,      truePart->ParentPDG       );
    output.FillVar               (truegpar_trueparentid,       truePart->ParentID        );
    output.FillVar               (truegpar_trueproc,           NewFunction(truePart));
    output.FillVar               (truegpar_trueendproc,        truePart->ProcessEnd      );
    output.FillVar               (truegpar_truegeneration,     truePart->Generation      );
    output.FillVar               (truegpar_truendau,    (Int_t)truePart->Daughters.size());
    output.FillVectorVarFromArray(truegpar_truepos,            truePart->Position,      4);
    output.FillVectorVarFromArray(truegpar_trueendpos,         truePart->PositionEnd,   4);
    output.FillVectorVarFromArray(truegpar_truedir,            truePart->Direction,     3);
    output.FillVectorVarFromArray(truegpar_trueenddir,         truePart->DirectionEnd,  3);
}

AnaTrueParticleB::ProcessEnum neutralKaonTree::NewFunction(const AnaTrueParticlePD* truePart)
{
  return truePart->ProcessStart;
}
