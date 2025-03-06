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

//********************************************************************
void neutralKaonTree::FillNeutralKaonVariables_TrueNeutralKaonCandidates(OutputManager& output, const AnaTrueParticlePD* truePart){
    //********************************************************************

      if(!truePart)return;

      output.FillVar               (truekaon_truemom,            truePart->Momentum        );
      output.FillVar               (truekaon_trueendmom,         truePart->MomentumEnd     );
      output.FillVar               (truekaon_truepdg,            truePart->PDG             );
      output.FillVar               (truekaon_trueparentpdg,      truePart->ParentPDG       );
      output.FillVar               (truekaon_trueparentid,       truePart->ParentID        );
      output.FillVar               (truekaon_trueproc,           truePart->ProcessStart    );
      output.FillVar               (truekaon_trueendproc,        truePart->ProcessEnd      );
      output.FillVar               (truekaon_truegeneration,     truePart->Generation      );
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