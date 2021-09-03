#include "PDSPAnalyzerTreeConverter.hxx"
#include "InputManager.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"



//********************************************************************
PDSPAnalyzerTreeConverter::PDSPAnalyzerTreeConverter():pdBaseConverter("pduneana/beamana"){
//********************************************************************

  //true matching flag. True if byHits matching desired, false if byE matching desired
  _byHits = true;

  _softwareVersion = "v09_24"; //not sure about this but not important
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillEventInfo(AnaEventInfo* info){
//*****************************************************************************

  info->Run    = run;
  info->SubRun = subrun;
  info->Event  = event;
  info->IsMC   = MC;
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillTrueInfo(AnaSpill* spill){
//*****************************************************************************

  // Clear the true particles vector
  spill->TrueParticles.clear();
  
  // Add the information of the MC truth particle that generated the event
  AnaTrueParticlePD* beamTruePart = MakeTrueParticle();
  FillTrueBeamTrueParticleInfo(beamTruePart);
  spill->TrueParticles.push_back(beamTruePart);    

  // Add the true daughters of the true beam particle
  for(int ipart = 0; ipart < (int)true_beam_daughter_ID->size(); ipart++){
    AnaTrueParticlePD* truePart = MakeTrueParticle();
    FillTrueBeamDaughterTrueParticleInfo(ipart, truePart, beamTruePart);
    beamTruePart->Daughters.push_back((*true_beam_daughter_ID)[ipart]);
    spill->TrueParticles.push_back(truePart);    
  }

  // Add the true granddaughters of the true beam particle
  for(int ipart = 0; ipart < (int)true_beam_grand_daughter_ID->size(); ipart++){
    AnaTrueParticlePD* truePart = MakeTrueParticle();
    FillTrueBeamGrandDaughterTrueParticleInfo(ipart, truePart);
    AnaTrueParticle* parent = pdAnaUtils::GetTrueParticle(spill->TrueParticles,(*true_beam_grand_daughter_parID)[ipart]);
    if(parent)parent->Daughters.push_back((*true_beam_grand_daughter_ID)[ipart]);
    spill->TrueParticles.push_back(truePart);    
  }  
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillTrueBeamTrueParticleInfo(AnaTrueParticlePD* truePart){
//*****************************************************************************

  truePart->ID = true_beam_ID;
  truePart->PDG = true_beam_PDG;

  truePart->Momentum    = true_beam_startP;
  truePart->MomentumEnd = true_beam_endP;

  if(true_beam_startP){
    truePart->Direction[0] = true_beam_startPx/true_beam_startP;
    truePart->Direction[1] = true_beam_startPy/true_beam_startP;
    truePart->Direction[2] = true_beam_startPz/true_beam_startP;
  }
  if(true_beam_endP){
    truePart->DirectionEnd[0] = true_beam_endPx/true_beam_endP;
    truePart->DirectionEnd[1] = true_beam_endPy/true_beam_endP;
    truePart->DirectionEnd[2] = true_beam_endPz/true_beam_endP;
  }
  
  truePart->Position[0] = true_beam_startX;
  truePart->Position[1] = true_beam_startY;
  truePart->Position[2] = true_beam_startZ;
  
  truePart->PositionEnd[0] = true_beam_endX;
  truePart->PositionEnd[1] = true_beam_endY;
  truePart->PositionEnd[2] = true_beam_endZ;
    
  truePart->ProcessEnd   = truePart->ConvertProcess(*true_beam_endProcess);  
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillTrueBeamDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart, AnaTrueParticlePD* parent){
//*****************************************************************************

  truePart->ID  = (*true_beam_daughter_ID)[ipart];
  truePart->PDG = (*true_beam_daughter_PDG)[ipart];
  
  truePart->ParentID   = parent->ID;
  truePart->ParentPDG  = parent->PDG;
  truePart->GParentPDG = parent->ParentPDG;

  truePart->ProcessStart = truePart->ConvertProcess((*true_beam_daughter_Process)[ipart]);
  truePart->ProcessEnd   = truePart->ConvertProcess((*true_beam_daughter_endProcess)[ipart]);      
  
  truePart->Momentum    = (*true_beam_daughter_startP)[ipart];

  if ((*true_beam_daughter_startP)[ipart]){
    truePart->Direction[0] = (*true_beam_daughter_startPx)[ipart]/(*true_beam_daughter_startP)[ipart];
    truePart->Direction[1] = (*true_beam_daughter_startPy)[ipart]/(*true_beam_daughter_startP)[ipart];
    truePart->Direction[2] = (*true_beam_daughter_startPz)[ipart]/(*true_beam_daughter_startP)[ipart];
  }

  truePart->Position[0] = (*true_beam_daughter_startX)[ipart];
  truePart->Position[1] = (*true_beam_daughter_startY)[ipart];
  truePart->Position[2] = (*true_beam_daughter_startZ)[ipart];
  
  truePart->PositionEnd[0] = (*true_beam_daughter_endX)[ipart];
  truePart->PositionEnd[1] = (*true_beam_daughter_endY)[ipart];
  truePart->PositionEnd[2] = (*true_beam_daughter_endZ)[ipart];

  truePart->Length = (*true_beam_daughter_len)[ipart];  
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillTrueBeamGrandDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart){
//*****************************************************************************

  truePart->ID = (*true_beam_grand_daughter_ID)[ipart];
  truePart->PDG = (*true_beam_grand_daughter_PDG)[ipart];

  truePart->ParentID = (*true_beam_grand_daughter_parID)[ipart];
  
  truePart->ProcessStart = truePart->ConvertProcess((*true_beam_grand_daughter_Process)[ipart]);
  truePart->ProcessEnd   = truePart->ConvertProcess((*true_beam_grand_daughter_endProcess)[ipart]);      
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillBeamInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBeamPD* beam){
//*****************************************************************************
  
  beam->GoodSpill = true;
  
  //Fill beam info
  beam->nMomenta = beam_inst_nMomenta;
  beam->nFibers[0] = beam_inst_nFibersP1;
  beam->nFibers[1] = beam_inst_nFibersP2;
  beam->nFibers[2] = beam_inst_nFibersP3;    
  beam->TOF = (*beam_inst_TOF)[0];
  
  // Create the BeamParticle object
  beam->BeamParticle = MakeParticle();
  
  //Fill beam particle info
  beam->BeamParticle->PositionEnd[0] = beam_inst_X;
  beam->BeamParticle->PositionEnd[1] = beam_inst_Y;
  beam->BeamParticle->PositionEnd[2] = beam_inst_Z;
  beam->BeamParticle->DirectionEnd[0] = beam_inst_dirX;
  beam->BeamParticle->DirectionEnd[1] = beam_inst_dirY;
  beam->BeamParticle->DirectionEnd[2] = beam_inst_dirZ;
  
  beam->BeamParticle->Momentum = beam_inst_P;
  
  for(int i = 0; i < (int)beam_inst_PDG_candidates->size(); i++)
    beam->PDGs.push_back(beam_inst_PDG_candidates->at(i));

  // ------------- True info for MC ----------------
  // TODO
  if(_isMC){
    beam->BeamParticle->TrueObject = pdAnaUtils::GetTrueParticle(trueParticles, true_beam_ID);
  }
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillBunchInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBunch* bunch, AnaBeamPD* beam){
//*****************************************************************************

  bunch->Bunch  = 1;
  bunch->Weight = 1;
  bunch->Particles.clear();
  bunch->Vertices.clear();
  
  // The beam particle
  AnaParticlePD* part = MakeParticle();
  FillBeamParticleInfo(trueParticles, part, beam);  
  bunch->Particles.push_back(part);

  // IMPORTANT !!!!! All daughters are reconstructed as both tracks and showers and added to the particle vector.
  //                So the actual number of particles is half of the size of the vector

  // The daughter tracks
  for(size_t i = 0; i < reco_daughter_allTrack_ID->size(); i++){
    AnaParticlePD* dautrk = MakeParticle();
    FillDaughterPFPInfo(i, dautrk);
    FillDaughterParticleTrackInfo(trueParticles, i, dautrk);
    bunch->Particles.push_back(dautrk);
    part->Daughters.push_back(dautrk);
  }

  // The daughter showers
  for(size_t i = 0; i < reco_daughter_allShower_ID->size(); i++){
    AnaParticlePD* dausho = MakeParticle();
    FillDaughterPFPInfo(i, dausho);
    FillDaughterParticleShowerInfo(trueParticles, i, dausho);
    bunch->Particles.push_back(dausho);
    part->Daughters.push_back(dausho);
  }
}


//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillBeamParticleInfo(std::vector<AnaTrueParticleB*>& trueParticles,
                                                     AnaParticlePD* part, AnaBeamPD* beam){
//*****************************************************************************

  // This is the beam particle according to pandora
  part->isPandora = true;
  part->isBeamPart = reco_beam_passes_beam_cuts;

  part->UniqueID  = reco_beam_PFP_ID;

  if      (reco_beam_type == 11) part->Type = AnaParticlePD::kShower;
  else if (reco_beam_type == 13) part->Type = AnaParticlePD::kTrack;
  else                           part->Type = AnaParticlePD::kUnknown;
  
  part->PositionEnd[0]  = reco_beam_endX;
  part->PositionEnd[1]  = reco_beam_endY;
  part->PositionEnd[2]  = reco_beam_endZ;
  
  part->PositionStart[0]  = reco_beam_startX;
  part->PositionStart[1]  = reco_beam_startY;
  part->PositionStart[2]  = reco_beam_startZ;
  
  part->DirectionEnd[0] = reco_beam_trackDirX;
  part->DirectionEnd[1] = reco_beam_trackDirY;
  part->DirectionEnd[2] = reco_beam_trackDirZ;
  
  part->DirectionStart[0] = reco_beam_trackDirX; 
  part->DirectionStart[1] = reco_beam_trackDirY;
  part->DirectionStart[2] = reco_beam_trackDirZ;
    
  part->Length     = reco_beam_len;
  part->Length_alt = reco_beam_alt_len;

  part->RangeMomentum[0] = reco_beam_momByRange_proton;
  part->RangeMomentum[1] = reco_beam_momByRange_muon;

  part->RangeMomentum_alt[0] = reco_beam_momByRange_alt_proton;
  part->RangeMomentum_alt[1] = reco_beam_momByRange_alt_muon;

  part->vtx_CNN_michelscore = reco_beam_vertex_michel_score;
  part->vtx_CNN_NHits = reco_beam_vertex_nHits;

  //hit information
  //TVector3 point;
  //point.SetXYZ(0,0,0);//no spatial information on the input tree

  for(size_t plane = 2; plane < 3; plane++){   // only the collection plane 
    for(size_t j = 0; j < reco_beam_dEdX_SCE->size(); j++){
      // Add hits
      AnaHitPD hit;//(plane,0,0,0, point);

      hit.dEdx          = (*reco_beam_dEdX_SCE)[j];
      hit.dQdx          = (*reco_beam_dQdX_SCE)[j];
      hit.dEdx_calib    = (*reco_beam_calibrated_dEdX_SCE)[j];      
      hit.ResidualRange = (*reco_beam_resRange_SCE)[j];
           
      part->Hits[plane].push_back(hit);
    }
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,(*reco_beam_calibrated_dEdX_SCE));
  }

  // --------- reco_beam_PFP ------------------------

  part->CNNscore[0] = reco_beam_PFP_trackScore_collection;
  part->CNNscore[1] = reco_beam_PFP_emScore_collection;
  part->CNNscore[2] = reco_beam_PFP_michelScore_collection;

  part->NHits = reco_beam_PFP_nHits;
  
  part->Chi2Proton = reco_beam_Chi2_proton;
  part->Chi2ndf    = reco_beam_Chi2_ndof;

  //------ Truth association ------- //TODO

  if (_isMC){
    AnaTrueParticlePD* truePart = MakeTrueParticle();
    FillBeamTrueParticleInfo(truePart);
    trueParticles.push_back(truePart);    
    part->TrueObject = truePart;
    part->TruePur = reco_beam_true_byHits_purity;
  }
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillBeamTrueParticleInfo(AnaTrueParticlePD* truePart){
//*****************************************************************************

  //TODO a better way of chosing is needed
  if (_byHits){
    truePart->ID  = reco_beam_true_byHits_ID;
    truePart->PDG = reco_beam_true_byHits_PDG;

    truePart->Origin  = reco_beam_true_byHits_origin;
    truePart->Matched = reco_beam_true_byHits_matched;

    truePart->ProcessStart = truePart->ConvertProcess(*reco_beam_true_byHits_process);
    truePart->ProcessEnd   = truePart->ConvertProcess(*reco_beam_true_byHits_endProcess);      
    
    truePart->Momentum    = reco_beam_true_byHits_startP;
    truePart->MomentumEnd = reco_beam_true_byHits_endP;
    
    if (reco_beam_true_byHits_startP){
      truePart->Direction[0] = reco_beam_true_byHits_startPx/reco_beam_true_byHits_startP;
      truePart->Direction[1] = reco_beam_true_byHits_startPy/reco_beam_true_byHits_startP;
      truePart->Direction[2] = reco_beam_true_byHits_startPz/reco_beam_true_byHits_startP;
    }
    if (reco_beam_true_byHits_endP){
      truePart->DirectionEnd[0] = reco_beam_true_byHits_endPx/reco_beam_true_byHits_endP;
      truePart->DirectionEnd[1] = reco_beam_true_byHits_endPy/reco_beam_true_byHits_endP;
      truePart->DirectionEnd[2] = reco_beam_true_byHits_endPz/reco_beam_true_byHits_endP;
    }
    
  }
  else{
    truePart->ID  = reco_beam_true_byE_ID;
    truePart->PDG = reco_beam_true_byE_PDG;
    
    truePart->ProcessStart = truePart->ConvertProcess(*reco_beam_true_byE_process);
    truePart->ProcessEnd   = truePart->ConvertProcess(*reco_beam_true_byE_endProcess);      
    
    truePart->Momentum    = reco_beam_true_byE_startP;
    truePart->MomentumEnd = reco_beam_true_byE_endP;
    
    if (reco_beam_true_byE_startP){
      truePart->Direction[0] = reco_beam_true_byE_startPx/reco_beam_true_byE_startP;
      truePart->Direction[1] = reco_beam_true_byE_startPy/reco_beam_true_byE_startP;
      truePart->Direction[2] = reco_beam_true_byE_startPz/reco_beam_true_byE_startP;
    }
    if (reco_beam_true_byE_endP){
      truePart->DirectionEnd[0] = reco_beam_true_byE_endPx/reco_beam_true_byE_endP;
      truePart->DirectionEnd[1] = reco_beam_true_byE_endPy/reco_beam_true_byE_endP;
      truePart->DirectionEnd[2] = reco_beam_true_byE_endPz/reco_beam_true_byE_endP;
    }
  }
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillDaughterPFPInfo(Int_t itrk, AnaParticlePD* part){
//*****************************************************************************

  part->UniqueID  = (*reco_daughter_PFP_ID)[itrk];
  
  part->CNNscore[0] = (*reco_daughter_PFP_trackScore_collection)[itrk];
  part->CNNscore[1] = (*reco_daughter_PFP_emScore_collection)[itrk];
  part->CNNscore[2] = (*reco_daughter_PFP_michelScore_collection)[itrk];

  part->NHits = (*reco_daughter_PFP_nHits)[itrk];

}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillDaughterParticleTrackInfo(std::vector<AnaTrueParticleB*>& trueParticles,
							      Int_t itrk, AnaParticlePD* part){
//*****************************************************************************

  //fill forced track info
  part->TrackID  = (*reco_daughter_allTrack_ID)[itrk];
  
  part->PositionEnd[0]  = (*reco_daughter_allTrack_endX)[itrk];
  part->PositionEnd[1]  = (*reco_daughter_allTrack_endY)[itrk];
  part->PositionEnd[2]  = (*reco_daughter_allTrack_endZ)[itrk];

  part->PositionStart[0]  = (*reco_daughter_allTrack_startX)[itrk];
  part->PositionStart[1]  = (*reco_daughter_allTrack_startY)[itrk];
  part->PositionStart[2]  = (*reco_daughter_allTrack_startZ)[itrk];

  /* Not available in input tree
  part->DirectionEnd[0] = (*reco_daughter_allTrack_trackDirX)[itrk];
  part->DirectionEnd[1] = (*reco_daughter_allTrack_trackDirY)[itrk];
  part->DirectionEnd[2] = (*reco_daughter_allTrack_trackDirZ)[itrk];

  part->DirectionStart[0] = (*reco_daughter_allTrack_trackDirX)[itrk]; 
  part->DirectionStart[1] = (*reco_daughter_allTrack_trackDirY)[itrk];
  part->DirectionStart[2] = (*reco_daughter_allTrack_trackDirZ)[itrk];
  */

  part->Length     = (*reco_daughter_allTrack_len)[itrk];
  part->Length_alt = (*reco_daughter_allTrack_alt_len)[itrk];

  part->RangeMomentum[0] = (*reco_daughter_allTrack_momByRange_proton)[itrk];
  part->RangeMomentum[1] = (*reco_daughter_allTrack_momByRange_muon)[itrk];

  part->RangeMomentum_alt[0] = (*reco_daughter_allTrack_momByRange_alt_proton)[itrk];
  part->RangeMomentum_alt[1] = (*reco_daughter_allTrack_momByRange_alt_muon)[itrk];

  part->Chi2Proton = (*reco_daughter_allTrack_Chi2_proton)[itrk];
  part->Chi2ndf    = (*reco_daughter_allTrack_Chi2_ndof)[itrk];

  part->vtx_CNN_michelscore = (*reco_daughter_allTrack_vertex_michel_score)[itrk];
  part->vtx_CNN_NHits       = (*reco_daughter_allTrack_vertex_nHits)[itrk];

  //fill hits info
  Float_t dedx, dqdx, dedx_cal, resRange;
  part->AveragedEdx = 0;
  part->AveragedQdx = 0;
  Int_t ncontrib=0;
  for(size_t plane = 2; plane < 3; plane++){   // only the last slice 
    UInt_t nHits = std::min((int)NMAXHITSPERPLANE,(int)(*reco_daughter_allTrack_calibrated_dEdX_SCE)[itrk].size());
    for(size_t j = 0; j < nHits; j++){   
      dedx = (*reco_daughter_allTrack_dEdX_SCE)[itrk][j];
      dqdx = (*reco_daughter_allTrack_dQdX_SCE)[itrk][j];
      dedx_cal = (*reco_daughter_allTrack_calibrated_dEdX_SCE)[itrk][j];
      resRange = (*reco_daughter_allTrack_resRange_SCE)[itrk][j];  

      part->AveragedEdx += dedx;
      part->AveragedQdx += dqdx;
      
      // Add hits
      //no spatial information on the input tree
      AnaHitPD hit;//(plane,0,0,0,point);
      
      hit.dEdx         = dedx;
      hit.dQdx         = dqdx;
      hit.dEdx_calib   = dedx_cal;
      hit.ResidualRange= resRange;
      
      part->Hits[plane].push_back(hit);
      
      ncontrib++;
    }
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,(*reco_daughter_allTrack_calibrated_dEdX_SCE)[itrk]);
  }

  if (ncontrib!=0){
    part->AveragedEdx /= 1.*ncontrib;
    part->AveragedQdx /= 1.*ncontrib;
  }

  //------ Truth association ------- //TODO
  if(_isMC){
    // Search for the true-reco association within the vector of TrueParticles
    part->TrueObject = pdAnaUtils::GetTrueParticle(trueParticles, (*reco_daughter_PFP_true_byHits_ID)[itrk]);

    // If not found create a new TrueParticle, fill it, and add it to the vector of TrueParticles
    if(!part->TrueObject){
      AnaTrueParticlePD* truePart = MakeTrueParticle();
      FillDaughterTrueParticleInfo(itrk, truePart);
      trueParticles.push_back(truePart);    
      
      part->TrueObject = truePart;
      
      part->TruePur = (*reco_daughter_PFP_true_byHits_purity)[itrk];
    }
  }

}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillDaughterParticleShowerInfo(std::vector<AnaTrueParticleB*>& trueParticles,
                                                               Int_t itrk, AnaParticlePD* part){
//*****************************************************************************

  part->ShowerID  = (*reco_daughter_allShower_ID)[itrk];

  part->PositionStart[0]  = (*reco_daughter_allShower_startX)[itrk];
  part->PositionStart[1]  = (*reco_daughter_allShower_startY)[itrk];
  part->PositionStart[2]  = (*reco_daughter_allShower_startZ)[itrk];

  part->DirectionStart[0]  = (*reco_daughter_allShower_dirX)[itrk];
  part->DirectionStart[1]  = (*reco_daughter_allShower_dirY)[itrk];
  part->DirectionStart[2]  = (*reco_daughter_allShower_dirZ)[itrk];
  
  part->Length = (*reco_daughter_allShower_len)[itrk];
  
  if (_isMC){
    // Search for the true-reco association within the vector of TrueParticles
    part->TrueObject = pdAnaUtils::GetTrueParticle(trueParticles, (*reco_daughter_PFP_true_byHits_ID)[itrk]);
    
    // If not found create a new TrueParticle, fill it, and add it to the vector of TrueParticles
    // It never enters here because the true particle corresponding to reco_daughter_PFP_true_byHits_ID)[itrk]
    // has been previously created in filldaughterparticletrackinfo
    if (!part->TrueObject){
      AnaTrueParticlePD* truePart = MakeTrueParticle();
      FillDaughterTrueParticleInfo(itrk, truePart);
      trueParticles.push_back(truePart);    
        
      part->TrueObject = truePart;

      part->TruePur = (*reco_daughter_PFP_true_byHits_purity)[itrk];
    }
  }  
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::FillDaughterTrueParticleInfo(Int_t ipart, AnaTrueParticlePD* truePart){
//*****************************************************************************
  
  truePart->ID  = (*reco_daughter_PFP_true_byHits_ID)[ipart];
  truePart->PDG = (*reco_daughter_PFP_true_byHits_PDG)[ipart];
  
  truePart->Origin  = (*reco_daughter_PFP_true_byHits_origin)[ipart];
  
  truePart->ParentID  = (*reco_daughter_PFP_true_byHits_parID)[ipart];
  truePart->ParentPDG = (*reco_daughter_PFP_true_byHits_parPDG)[ipart];
  
  truePart->ProcessStart = truePart->ConvertProcess((*reco_daughter_PFP_true_byHits_process)[ipart]);
  truePart->ProcessEnd   = truePart->ConvertProcess((*reco_daughter_PFP_true_byHits_endProcess)[ipart]);      
  
  truePart->Momentum    = (*reco_daughter_PFP_true_byHits_startP)[ipart];

  if((*reco_daughter_PFP_true_byHits_startP)[ipart]){
    truePart->Direction[0] = (*reco_daughter_PFP_true_byHits_startPx)[ipart]/(*reco_daughter_PFP_true_byHits_startP)[ipart];;
    truePart->Direction[1] = (*reco_daughter_PFP_true_byHits_startPy)[ipart]/(*reco_daughter_PFP_true_byHits_startP)[ipart];;
    truePart->Direction[2] = (*reco_daughter_PFP_true_byHits_startPz)[ipart]/(*reco_daughter_PFP_true_byHits_startP)[ipart];;
  }

  /* not available
    if ((*reco_daughter_PFP_true_byHits_endP)[ipart]){
    truePart->DirectionEnd[0] = (*reco_daughter_PFP_true_byHits_endPx)[ipart]/(*reco_daughter_PFP_true_byHits_endP)[ipart];;
    truePart->DirectionEnd[1] = (*reco_daughter_PFP_true_byHits_endPy)[ipart]/(*reco_daughter_PFP_true_byHits_endP)[ipart];;
    truePart->DirectionEnd[2] = (*reco_daughter_PFP_true_byHits_endPz)[ipart]/(*reco_daughter_PFP_true_byHits_endP)[ipart];;
    }
  */

  truePart->Position[0] = (*reco_daughter_PFP_true_byHits_startX)[ipart];
  truePart->Position[1] = (*reco_daughter_PFP_true_byHits_startY)[ipart];
  truePart->Position[2] = (*reco_daughter_PFP_true_byHits_startZ)[ipart];
  
  truePart->PositionEnd[0] = (*reco_daughter_PFP_true_byHits_endX)[ipart];
  truePart->PositionEnd[1] = (*reco_daughter_PFP_true_byHits_endY)[ipart];
  truePart->PositionEnd[2] = (*reco_daughter_PFP_true_byHits_endZ)[ipart];
  
  truePart->Length =   (*reco_daughter_PFP_true_byHits_len)[ipart];
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::InitializeVariables(){
//*****************************************************************************

  run = 0;
  subrun = 0;
  event = 0;
  MC = 0;
  true_beam_PDG = 0;
  true_beam_mass = 0;
  true_beam_ID = 0;
  true_beam_endProcess = 0;
  true_beam_startX = 0;
  true_beam_startY = 0;
  true_beam_startZ = 0;
  true_beam_endX = 0;
  true_beam_endY = 0;
  true_beam_endZ = 0;
  true_beam_startP = 0;
  true_beam_startPx = 0;
  true_beam_startPy = 0;
  true_beam_startPz = 0;
  true_beam_endP = 0;
  true_beam_endPx = 0;
  true_beam_endPy = 0;
  true_beam_endPz = 0;
  true_beam_last_len = 0;
  true_beam_startDirX = 0;
  true_beam_startDirY = 0;
  true_beam_startDirZ = 0;
  true_beam_nElasticScatters = 0;
  true_beam_elastic_costheta = 0;
  true_beam_elastic_X = 0;
  true_beam_elastic_Y = 0;
  true_beam_elastic_Z = 0;
  true_beam_elastic_deltaE = 0;
  true_beam_elastic_IDE_edep = 0;
  true_beam_IDE_totalDep = 0;
  true_daughter_nPi0 = 0;
  true_daughter_nPiPlus = 0;
  true_daughter_nProton = 0;
  true_daughter_nNeutron = 0;
  true_daughter_nPiMinus = 0;
  true_daughter_nNucleus = 0;
  true_beam_processes = 0;
  true_beam_incidentEnergies = 0;
  true_beam_interactingEnergy = 0;
  true_beam_slices = 0;
  true_beam_slices_deltaE = 0;
  true_beam_traj_X = 0;
  true_beam_traj_Y = 0;
  true_beam_traj_Z = 0;
  true_beam_traj_KE = 0;
  reco_beam_type = 0;
  reco_beam_trackID = 0;
  reco_beam_startX = 0;
  reco_beam_startY = 0;
  reco_beam_startZ = 0;
  reco_beam_endX = 0;
  reco_beam_endY = 0;
  reco_beam_endZ = 0;
  reco_beam_trackDirX = 0;
  reco_beam_trackDirY = 0;
  reco_beam_trackDirZ = 0;
  reco_beam_trackEndDirX = 0;
  reco_beam_trackEndDirY = 0;
  reco_beam_trackEndDirZ = 0;
  reco_beam_flipped = 0;
  reco_beam_len = 0;
  reco_beam_alt_len = 0;
  reco_beam_momByRange_proton = 0;
  reco_beam_momByRange_muon = 0;
  reco_beam_momByRange_alt_proton = 0;
  reco_beam_momByRange_alt_muon = 0;
  reco_beam_calo_startX = 0;
  reco_beam_calo_startY = 0;
  reco_beam_calo_startZ = 0;
  reco_beam_calo_endX = 0;
  reco_beam_calo_endY = 0;
  reco_beam_calo_endZ = 0;
  reco_beam_calo_startDirX = 0;
  reco_beam_calo_startDirY = 0;
  reco_beam_calo_startDirZ = 0;
  reco_beam_calo_endDirX = 0;
  reco_beam_calo_endDirY = 0;
  reco_beam_calo_endDirZ = 0;
  reco_beam_vertex_michel_score = 0;
  reco_beam_vertex_nHits = 0;
  reco_beam_dQdX_SCE = 0;
  reco_beam_dQ = 0;
  reco_beam_dEdX_SCE = 0;
  reco_beam_calibrated_dEdX_SCE = 0;
  reco_beam_resRange_SCE = 0;
  reco_beam_TrkPitch_SCE = 0;
  reco_beam_dQdX_NoSCE = 0;
  reco_beam_dEdX_NoSCE = 0;
  reco_beam_calibrated_dEdX_NoSCE = 0;
  reco_beam_resRange_NoSCE = 0;
  reco_beam_TrkPitch_NoSCE = 0;
  reco_beam_calo_wire = 0;
  reco_beam_calo_wire_z = 0;
  reco_beam_calo_tick = 0;
  reco_beam_calo_TPC = 0;
  reco_beam_passes_beam_cuts = 0;
  reco_beam_Chi2_proton = 0;
  reco_beam_Chi2_ndof = 0;
  reco_beam_PFP_ID = 0;
  reco_beam_PFP_nHits = 0;
  reco_beam_PFP_trackScore = 0;
  reco_beam_PFP_emScore = 0;
  reco_beam_PFP_michelScore = 0;
  reco_beam_PFP_trackScore_collection = 0;
  reco_beam_PFP_emScore_collection = 0;
  reco_beam_PFP_michelScore_collection = 0;
  reco_beam_allTrack_ID = 0;
  reco_beam_allTrack_beam_cuts = 0;
  reco_beam_allTrack_flipped = 0;
  reco_beam_allTrack_len = 0;
  reco_beam_allTrack_startX = 0;
  reco_beam_allTrack_startY = 0;
  reco_beam_allTrack_startZ = 0;
  reco_beam_allTrack_endX = 0;
  reco_beam_allTrack_endY = 0;
  reco_beam_allTrack_endZ = 0;
  reco_beam_allTrack_trackDirX = 0;
  reco_beam_allTrack_trackDirY = 0;
  reco_beam_allTrack_trackDirZ = 0;
  reco_beam_allTrack_trackEndDirX = 0;
  reco_beam_allTrack_trackEndDirY = 0;
  reco_beam_allTrack_trackEndDirZ = 0;
  reco_beam_allTrack_resRange = 0;
  reco_beam_allTrack_calibrated_dEdX = 0;
  reco_beam_allTrack_Chi2_proton = 0;
  reco_beam_allTrack_Chi2_ndof = 0;
  reco_beam_true_byHits_endProcess = 0;
  reco_beam_true_byHits_process = 0;
  reco_beam_true_byHits_origin = 0;
  reco_beam_true_byHits_PDG = 0;
  reco_beam_true_byHits_ID = 0;
  reco_beam_true_byHits_matched = 0;
  reco_beam_true_byHits_purity = 0;
  reco_beam_true_byHits_startE = 0;
  reco_beam_true_byHits_startP = 0;
  reco_beam_true_byHits_startPx = 0;
  reco_beam_true_byHits_startPy = 0;
  reco_beam_true_byHits_startPz = 0;
  reco_beam_true_byHits_endE = 0;
  reco_beam_true_byHits_endP = 0;
  reco_beam_true_byHits_endPx = 0;
  reco_beam_true_byHits_endPy = 0;
  reco_beam_true_byHits_endPz = 0;
  reco_beam_true_byE_startE = 0;
  reco_beam_true_byE_startP = 0;
  reco_beam_true_byE_startPx = 0;
  reco_beam_true_byE_startPy = 0;
  reco_beam_true_byE_startPz = 0;
  reco_beam_true_byE_endE = 0;
  reco_beam_true_byE_endP = 0;
  reco_beam_true_byE_endPx = 0;
  reco_beam_true_byE_endPy = 0;
  reco_beam_true_byE_endPz = 0;
  reco_beam_true_byE_endProcess = 0;
  reco_beam_true_byE_process = 0;
  reco_beam_true_byE_origin = 0;
  reco_beam_true_byE_PDG = 0;
  reco_beam_true_byE_ID = 0;
  reco_beam_true_byE_matched = 0;
  reco_beam_incidentEnergies = 0;
  reco_beam_interactingEnergy = 0;
  reco_track_startX = 0;
  reco_track_startY = 0;
  reco_track_startZ = 0;
  reco_track_endX = 0;
  reco_track_endY = 0;
  reco_track_endZ = 0;
  reco_track_michel_score = 0;
  reco_track_nHits = 0;
  reco_track_ID = 0;
  true_beam_daughter_ID = 0;
  true_beam_daughter_PDG = 0;
  true_beam_daughter_len = 0;
  true_beam_daughter_startX = 0;
  true_beam_daughter_startY = 0;
  true_beam_daughter_startZ = 0;
  true_beam_daughter_endX = 0;
  true_beam_daughter_endY = 0;
  true_beam_daughter_endZ = 0;
  true_beam_daughter_startP = 0;
  true_beam_daughter_startPx = 0;
  true_beam_daughter_startPy = 0;
  true_beam_daughter_startPz = 0;
  true_beam_daughter_Process = 0;
  true_beam_daughter_endProcess = 0;
  true_beam_grand_daughter_ID = 0;
  true_beam_grand_daughter_parID = 0;
  true_beam_grand_daughter_PDG = 0;
  true_beam_grand_daughter_Process = 0;
  true_beam_grand_daughter_endProcess = 0;
  true_beam_Pi0_decay_ID = 0;
  true_beam_Pi0_decay_parID = 0;
  true_beam_Pi0_decay_PDG = 0;
  true_beam_Pi0_decay_startX = 0;
  true_beam_Pi0_decay_startY = 0;
  true_beam_Pi0_decay_startZ = 0;
  true_beam_Pi0_decay_startP = 0;
  true_beam_Pi0_decay_startPx = 0;
  true_beam_Pi0_decay_startPy = 0;
  true_beam_Pi0_decay_startPz = 0;
  true_beam_Pi0_decay_len = 0;
  reco_daughter_PFP_ID = 0;
  reco_daughter_PFP_nHits = 0;
  reco_daughter_PFP_trackScore = 0;
  reco_daughter_PFP_emScore = 0;
  reco_daughter_PFP_michelScore = 0;
  reco_daughter_PFP_nHits_collection = 0;
  reco_daughter_PFP_trackScore_collection = 0;
  reco_daughter_PFP_emScore_collection = 0;
  reco_daughter_PFP_michelScore_collection = 0;
  reco_daughter_allTrack_ID = 0;
  reco_daughter_allTrack_Theta = 0;
  reco_daughter_allTrack_Phi = 0;
  reco_daughter_allTrack_len = 0;
  reco_daughter_allTrack_alt_len = 0;
  reco_daughter_allTrack_momByRange_proton = 0;
  reco_daughter_allTrack_momByRange_muon = 0;
  reco_daughter_allTrack_momByRange_alt_proton = 0;
  reco_daughter_allTrack_momByRange_alt_muon = 0;
  reco_daughter_allTrack_startX = 0;
  reco_daughter_allTrack_startY = 0;
  reco_daughter_allTrack_startZ = 0;
  reco_daughter_allTrack_endX = 0;
  reco_daughter_allTrack_endY = 0;
  reco_daughter_allTrack_endZ = 0;
  reco_daughter_allTrack_vertex_michel_score = 0;
  reco_daughter_allTrack_vertex_nHits = 0;
  reco_daughter_allTrack_dQdX_SCE = 0;
  reco_daughter_allTrack_dEdX_SCE = 0;
  reco_daughter_allTrack_calibrated_dEdX_SCE = 0;
  reco_daughter_allTrack_resRange_SCE = 0;
  reco_daughter_allTrack_Chi2_proton = 0;
  reco_daughter_allTrack_Chi2_ndof = 0;
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane0 = 0;
  reco_daughter_allTrack_resRange_plane0 = 0;
  reco_daughter_allTrack_Chi2_proton_plane0 = 0;
  reco_daughter_allTrack_Chi2_ndof_plane0 = 0;
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane1 = 0;
  reco_daughter_allTrack_resRange_plane1 = 0;
  reco_daughter_allTrack_Chi2_proton_plane1 = 0;
  reco_daughter_allTrack_Chi2_ndof_plane1 = 0;
  reco_daughter_allShower_ID = 0;
  reco_daughter_allShower_len = 0;
  reco_daughter_allShower_startX = 0;
  reco_daughter_allShower_startY = 0;
  reco_daughter_allShower_startZ = 0;
  reco_daughter_allShower_dirX = 0;
  reco_daughter_allShower_dirY = 0;
  reco_daughter_allShower_dirZ = 0;
  reco_daughter_allShower_energy = 0;
  reco_daughter_PFP_true_byHits_PDG = 0;
  reco_daughter_PFP_true_byHits_ID = 0;
  reco_daughter_PFP_true_byHits_origin = 0;
  reco_daughter_PFP_true_byHits_parID = 0;
  reco_daughter_PFP_true_byHits_parPDG = 0;
  reco_daughter_PFP_true_byHits_process = 0;
  reco_daughter_PFP_true_byHits_sharedHits = 0;
  reco_daughter_PFP_true_byHits_emHits = 0;
  reco_daughter_PFP_true_byHits_len = 0;
  reco_daughter_PFP_true_byHits_startX = 0;
  reco_daughter_PFP_true_byHits_startY = 0;
  reco_daughter_PFP_true_byHits_startZ = 0;
  reco_daughter_PFP_true_byHits_endX = 0;
  reco_daughter_PFP_true_byHits_endY = 0;
  reco_daughter_PFP_true_byHits_endZ = 0;
  reco_daughter_PFP_true_byHits_startPx = 0;
  reco_daughter_PFP_true_byHits_startPy = 0;
  reco_daughter_PFP_true_byHits_startPz = 0;
  reco_daughter_PFP_true_byHits_startP = 0;
  reco_daughter_PFP_true_byHits_startE = 0;
  reco_daughter_PFP_true_byHits_endProcess = 0;
  reco_daughter_PFP_true_byHits_purity = 0;
  reco_daughter_PFP_true_byHits_completeness = 0;
  reco_daughter_PFP_true_byE_PDG = 0;
  reco_daughter_PFP_true_byE_len = 0;
  reco_daughter_PFP_true_byE_purity = 0;
  reco_daughter_PFP_true_byE_completeness = 0; 
  beam_inst_P = 0;
  beam_inst_TOF = 0;
  beam_inst_TOF_Chan = 0;
  beam_inst_X = 0;
  beam_inst_Y = 0;
  beam_inst_Z = 0;
  beam_inst_dirX = 0;
  beam_inst_dirY = 0;
  beam_inst_dirZ = 0;
  beam_inst_PDG_candidates = 0;
  beam_inst_nFibersP1 = 0;
  beam_inst_nFibersP2 = 0;
  beam_inst_nFibersP3 = 0;
  beam_inst_nMomenta = 0;
  beam_inst_valid = 0;                                                                
}

//*****************************************************************************
void PDSPAnalyzerTreeConverter::SetBranchAddresses(){
//*****************************************************************************

  fChain->SetBranchAddress("run",&run,&b_run);
  fChain->SetBranchAddress("subrun",&subrun,&b_subrun);
  fChain->SetBranchAddress("event",&event,&b_event);
  fChain->SetBranchAddress("MC",&MC,&b_MC);
  fChain->SetBranchAddress("true_beam_PDG",&true_beam_PDG,&b_true_beam_PDG);
  fChain->SetBranchAddress("true_beam_mass",&true_beam_mass,&b_true_beam_mass);
  fChain->SetBranchAddress("true_beam_ID",&true_beam_ID,&b_true_beam_ID);
  fChain->SetBranchAddress("true_beam_endProcess",&true_beam_endProcess,&b_true_beam_endProcess);
  fChain->SetBranchAddress("true_beam_startX",&true_beam_startX,&b_true_beam_startX);
  fChain->SetBranchAddress("true_beam_startY",&true_beam_startY,&b_true_beam_startY);
  fChain->SetBranchAddress("true_beam_startZ",&true_beam_startZ,&b_true_beam_startZ);
  fChain->SetBranchAddress("true_beam_endX",&true_beam_endX,&b_true_beam_endX);
  fChain->SetBranchAddress("true_beam_endY",&true_beam_endY,&b_true_beam_endY);
  fChain->SetBranchAddress("true_beam_endZ",&true_beam_endZ,&b_true_beam_endZ);
  fChain->SetBranchAddress("true_beam_startP",&true_beam_startP,&b_true_beam_startP);
  fChain->SetBranchAddress("true_beam_startPx",&true_beam_startPx,&b_true_beam_startPx);
  fChain->SetBranchAddress("true_beam_startPy",&true_beam_startPy,&b_true_beam_startPy);
  fChain->SetBranchAddress("true_beam_startPz",&true_beam_startPz,&b_true_beam_startPz);
  fChain->SetBranchAddress("true_beam_endP",&true_beam_endP,&b_true_beam_endP);
  fChain->SetBranchAddress("true_beam_endPx",&true_beam_endPx,&b_true_beam_endPx);
  fChain->SetBranchAddress("true_beam_endPy",&true_beam_endPy,&b_true_beam_endPy);
  fChain->SetBranchAddress("true_beam_endPz",&true_beam_endPz,&b_true_beam_endPz);
  fChain->SetBranchAddress("true_beam_last_len",&true_beam_last_len,&b_true_beam_last_len);
  fChain->SetBranchAddress("true_beam_startDirX",&true_beam_startDirX,&b_true_beam_startDirX);
  fChain->SetBranchAddress("true_beam_startDirY",&true_beam_startDirY,&b_true_beam_startDirY);
  fChain->SetBranchAddress("true_beam_startDirZ",&true_beam_startDirZ,&b_true_beam_startDirZ);
  fChain->SetBranchAddress("true_beam_nElasticScatters",&true_beam_nElasticScatters,&b_true_beam_nElasticScatters);
  fChain->SetBranchAddress("true_beam_elastic_costheta",&true_beam_elastic_costheta,&b_true_beam_elastic_costheta);
  fChain->SetBranchAddress("true_beam_elastic_X",&true_beam_elastic_X,&b_true_beam_elastic_X);
  fChain->SetBranchAddress("true_beam_elastic_Y",&true_beam_elastic_Y,&b_true_beam_elastic_Y);
  fChain->SetBranchAddress("true_beam_elastic_Z",&true_beam_elastic_Z,&b_true_beam_elastic_Z);
  fChain->SetBranchAddress("true_beam_elastic_deltaE",&true_beam_elastic_deltaE,&b_true_beam_elastic_deltaE);
  fChain->SetBranchAddress("true_beam_elastic_IDE_edep",&true_beam_elastic_IDE_edep,&b_true_beam_elastic_IDE_edep);
  fChain->SetBranchAddress("true_beam_IDE_totalDep",&true_beam_IDE_totalDep,&b_true_beam_IDE_totalDep);
  fChain->SetBranchAddress("true_daughter_nPi0",&true_daughter_nPi0,&b_true_daughter_nPi0);
  fChain->SetBranchAddress("true_daughter_nPiPlus",&true_daughter_nPiPlus,&b_true_daughter_nPiPlus);
  fChain->SetBranchAddress("true_daughter_nProton",&true_daughter_nProton,&b_true_daughter_nProton);
  fChain->SetBranchAddress("true_daughter_nNeutron",&true_daughter_nNeutron,&b_true_daughter_nNeutron);
  fChain->SetBranchAddress("true_daughter_nPiMinus",&true_daughter_nPiMinus,&b_true_daughter_nPiMinus);
  fChain->SetBranchAddress("true_daughter_nNucleus",&true_daughter_nNucleus,&b_true_daughter_nNucleus);
  fChain->SetBranchAddress("true_beam_processes",&true_beam_processes,&b_true_beam_processes);
  fChain->SetBranchAddress("true_beam_incidentEnergies",&true_beam_incidentEnergies,&b_true_beam_incidentEnergies);
  fChain->SetBranchAddress("true_beam_interactingEnergy",&true_beam_interactingEnergy,&b_true_beam_interactingEnergy);
  fChain->SetBranchAddress("true_beam_slices",&true_beam_slices,&b_true_beam_slices);
  fChain->SetBranchAddress("true_beam_slices_deltaE",&true_beam_slices_deltaE,&b_true_beam_slices_deltaE);
  fChain->SetBranchAddress("true_beam_traj_X",&true_beam_traj_X,&b_true_beam_traj_X);
  fChain->SetBranchAddress("true_beam_traj_Y",&true_beam_traj_Y,&b_true_beam_traj_Y);
  fChain->SetBranchAddress("true_beam_traj_Z",&true_beam_traj_Z,&b_true_beam_traj_Z);
  fChain->SetBranchAddress("true_beam_traj_KE",&true_beam_traj_KE,&b_true_beam_traj_KE);
  fChain->SetBranchAddress("reco_beam_type",&reco_beam_type,&b_reco_beam_type);
  fChain->SetBranchAddress("reco_beam_trackID",&reco_beam_trackID,&b_reco_beam_trackID);
  fChain->SetBranchAddress("reco_beam_startX",&reco_beam_startX,&b_reco_beam_startX);
  fChain->SetBranchAddress("reco_beam_startY",&reco_beam_startY,&b_reco_beam_startY);
  fChain->SetBranchAddress("reco_beam_startZ",&reco_beam_startZ,&b_reco_beam_startZ);
  fChain->SetBranchAddress("reco_beam_endX",&reco_beam_endX,&b_reco_beam_endX);
  fChain->SetBranchAddress("reco_beam_endY",&reco_beam_endY,&b_reco_beam_endY);
  fChain->SetBranchAddress("reco_beam_endZ",&reco_beam_endZ,&b_reco_beam_endZ);
  fChain->SetBranchAddress("reco_beam_trackDirX",&reco_beam_trackDirX,&b_reco_beam_trackDirX);
  fChain->SetBranchAddress("reco_beam_trackDirY",&reco_beam_trackDirY,&b_reco_beam_trackDirY);
  fChain->SetBranchAddress("reco_beam_trackDirZ",&reco_beam_trackDirZ,&b_reco_beam_trackDirZ);
  fChain->SetBranchAddress("reco_beam_trackEndDirX",&reco_beam_trackEndDirX,&b_reco_beam_trackEndDirX);
  fChain->SetBranchAddress("reco_beam_trackEndDirY",&reco_beam_trackEndDirY,&b_reco_beam_trackEndDirY);
  fChain->SetBranchAddress("reco_beam_trackEndDirZ",&reco_beam_trackEndDirZ,&b_reco_beam_trackEndDirZ);
  fChain->SetBranchAddress("reco_beam_flipped",&reco_beam_flipped,&b_reco_beam_flipped);
  fChain->SetBranchAddress("reco_beam_len",&reco_beam_len,&b_reco_beam_len);
  fChain->SetBranchAddress("reco_beam_alt_len",&reco_beam_alt_len,&b_reco_beam_alt_len);
  fChain->SetBranchAddress("reco_beam_momByRange_proton",&reco_beam_momByRange_proton,&b_reco_beam_momByRange_proton);
  fChain->SetBranchAddress("reco_beam_momByRange_muon",&reco_beam_momByRange_muon,&b_reco_beam_momByRange_muon);
  fChain->SetBranchAddress("reco_beam_momByRange_alt_proton",&reco_beam_momByRange_alt_proton,&b_reco_beam_momByRange_alt_proton);
  fChain->SetBranchAddress("reco_beam_momByRange_alt_muon",&reco_beam_momByRange_alt_muon,&b_reco_beam_momByRange_alt_muon);
  fChain->SetBranchAddress("reco_beam_calo_startX",&reco_beam_calo_startX,&b_reco_beam_calo_startX);
  fChain->SetBranchAddress("reco_beam_calo_startY",&reco_beam_calo_startY,&b_reco_beam_calo_startY);
  fChain->SetBranchAddress("reco_beam_calo_startZ",&reco_beam_calo_startZ,&b_reco_beam_calo_startZ);
  fChain->SetBranchAddress("reco_beam_calo_endX",&reco_beam_calo_endX,&b_reco_beam_calo_endX);
  fChain->SetBranchAddress("reco_beam_calo_endY",&reco_beam_calo_endY,&b_reco_beam_calo_endY);
  fChain->SetBranchAddress("reco_beam_calo_endZ",&reco_beam_calo_endZ,&b_reco_beam_calo_endZ);
  fChain->SetBranchAddress("reco_beam_calo_startDirX",&reco_beam_calo_startDirX,&b_reco_beam_calo_startDirX);
  fChain->SetBranchAddress("reco_beam_calo_startDirY",&reco_beam_calo_startDirY,&b_reco_beam_calo_startDirY);
  fChain->SetBranchAddress("reco_beam_calo_startDirZ",&reco_beam_calo_startDirZ,&b_reco_beam_calo_startDirZ);
  fChain->SetBranchAddress("reco_beam_calo_endDirX",&reco_beam_calo_endDirX,&b_reco_beam_calo_endDirX);
  fChain->SetBranchAddress("reco_beam_calo_endDirY",&reco_beam_calo_endDirY,&b_reco_beam_calo_endDirY);
  fChain->SetBranchAddress("reco_beam_calo_endDirZ",&reco_beam_calo_endDirZ,&b_reco_beam_calo_endDirZ);
  fChain->SetBranchAddress("reco_beam_vertex_michel_score",&reco_beam_vertex_michel_score,&b_reco_beam_vertex_michel_score);
  fChain->SetBranchAddress("reco_beam_vertex_nHits",&reco_beam_vertex_nHits,&b_reco_beam_vertex_nHits);
  fChain->SetBranchAddress("reco_beam_dQdX_SCE",&reco_beam_dQdX_SCE,&b_reco_beam_dQdX_SCE);
  fChain->SetBranchAddress("reco_beam_dQ",&reco_beam_dQ,&b_reco_beam_dQ);
  fChain->SetBranchAddress("reco_beam_dEdX_SCE",&reco_beam_dEdX_SCE,&b_reco_beam_dEdX_SCE);
  fChain->SetBranchAddress("reco_beam_calibrated_dEdX_SCE",&reco_beam_calibrated_dEdX_SCE,&b_reco_beam_calibrated_dEdX_SCE);
  fChain->SetBranchAddress("reco_beam_resRange_SCE",&reco_beam_resRange_SCE,&b_reco_beam_resRange_SCE);
  fChain->SetBranchAddress("reco_beam_TrkPitch_SCE",&reco_beam_TrkPitch_SCE,&b_reco_beam_TrkPitch_SCE);
  fChain->SetBranchAddress("reco_beam_dQdX_NoSCE",&reco_beam_dQdX_NoSCE,&b_reco_beam_dQdX_NoSCE);
  fChain->SetBranchAddress("reco_beam_dEdX_NoSCE",&reco_beam_dEdX_NoSCE,&b_reco_beam_dEdX_NoSCE);
  fChain->SetBranchAddress("reco_beam_calibrated_dEdX_NoSCE",&reco_beam_calibrated_dEdX_NoSCE,&b_reco_beam_calibrated_dEdX_NoSCE);
  fChain->SetBranchAddress("reco_beam_resRange_NoSCE",&reco_beam_resRange_NoSCE,&b_reco_beam_resRange_NoSCE);
  fChain->SetBranchAddress("reco_beam_TrkPitch_NoSCE",&reco_beam_TrkPitch_NoSCE,&b_reco_beam_TrkPitch_NoSCE);
  fChain->SetBranchAddress("reco_beam_calo_wire",&reco_beam_calo_wire,&b_reco_beam_calo_wire);
  fChain->SetBranchAddress("reco_beam_calo_wire_z",&reco_beam_calo_wire_z,&b_reco_beam_calo_wire_z);
  fChain->SetBranchAddress("reco_beam_calo_tick",&reco_beam_calo_tick,&b_reco_beam_calo_tick);
  fChain->SetBranchAddress("reco_beam_calo_TPC",&reco_beam_calo_TPC,&b_reco_beam_calo_TPC);
  fChain->SetBranchAddress("reco_beam_passes_beam_cuts",&reco_beam_passes_beam_cuts,&b_reco_beam_passes_beam_cuts);
  fChain->SetBranchAddress("reco_beam_Chi2_proton",&reco_beam_Chi2_proton,&b_reco_beam_Chi2_proton);
  fChain->SetBranchAddress("reco_beam_Chi2_ndof",&reco_beam_Chi2_ndof,&b_reco_beam_Chi2_ndof);
  fChain->SetBranchAddress("reco_beam_PFP_ID",&reco_beam_PFP_ID,&b_reco_beam_PFP_ID);
  fChain->SetBranchAddress("reco_beam_PFP_nHits",&reco_beam_PFP_nHits,&b_reco_beam_PFP_nHits);
  fChain->SetBranchAddress("reco_beam_PFP_trackScore",&reco_beam_PFP_trackScore,&b_reco_beam_PFP_trackScore);
  fChain->SetBranchAddress("reco_beam_PFP_emScore",&reco_beam_PFP_emScore,&b_reco_beam_PFP_emScore);
  fChain->SetBranchAddress("reco_beam_PFP_michelScore",&reco_beam_PFP_michelScore,&b_reco_beam_PFP_michelScore);
  fChain->SetBranchAddress("reco_beam_PFP_trackScore_collection",&reco_beam_PFP_trackScore_collection,&b_reco_beam_PFP_trackScore_collection);
  fChain->SetBranchAddress("reco_beam_PFP_emScore_collection",&reco_beam_PFP_emScore_collection,&b_reco_beam_PFP_emScore_collection);
  fChain->SetBranchAddress("reco_beam_PFP_michelScore_collection",&reco_beam_PFP_michelScore_collection,&b_reco_beam_PFP_michelScore_collection);
  fChain->SetBranchAddress("reco_beam_allTrack_ID",&reco_beam_allTrack_ID,&b_reco_beam_allTrack_ID);
  fChain->SetBranchAddress("reco_beam_allTrack_beam_cuts",&reco_beam_allTrack_beam_cuts,&b_reco_beam_allTrack_beam_cuts);
  fChain->SetBranchAddress("reco_beam_allTrack_flipped",&reco_beam_allTrack_flipped,&b_reco_beam_allTrack_flipped);
  fChain->SetBranchAddress("reco_beam_allTrack_len",&reco_beam_allTrack_len,&b_reco_beam_allTrack_len);
  fChain->SetBranchAddress("reco_beam_allTrack_startX",&reco_beam_allTrack_startX,&b_reco_beam_allTrack_startX);
  fChain->SetBranchAddress("reco_beam_allTrack_startY",&reco_beam_allTrack_startY,&b_reco_beam_allTrack_startY);
  fChain->SetBranchAddress("reco_beam_allTrack_startZ",&reco_beam_allTrack_startZ,&b_reco_beam_allTrack_startZ);
  fChain->SetBranchAddress("reco_beam_allTrack_endX",&reco_beam_allTrack_endX,&b_reco_beam_allTrack_endX);
  fChain->SetBranchAddress("reco_beam_allTrack_endY",&reco_beam_allTrack_endY,&b_reco_beam_allTrack_endY);
  fChain->SetBranchAddress("reco_beam_allTrack_endZ",&reco_beam_allTrack_endZ,&b_reco_beam_allTrack_endZ);
  fChain->SetBranchAddress("reco_beam_allTrack_trackDirX",&reco_beam_allTrack_trackDirX,&b_reco_beam_allTrack_trackDirX);
  fChain->SetBranchAddress("reco_beam_allTrack_trackDirY",&reco_beam_allTrack_trackDirY,&b_reco_beam_allTrack_trackDirY);
  fChain->SetBranchAddress("reco_beam_allTrack_trackDirZ",&reco_beam_allTrack_trackDirZ,&b_reco_beam_allTrack_trackDirZ);
  fChain->SetBranchAddress("reco_beam_allTrack_trackEndDirX",&reco_beam_allTrack_trackEndDirX,&b_reco_beam_allTrack_trackEndDirX);
  fChain->SetBranchAddress("reco_beam_allTrack_trackEndDirY",&reco_beam_allTrack_trackEndDirY,&b_reco_beam_allTrack_trackEndDirY);
  fChain->SetBranchAddress("reco_beam_allTrack_trackEndDirZ",&reco_beam_allTrack_trackEndDirZ,&b_reco_beam_allTrack_trackEndDirZ);
  fChain->SetBranchAddress("reco_beam_allTrack_resRange",&reco_beam_allTrack_resRange,&b_reco_beam_allTrack_resRange);
  fChain->SetBranchAddress("reco_beam_allTrack_calibrated_dEdX",&reco_beam_allTrack_calibrated_dEdX,&b_reco_beam_allTrack_calibrated_dEdX);
  fChain->SetBranchAddress("reco_beam_allTrack_Chi2_proton",&reco_beam_allTrack_Chi2_proton,&b_reco_beam_allTrack_Chi2_proton);
  fChain->SetBranchAddress("reco_beam_allTrack_Chi2_ndof",&reco_beam_allTrack_Chi2_ndof,&b_reco_beam_allTrack_Chi2_ndof);
  fChain->SetBranchAddress("reco_beam_true_byHits_endProcess",&reco_beam_true_byHits_endProcess,&b_reco_beam_true_byHits_endProcess);
  fChain->SetBranchAddress("reco_beam_true_byHits_process",&reco_beam_true_byHits_process,&b_reco_beam_true_byHits_process);
  fChain->SetBranchAddress("reco_beam_true_byHits_origin",&reco_beam_true_byHits_origin,&b_reco_beam_true_byHits_origin);
  fChain->SetBranchAddress("reco_beam_true_byHits_PDG",&reco_beam_true_byHits_PDG,&b_reco_beam_true_byHits_PDG);
  fChain->SetBranchAddress("reco_beam_true_byHits_ID",&reco_beam_true_byHits_ID,&b_reco_beam_true_byHits_ID);
  fChain->SetBranchAddress("reco_beam_true_byHits_matched",&reco_beam_true_byHits_matched,&b_reco_beam_true_byHits_matched);
  fChain->SetBranchAddress("reco_beam_true_byHits_purity",&reco_beam_true_byHits_purity,&b_reco_beam_true_byHits_purity);
  fChain->SetBranchAddress("reco_beam_true_byHits_startE",&reco_beam_true_byHits_startE,&b_reco_beam_true_byHits_startE);
  fChain->SetBranchAddress("reco_beam_true_byHits_startP",&reco_beam_true_byHits_startP,&b_reco_beam_true_byHits_startP);
  fChain->SetBranchAddress("reco_beam_true_byHits_startPx",&reco_beam_true_byHits_startPx,&b_reco_beam_true_byHits_startPx);
  fChain->SetBranchAddress("reco_beam_true_byHits_startPy",&reco_beam_true_byHits_startPy,&b_reco_beam_true_byHits_startPy);
  fChain->SetBranchAddress("reco_beam_true_byHits_startPz",&reco_beam_true_byHits_startPz,&b_reco_beam_true_byHits_startPz);
  fChain->SetBranchAddress("reco_beam_true_byHits_endE",&reco_beam_true_byHits_endE,&b_reco_beam_true_byHits_endE);
  fChain->SetBranchAddress("reco_beam_true_byHits_endP",&reco_beam_true_byHits_endP,&b_reco_beam_true_byHits_endP);
  fChain->SetBranchAddress("reco_beam_true_byHits_endPx",&reco_beam_true_byHits_endPx,&b_reco_beam_true_byHits_endPx);
  fChain->SetBranchAddress("reco_beam_true_byHits_endPy",&reco_beam_true_byHits_endPy,&b_reco_beam_true_byHits_endPy);
  fChain->SetBranchAddress("reco_beam_true_byHits_endPz",&reco_beam_true_byHits_endPz,&b_reco_beam_true_byHits_endPz);
  fChain->SetBranchAddress("reco_beam_true_byE_startE",&reco_beam_true_byE_startE,&b_reco_beam_true_byE_startE);
  fChain->SetBranchAddress("reco_beam_true_byE_startP",&reco_beam_true_byE_startP,&b_reco_beam_true_byE_startP);
  fChain->SetBranchAddress("reco_beam_true_byE_startPx",&reco_beam_true_byE_startPx,&b_reco_beam_true_byE_startPx);
  fChain->SetBranchAddress("reco_beam_true_byE_startPy",&reco_beam_true_byE_startPy,&b_reco_beam_true_byE_startPy);
  fChain->SetBranchAddress("reco_beam_true_byE_startPz",&reco_beam_true_byE_startPz,&b_reco_beam_true_byE_startPz);
  fChain->SetBranchAddress("reco_beam_true_byE_endE",&reco_beam_true_byE_endE,&b_reco_beam_true_byE_endE);
  fChain->SetBranchAddress("reco_beam_true_byE_endP",&reco_beam_true_byE_endP,&b_reco_beam_true_byE_endP);
  fChain->SetBranchAddress("reco_beam_true_byE_endPx",&reco_beam_true_byE_endPx,&b_reco_beam_true_byE_endPx);
  fChain->SetBranchAddress("reco_beam_true_byE_endPy",&reco_beam_true_byE_endPy,&b_reco_beam_true_byE_endPy);
  fChain->SetBranchAddress("reco_beam_true_byE_endPz",&reco_beam_true_byE_endPz,&b_reco_beam_true_byE_endPz);
  fChain->SetBranchAddress("reco_beam_true_byE_endProcess",&reco_beam_true_byE_endProcess,&b_reco_beam_true_byE_endProcess);
  fChain->SetBranchAddress("reco_beam_true_byE_process",&reco_beam_true_byE_process,&b_reco_beam_true_byE_process);
  fChain->SetBranchAddress("reco_beam_true_byE_origin",&reco_beam_true_byE_origin,&b_reco_beam_true_byE_origin);
  fChain->SetBranchAddress("reco_beam_true_byE_PDG",&reco_beam_true_byE_PDG,&b_reco_beam_true_byE_PDG);
  fChain->SetBranchAddress("reco_beam_true_byE_ID",&reco_beam_true_byE_ID,&b_reco_beam_true_byE_ID);
  fChain->SetBranchAddress("reco_beam_true_byE_matched",&reco_beam_true_byE_matched,&b_reco_beam_true_byE_matched);
  fChain->SetBranchAddress("reco_beam_incidentEnergies",&reco_beam_incidentEnergies,&b_reco_beam_incidentEnergies);
  fChain->SetBranchAddress("reco_beam_interactingEnergy",&reco_beam_interactingEnergy,&b_reco_beam_interactingEnergy);
  fChain->SetBranchAddress("reco_track_startX",&reco_track_startX,&b_reco_track_startX);
  fChain->SetBranchAddress("reco_track_startY",&reco_track_startY,&b_reco_track_startY);
  fChain->SetBranchAddress("reco_track_startZ",&reco_track_startZ,&b_reco_track_startZ);
  fChain->SetBranchAddress("reco_track_endX",&reco_track_endX,&b_reco_track_endX);
  fChain->SetBranchAddress("reco_track_endY",&reco_track_endY,&b_reco_track_endY);
  fChain->SetBranchAddress("reco_track_endZ",&reco_track_endZ,&b_reco_track_endZ);
  fChain->SetBranchAddress("reco_track_michel_score",&reco_track_michel_score,&b_reco_track_michel_score);
  fChain->SetBranchAddress("reco_track_nHits",&reco_track_nHits,&b_reco_track_nHits);
  fChain->SetBranchAddress("reco_track_ID",&reco_track_ID,&b_reco_track_ID);
  fChain->SetBranchAddress("true_beam_daughter_ID",&true_beam_daughter_ID,&b_true_beam_daughter_ID);
  fChain->SetBranchAddress("true_beam_daughter_PDG",&true_beam_daughter_PDG,&b_true_beam_daughter_PDG);
  fChain->SetBranchAddress("true_beam_daughter_len",&true_beam_daughter_len,&b_true_beam_daughter_len);
  fChain->SetBranchAddress("true_beam_daughter_startX",&true_beam_daughter_startX,&b_true_beam_daughter_startX);
  fChain->SetBranchAddress("true_beam_daughter_startY",&true_beam_daughter_startY,&b_true_beam_daughter_startY);
  fChain->SetBranchAddress("true_beam_daughter_startZ",&true_beam_daughter_startZ,&b_true_beam_daughter_startZ);
  fChain->SetBranchAddress("true_beam_daughter_endX",&true_beam_daughter_endX,&b_true_beam_daughter_endX);
  fChain->SetBranchAddress("true_beam_daughter_endY",&true_beam_daughter_endY,&b_true_beam_daughter_endY);
  fChain->SetBranchAddress("true_beam_daughter_endZ",&true_beam_daughter_endZ,&b_true_beam_daughter_endZ);
  fChain->SetBranchAddress("true_beam_daughter_startP",&true_beam_daughter_startP,&b_true_beam_daughter_startP);
  fChain->SetBranchAddress("true_beam_daughter_startPx",&true_beam_daughter_startPx,&b_true_beam_daughter_startPx);
  fChain->SetBranchAddress("true_beam_daughter_startPy",&true_beam_daughter_startPy,&b_true_beam_daughter_startPy);
  fChain->SetBranchAddress("true_beam_daughter_startPz",&true_beam_daughter_startPz,&b_true_beam_daughter_startPz);
  fChain->SetBranchAddress("true_beam_daughter_Process",&true_beam_daughter_Process,&b_true_beam_daughter_Process);
  fChain->SetBranchAddress("true_beam_daughter_endProcess",&true_beam_daughter_endProcess,&b_true_beam_daughter_endProcess);
  fChain->SetBranchAddress("true_beam_grand_daughter_ID",&true_beam_grand_daughter_ID,&b_true_beam_grand_daughter_ID);
  fChain->SetBranchAddress("true_beam_grand_daughter_parID",&true_beam_grand_daughter_parID,&b_true_beam_grand_daughter_parID);
  fChain->SetBranchAddress("true_beam_grand_daughter_PDG",&true_beam_grand_daughter_PDG,&b_true_beam_grand_daughter_PDG);
  fChain->SetBranchAddress("true_beam_grand_daughter_Process",&true_beam_grand_daughter_Process,&b_true_beam_grand_daughter_Process);
  fChain->SetBranchAddress("true_beam_grand_daughter_endProcess",&true_beam_grand_daughter_endProcess,&b_true_beam_grand_daughter_endProcess);
  fChain->SetBranchAddress("true_beam_Pi0_decay_ID",&true_beam_Pi0_decay_ID,&b_true_beam_Pi0_decay_ID);
  fChain->SetBranchAddress("true_beam_Pi0_decay_parID",&true_beam_Pi0_decay_parID,&b_true_beam_Pi0_decay_parID);
  fChain->SetBranchAddress("true_beam_Pi0_decay_PDG",&true_beam_Pi0_decay_PDG,&b_true_beam_Pi0_decay_PDG);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startX",&true_beam_Pi0_decay_startX,&b_true_beam_Pi0_decay_startX);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startY",&true_beam_Pi0_decay_startY,&b_true_beam_Pi0_decay_startY);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startZ",&true_beam_Pi0_decay_startZ,&b_true_beam_Pi0_decay_startZ);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startP",&true_beam_Pi0_decay_startP,&b_true_beam_Pi0_decay_startP);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startPx",&true_beam_Pi0_decay_startPx,&b_true_beam_Pi0_decay_startPx);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startPy",&true_beam_Pi0_decay_startPy,&b_true_beam_Pi0_decay_startPy);
  fChain->SetBranchAddress("true_beam_Pi0_decay_startPz",&true_beam_Pi0_decay_startPz,&b_true_beam_Pi0_decay_startPz);
  fChain->SetBranchAddress("true_beam_Pi0_decay_len",&true_beam_Pi0_decay_len,&b_true_beam_Pi0_decay_len);
  fChain->SetBranchAddress("reco_daughter_PFP_ID",&reco_daughter_PFP_ID,&b_reco_daughter_PFP_ID);
  fChain->SetBranchAddress("reco_daughter_PFP_nHits",&reco_daughter_PFP_nHits,&b_reco_daughter_PFP_nHits);
  fChain->SetBranchAddress("reco_daughter_PFP_trackScore",&reco_daughter_PFP_trackScore,&b_reco_daughter_PFP_trackScore);
  fChain->SetBranchAddress("reco_daughter_PFP_emScore",&reco_daughter_PFP_emScore,&b_reco_daughter_PFP_emScore);
  fChain->SetBranchAddress("reco_daughter_PFP_michelScore",&reco_daughter_PFP_michelScore,&b_reco_daughter_PFP_michelScore);
  fChain->SetBranchAddress("reco_daughter_PFP_nHits_collection",&reco_daughter_PFP_nHits_collection,&b_reco_daughter_PFP_nHits_collection);
  fChain->SetBranchAddress("reco_daughter_PFP_trackScore_collection",&reco_daughter_PFP_trackScore_collection,&b_reco_daughter_PFP_trackScore_collection);
  fChain->SetBranchAddress("reco_daughter_PFP_emScore_collection",&reco_daughter_PFP_emScore_collection,&b_reco_daughter_PFP_emScore_collection);
  fChain->SetBranchAddress("reco_daughter_PFP_michelScore_collection",&reco_daughter_PFP_michelScore_collection,&b_reco_daughter_PFP_michelScore_collection);
  fChain->SetBranchAddress("reco_daughter_allTrack_ID",&reco_daughter_allTrack_ID,&b_reco_daughter_allTrack_ID);
  fChain->SetBranchAddress("reco_daughter_allTrack_Theta",&reco_daughter_allTrack_Theta,&b_reco_daughter_allTrack_Theta);
  fChain->SetBranchAddress("reco_daughter_allTrack_Phi",&reco_daughter_allTrack_Phi,&b_reco_daughter_allTrack_Phi);
  fChain->SetBranchAddress("reco_daughter_allTrack_len",&reco_daughter_allTrack_len,&b_reco_daughter_allTrack_len);
  fChain->SetBranchAddress("reco_daughter_allTrack_alt_len",&reco_daughter_allTrack_alt_len,&b_reco_daughter_allTrack_alt_len);
  fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_proton",&reco_daughter_allTrack_momByRange_proton,&b_reco_daughter_allTrack_momByRange_proton);
  fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_muon",&reco_daughter_allTrack_momByRange_muon,&b_reco_daughter_allTrack_momByRange_muon);
  fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_alt_proton",&reco_daughter_allTrack_momByRange_alt_proton,&b_reco_daughter_allTrack_momByRange_alt_proton);
  fChain->SetBranchAddress("reco_daughter_allTrack_momByRange_alt_muon",&reco_daughter_allTrack_momByRange_alt_muon,&b_reco_daughter_allTrack_momByRange_alt_muon);
  fChain->SetBranchAddress("reco_daughter_allTrack_startX",&reco_daughter_allTrack_startX,&b_reco_daughter_allTrack_startX);
  fChain->SetBranchAddress("reco_daughter_allTrack_startY",&reco_daughter_allTrack_startY,&b_reco_daughter_allTrack_startY);
  fChain->SetBranchAddress("reco_daughter_allTrack_startZ",&reco_daughter_allTrack_startZ,&b_reco_daughter_allTrack_startZ);
  fChain->SetBranchAddress("reco_daughter_allTrack_endX",&reco_daughter_allTrack_endX,&b_reco_daughter_allTrack_endX);
  fChain->SetBranchAddress("reco_daughter_allTrack_endY",&reco_daughter_allTrack_endY,&b_reco_daughter_allTrack_endY);
  fChain->SetBranchAddress("reco_daughter_allTrack_endZ",&reco_daughter_allTrack_endZ,&b_reco_daughter_allTrack_endZ);
  fChain->SetBranchAddress("reco_daughter_allTrack_vertex_michel_score",&reco_daughter_allTrack_vertex_michel_score,&b_reco_daughter_allTrack_vertex_michel_score);
  fChain->SetBranchAddress("reco_daughter_allTrack_vertex_nHits",&reco_daughter_allTrack_vertex_nHits,&b_reco_daughter_allTrack_vertex_nHits);
  fChain->SetBranchAddress("reco_daughter_allTrack_dQdX_SCE",&reco_daughter_allTrack_dQdX_SCE,&b_reco_daughter_allTrack_dQdX_SCE);
  fChain->SetBranchAddress("reco_daughter_allTrack_dEdX_SCE",&reco_daughter_allTrack_dEdX_SCE,&b_reco_daughter_allTrack_dEdX_SCE);
  fChain->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE",&reco_daughter_allTrack_calibrated_dEdX_SCE,&b_reco_daughter_allTrack_calibrated_dEdX_SCE);
  fChain->SetBranchAddress("reco_daughter_allTrack_resRange_SCE",&reco_daughter_allTrack_resRange_SCE,&b_reco_daughter_allTrack_resRange_SCE);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_proton",&reco_daughter_allTrack_Chi2_proton,&b_reco_daughter_allTrack_Chi2_proton);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof",&reco_daughter_allTrack_Chi2_ndof,&b_reco_daughter_allTrack_Chi2_ndof);
  fChain->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE_plane0",&reco_daughter_allTrack_calibrated_dEdX_SCE_plane0,&b_reco_daughter_allTrack_calibrated_dEdX_SCE_plane0);
  fChain->SetBranchAddress("reco_daughter_allTrack_resRange_plane0",&reco_daughter_allTrack_resRange_plane0,&b_reco_daughter_allTrack_resRange_plane0);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_proton_plane0",&reco_daughter_allTrack_Chi2_proton_plane0,&b_reco_daughter_allTrack_Chi2_proton_plane0);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof_plane0",&reco_daughter_allTrack_Chi2_ndof_plane0,&b_reco_daughter_allTrack_Chi2_ndof_plane0);
  fChain->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE_plane1",&reco_daughter_allTrack_calibrated_dEdX_SCE_plane1,&b_reco_daughter_allTrack_calibrated_dEdX_SCE_plane1);
  fChain->SetBranchAddress("reco_daughter_allTrack_resRange_plane1",&reco_daughter_allTrack_resRange_plane1,&b_reco_daughter_allTrack_resRange_plane1);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_proton_plane1",&reco_daughter_allTrack_Chi2_proton_plane1,&b_reco_daughter_allTrack_Chi2_proton_plane1);
  fChain->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof_plane1",&reco_daughter_allTrack_Chi2_ndof_plane1,&b_reco_daughter_allTrack_Chi2_ndof_plane1);
  fChain->SetBranchAddress("reco_daughter_allShower_ID",&reco_daughter_allShower_ID,&b_reco_daughter_allShower_ID);
  fChain->SetBranchAddress("reco_daughter_allShower_len",&reco_daughter_allShower_len,&b_reco_daughter_allShower_len);
  fChain->SetBranchAddress("reco_daughter_allShower_startX",&reco_daughter_allShower_startX,&b_reco_daughter_allShower_startX);
  fChain->SetBranchAddress("reco_daughter_allShower_startY",&reco_daughter_allShower_startY,&b_reco_daughter_allShower_startY);
  fChain->SetBranchAddress("reco_daughter_allShower_startZ",&reco_daughter_allShower_startZ,&b_reco_daughter_allShower_startZ);
  fChain->SetBranchAddress("reco_daughter_allShower_dirX",&reco_daughter_allShower_dirX,&b_reco_daughter_allShower_dirX);
  fChain->SetBranchAddress("reco_daughter_allShower_dirY",&reco_daughter_allShower_dirY,&b_reco_daughter_allShower_dirY);
  fChain->SetBranchAddress("reco_daughter_allShower_dirZ",&reco_daughter_allShower_dirZ,&b_reco_daughter_allShower_dirZ);
  fChain->SetBranchAddress("reco_daughter_allShower_energy",&reco_daughter_allShower_energy,&b_reco_daughter_allShower_energy);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_PDG",&reco_daughter_PFP_true_byHits_PDG,&b_reco_daughter_PFP_true_byHits_PDG);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_ID",&reco_daughter_PFP_true_byHits_ID,&b_reco_daughter_PFP_true_byHits_ID);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_origin",&reco_daughter_PFP_true_byHits_origin,&b_reco_daughter_PFP_true_byHits_origin);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_parID",&reco_daughter_PFP_true_byHits_parID,&b_reco_daughter_PFP_true_byHits_parID);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_parPDG",&reco_daughter_PFP_true_byHits_parPDG,&b_reco_daughter_PFP_true_byHits_parPDG);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_process",&reco_daughter_PFP_true_byHits_process,&b_reco_daughter_PFP_true_byHits_process);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_sharedHits",&reco_daughter_PFP_true_byHits_sharedHits,&b_reco_daughter_PFP_true_byHits_sharedHits);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_emHits",&reco_daughter_PFP_true_byHits_emHits,&b_reco_daughter_PFP_true_byHits_emHits);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_len",&reco_daughter_PFP_true_byHits_len,&b_reco_daughter_PFP_true_byHits_len);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startX",&reco_daughter_PFP_true_byHits_startX,&b_reco_daughter_PFP_true_byHits_startX);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startY",&reco_daughter_PFP_true_byHits_startY,&b_reco_daughter_PFP_true_byHits_startY);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startZ",&reco_daughter_PFP_true_byHits_startZ,&b_reco_daughter_PFP_true_byHits_startZ);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endX",&reco_daughter_PFP_true_byHits_endX,&b_reco_daughter_PFP_true_byHits_endX);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endY",&reco_daughter_PFP_true_byHits_endY,&b_reco_daughter_PFP_true_byHits_endY);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endZ",&reco_daughter_PFP_true_byHits_endZ,&b_reco_daughter_PFP_true_byHits_endZ);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPx",&reco_daughter_PFP_true_byHits_startPx,&b_reco_daughter_PFP_true_byHits_startPx);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPy",&reco_daughter_PFP_true_byHits_startPy,&b_reco_daughter_PFP_true_byHits_startPy);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startPz",&reco_daughter_PFP_true_byHits_startPz,&b_reco_daughter_PFP_true_byHits_startPz);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startP",&reco_daughter_PFP_true_byHits_startP,&b_reco_daughter_PFP_true_byHits_startP);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_startE",&reco_daughter_PFP_true_byHits_startE,&b_reco_daughter_PFP_true_byHits_startE);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_endProcess",&reco_daughter_PFP_true_byHits_endProcess,&b_reco_daughter_PFP_true_byHits_endProcess);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_purity",&reco_daughter_PFP_true_byHits_purity,&b_reco_daughter_PFP_true_byHits_purity);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byHits_completeness",&reco_daughter_PFP_true_byHits_completeness,&b_reco_daughter_PFP_true_byHits_completeness);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byE_PDG",&reco_daughter_PFP_true_byE_PDG,&b_reco_daughter_PFP_true_byE_PDG);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byE_len",&reco_daughter_PFP_true_byE_len,&b_reco_daughter_PFP_true_byE_len);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byE_purity",&reco_daughter_PFP_true_byE_purity,&b_reco_daughter_PFP_true_byE_purity);
  fChain->SetBranchAddress("reco_daughter_PFP_true_byE_completeness",&reco_daughter_PFP_true_byE_completeness,&b_reco_daughter_PFP_true_byE_completeness);
  fChain->SetBranchAddress("beam_inst_P",&beam_inst_P,&b_beam_inst_P);
  fChain->SetBranchAddress("beam_inst_TOF",&beam_inst_TOF,&b_beam_inst_TOF);
  fChain->SetBranchAddress("beam_inst_TOF_Chan",&beam_inst_TOF_Chan,&b_beam_inst_TOF_Chan);
  fChain->SetBranchAddress("beam_inst_X",&beam_inst_X,&b_beam_inst_X);
  fChain->SetBranchAddress("beam_inst_Y",&beam_inst_Y,&b_beam_inst_Y);
  fChain->SetBranchAddress("beam_inst_Z",&beam_inst_Z,&b_beam_inst_Z);
  fChain->SetBranchAddress("beam_inst_dirX",&beam_inst_dirX,&b_beam_inst_dirX);
  fChain->SetBranchAddress("beam_inst_dirY",&beam_inst_dirY,&b_beam_inst_dirY);
  fChain->SetBranchAddress("beam_inst_dirZ",&beam_inst_dirZ,&b_beam_inst_dirZ);
  fChain->SetBranchAddress("beam_inst_PDG_candidates",&beam_inst_PDG_candidates,&b_beam_inst_PDG_candidates);
  fChain->SetBranchAddress("beam_inst_nFibersP1",&beam_inst_nFibersP1,&b_beam_inst_nFibersP1);
  fChain->SetBranchAddress("beam_inst_nFibersP2",&beam_inst_nFibersP2,&b_beam_inst_nFibersP2);
  fChain->SetBranchAddress("beam_inst_nFibersP3",&beam_inst_nFibersP3,&b_beam_inst_nFibersP3);
  fChain->SetBranchAddress("beam_inst_nMomenta",&beam_inst_nMomenta,&b_beam_inst_nMomenta);
  fChain->SetBranchAddress("beam_inst_valid",&beam_inst_valid,&b_beam_inst_valid);
}
