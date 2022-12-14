#include "MichelRemovingTreeConverter.hxx"
#include "InputManager.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"


//********************************************************************
MichelRemovingTreeConverter::MichelRemovingTreeConverter(const std::string& name):pdBaseConverter(name){
//********************************************************************

  //nothing to do here for the moment
}

//*****************************************************************************
void MichelRemovingTreeConverter::FillEventInfo(AnaEventInfo* info){
//*****************************************************************************

  info->Run    = run;
  info->SubRun = subrun;
  info->Event  = event;
  info->IsMC   = !isData;
}

//*****************************************************************************
void MichelRemovingTreeConverter::FillBunchInfo(std::vector<AnaTrueParticleB*>& trueParticles, AnaBunch* bunch, AnaBeamPD* beam){
//*****************************************************************************

  bunch->Bunch  = 1;
  bunch->Weight = 1;
  bunch->Particles.clear();
  bunch->Vertices.clear();
  
  //loop over tracks
  for(int itrk = 0; itrk < cross_trks; itrk++){
    //create new AnaParticle
    AnaParticlePD* part = MakeParticle();
    //fill particle info
    FillParticleInfo(part,itrk);
    //add particle to vector of particles
    bunch->Particles.push_back(part);
  }
}

//*****************************************************************************
void MichelRemovingTreeConverter::FillParticleInfo(AnaParticlePD* part, const int itrk){
//*****************************************************************************

  part->UniqueID = TrkID[itrk];
  part->Length   = trklen[itrk];
  part->ThetaXZ  = trackthetaxz[itrk];
  part->ThetaYZ  = trackthetayz[itrk];
  part->PositionStart[0] = trkstartx[itrk];
  part->PositionStart[1] = trkstarty[itrk];
  part->PositionStart[2] = trkstartz[itrk];
  part->PositionStart[0] = trkendx[itrk];
  part->PositionStart[1] = trkendy[itrk];
  part->PositionStart[2] = trkendz[itrk];

  //hit info (only collection needed)
  part->NHitsPerPlane[2] = ntrkhits[itrk][2];
  //loop over hits
  for(int ihit = 0; ihit < _MaxHits; ihit++){
    //create an AnaHit
    if(!IsUsableHit(itrk,ihit))continue;
    AnaHitPD hit;
    //fill calo info
    hit.dQdx          = trkdqdx[itrk][2][ihit];
    hit.dEdx          = trkdedx[itrk][2][ihit];
    hit.ResidualRange = trkresrange[itrk][2][ihit];
    hit.Pitch         = trkpitch[itrk][2][ihit];
    hit.Position.SetXYZ(trkhitx[itrk][2][ihit],trkhity[itrk][2][ihit],trkhitz[itrk][2][ihit]);
    
    hit.dQdx_NoSCE          = trkdqdx[itrk][2][ihit];
    hit.Pitch_NoSCE         = trkpitch[itrk][2][ihit];
    hit.Position_NoSCE.SetXYZ(trkhitx[itrk][2][ihit],trkhity[itrk][2][ihit],trkhitz[itrk][2][ihit]);
    
    hit.TPCid               = pdAnaUtils::GetHitTPCid(hit);
    //add hit to vector of hits
    part->Hits[2].push_back(hit);
  }
  pdAnaUtils::EstimateHitsDirection(part);
}

//*****************************************************************************
bool MichelRemovingTreeConverter::IsUsableHit(const int itrk, const int ihit){
//*****************************************************************************

  bool ItIs = true;
  if(trkhitx[itrk][2][ihit] == -99999 || trkhity[itrk][2][ihit] == -99999 || 
     trkhitz[itrk][2][ihit] == -99999 || trkdqdx[itrk][2][ihit] == -99999) ItIs = false;
  
  return ItIs;
}

//*****************************************************************************
void MichelRemovingTreeConverter::InitializeVariables(){
//*****************************************************************************

  isData = false;
  run    = -999;                  
  subrun = -999;               
  event  = -999;
  
  cross_trks = _MaxTracks;
  for(int itrk = 0; itrk < _MaxTracks; itrk ++){
    trackthetaxz[itrk] = -999;
    trackthetayz[itrk] = -999;
    trkstartx[itrk]    = -999;
    trkstarty[itrk]    = -999;
    trkstartz[itrk]    = -999;
    trkendx[itrk]      = -999;
    trkendy[itrk]      = -999;
    trkendz[itrk]      = -999;
    trklen[itrk]       = -999;
    TrkID[itrk]        = -999;
    for(int ipln = 0; ipln < _NPlanes; ipln++){
      ntrkhits[itrk][ipln] = -999;
      for(int ihit = 0; ihit < _MaxHits; ihit++){
	trkdqdx[itrk][ipln][ihit]     = -999;
	trkdedx[itrk][ipln][ihit]     = -999;
	trkresrange[itrk][ipln][ihit] = -999;
	trkhitx[itrk][ipln][ihit]     = -999;
	trkhity[itrk][ipln][ihit]     = -999;
	trkhitz[itrk][ipln][ihit]     = -999;
	trkpitch[itrk][ipln][ihit]    = -999;
      }
    }
  }  
}

//*****************************************************************************
void MichelRemovingTreeConverter::SetBranchAddresses(){
//*****************************************************************************

  fChain->SetBranchAddress("isData" , &isData , &b_isData);
  fChain->SetBranchAddress("event"  , &event  , &b_event);                  
  fChain->SetBranchAddress("run"    , &run    , &b_run);
  fChain->SetBranchAddress("subrun" , &subrun , &b_subrun);

  fChain->SetBranchAddress("cross_trks",   &cross_trks  , &b_cross_trks);
  fChain->SetBranchAddress("trackthetaxz", &trackthetaxz, &b_trackthetaxz);
  fChain->SetBranchAddress("trackthetayz", &trackthetayz, &b_trackthetayz);
  fChain->SetBranchAddress("trkstartx"   , &trkstartx   , &b_trkstartx);
  fChain->SetBranchAddress("trkstarty"   , &trkstarty   , &b_trkstarty);
  fChain->SetBranchAddress("trkstartz"   , &trkstartz   , &b_trkstartz);
  fChain->SetBranchAddress("trkendx"     , &trkendx     , &b_trkendx);
  fChain->SetBranchAddress("trkendy"     , &trkendy     , &b_trkendy);
  fChain->SetBranchAddress("trkendz"     , &trkendz     , &b_trkendz);
  fChain->SetBranchAddress("trklen"      , &trklen      , &b_trklen);
  fChain->SetBranchAddress("TrkID"       , &TrkID       , &b_TrkID); 
  fChain->SetBranchAddress("ntrkhits"    , &ntrkhits    , &b_ntrkhits);
  fChain->SetBranchAddress("trkdqdx"     , &trkdqdx     , &b_trkdqdx);
  fChain->SetBranchAddress("trkdedx"     , &trkdedx     , &b_trkdedx);
  fChain->SetBranchAddress("trkresrange" , &trkresrange , &b_trkresrange);
  fChain->SetBranchAddress("trkhitx"     , &trkhitx     , &b_trkhitx);
  fChain->SetBranchAddress("trkhity"     , &trkhity     , &b_trkhity);
  fChain->SetBranchAddress("trkhitz"     , &trkhitz     , &b_trkhitz);
  fChain->SetBranchAddress("trkpitch"    , &trkpitch    , &b_trkpitch);
}
