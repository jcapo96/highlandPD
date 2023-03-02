#include "MichelRemovingTreeConverter.hxx"
#include "InputManager.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"


//********************************************************************
MichelRemovingTreeConverter::MichelRemovingTreeConverter(const std::string& name):pdBaseConverter(name){
//********************************************************************

  _XZmin = ND::params().GetParameterI("pdIO.MichelRemoving.ThetaXZ_min");
  _XZmax = ND::params().GetParameterI("pdIO.MichelRemoving.ThetaXZ_max");
  _YZmin = ND::params().GetParameterI("pdIO.MichelRemoving.ThetaYZ_min");
  _YZmax = ND::params().GetParameterI("pdIO.MichelRemoving.ThetaYZ_max");

  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "ANGLE BINNING FOR TRACK SELECTION" << std::endl;
  std::cout << _XZmin << " < thetaXZ <" << _XZmax << std::endl;
  std::cout << 180-_XZmin << " < thetaXZ <" << 180-_XZmax << std::endl;
  std::cout << _YZmin << " < thetaYZ <" << _YZmax << std::endl;
  std::cout << 180-_YZmin << " < thetaYZ <" << 180-_YZmax << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
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
    if(!AngleCut(itrk))continue;
    //create new AnaParticle
    AnaParticlePD* part = MakeParticle();
    //fill particle info
    FillParticleInfo(part,itrk);
    //add particle to vector of particles
    if(part->Hits[2].size()>0)
      bunch->Particles.push_back(part);
    else
      delete part;
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
    if(!IsUsableHit(itrk,ihit) || !IsInterestingHit(itrk,ihit))continue;
    AnaHitPD hit;
    //fill calo info
    hit.dQdx_elife    = trkdqdx[itrk][2][ihit];
    //hit.dEdx          = trkdedx[itrk][2][ihit];
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
bool MichelRemovingTreeConverter::IsInterestingHit(const int itrk, const int ihit){
//*****************************************************************************

  bool ItIs = false;

  if(trkhitx[itrk][2][ihit] > -360 && trkhitx[itrk][2][ihit] < 0 &&
     abs(trkhity[itrk][2][ihit]-300) < 300 &&
     abs(trkhitz[itrk][2][ihit]-350) < 350)
    ItIs = true;

  return ItIs;
}

//*****************************************************************************
bool MichelRemovingTreeConverter::AngleCut(const int itrk){
//*****************************************************************************

  double xz = abs(180/TMath::Pi()*trackthetaxz[itrk]);
  double yz = abs(180/TMath::Pi()*trackthetayz[itrk]);
  if(((xz>_XZmin && xz<_XZmax) || 
      (xz>180-_XZmax && xz<180-_XZmin)) 
     &&
     ((yz>_YZmin && yz<_YZmax) || 
      (yz>180-_YZmax && yz<180-_YZmin)))
    return true;
  else return false;
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
