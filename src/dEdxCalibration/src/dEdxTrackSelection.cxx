#include "dEdxTrackSelection.hxx"
#include "EventBoxdEdx.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
dEdxTrackSelection::dEdxTrackSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxdEdx) {
//********************************************************************

}

//********************************************************************
void dEdxTrackSelection::DefineSteps(){
//********************************************************************

  AddStep(StepBase::kAction, "get vector of stopping tracks", new GetNonStoppingTracksAction(), true);
  AddStep(StepBase::kAction, "angle track selection"        , new TracksAngleAction()         , true);
  AddStep(StepBase::kCut   , "check vector is not empty"    , new EventHasTracksCut()         , true);
  //AddStep(StepBase::kCut   , "check hits is volume"         , new TracksInDefinedVolumeCut()  , true);

  SetBranchAlias(0,"trunk");
}

//**************************************************
bool GetNonStoppingTracksAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxdEdx& box = *static_cast<ToyBoxdEdx*>(&boxB); 
    
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(!pdAnaUtils::IsStoppingInFV(part))box.Tracks.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool TracksAngleAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxdEdx& box = *static_cast<ToyBoxdEdx*>(&boxB);   
  
  //keep only tracks with a particular angular distribution
  for(int itrk = 0; itrk < (int)box.Tracks.size(); itrk++){
    if((abs(180/TMath::Pi()*box.Tracks[itrk]->ThetaXZ)>60 && 
	abs(180/TMath::Pi()*box.Tracks[itrk]->ThetaXZ)<120)
       ||
       (abs(180/TMath::Pi()*box.Tracks[itrk]->ThetaYZ)>80 && 
	abs(180/TMath::Pi()*box.Tracks[itrk]->ThetaYZ)<100)){
      box.Tracks.erase(box.Tracks.begin()+itrk);
      itrk--;
    }
  }

  return true;
}

//**************************************************
bool EventHasTracksCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxdEdx& box = *static_cast<ToyBoxdEdx*>(&boxB);   

  //if there is at least one track to calibrate, pass
  if (box.Tracks.size()>0) return true;
  else return false;
}

//**************************************************
bool TracksInDefinedVolumeCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  bool result = false;

  //cast the box
  ToyBoxdEdx& box = *static_cast<ToyBoxdEdx*>(&boxB);   

  for(int itrk = 0; itrk < (int)box.Tracks.size(); itrk++){
    AnaParticlePD* part = box.Tracks[itrk];
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      if(part->Hits[2][ihit].Position.X() > -360 && part->Hits[2][ihit].Position.X() < 0 &&
	 abs(part->Hits[2][ihit].Position.Y()-300) < 20 &&
	 abs(part->Hits[2][ihit].Position.Z()-345) < 20){
	result = true;
	break;
      }
    }
    if(result)break;
  }
  
  return result;
}

//**************************************************
void dEdxTrackSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxdEdx])
    event.EventBoxes[EventBoxId::kEventBoxdEdx] = new EventBoxdEdx();

  boxUtils::FillAllTracksdEdx(event);
}

