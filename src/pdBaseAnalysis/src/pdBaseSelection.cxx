#include "pdBaseSelection.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
pdBaseSelection::pdBaseSelection(bool forceBreak):SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************  
}

//********************************************************************
void pdBaseSelection::DefineSteps(){
//********************************************************************

  AddStep(StepBase::kCut,    "event quality",     new EventQualityCut(),           true);
  AddStep(StepBase::kAction, "find beam track",   new FindBeamTrackAction());  
  AddStep(StepBase::kCut,    "beam trakc exists", new BeamTrackExistsCut());

  SetBranchAlias(0,"dummy");
  
  // No preselection for the moment
  SetPreSelectionAccumLevel(-1);

}

//**************************************************
bool EventQualityCut::Apply(AnaEventC& eventC, ToyBoxB& box) const{
//**************************************************

  (void)box;

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  if(event.GetIsMC())  return true;              // This is MC, ignore DQ. 

  if (enableDQCut) {
    if(!event.DataQuality->GoodDaq ) return false;  // Bad Detector Data quality
  }

  if (enableBeamQualityCut) {
    if(!event.Beam->GoodSpill ) return false;      // Bad Spill
  }

  return true;
}

//**************************************************
bool FindBeamTrackAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 

  // This action fills box.MainTrack 
  box.MainTrack = NULL;
  
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //select the one tagged as beam particle by Pandora
    if(part->isPandora){
      box.MainTrack = part;
      break;
    }
  }
  return true;
}

//**************************************************
bool BeamTrackExistsCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  else return true;
}

//**************************************************
bool BeamPDGCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  bool PDG = false;

  //Get the beam from the event and the beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  if(event.GetIsMC()){
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    if(!beamPart)return false;
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(beamPart->TrueObject);
    if(!truePart)return false;
    if(abs(_equal_cut)==abs(truePart->PDG))PDG = true;
  }
  else{
    if(std::find(beam->PDGs.begin(),beam->PDGs.end(),_equal_cut) != beam->PDGs.end())PDG = true;
  }

  return PDG;
}


//**************************************************
void pdBaseSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();
  
  boxUtils::FillCandidateAndDaughters(event);
}

//**************************************************
bool BeamQualityCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  //get BeamQualityCuts depending on MC/data and nominal momentum
  double mean_x, mean_y, mean_z, sigma_x, sigma_y, sigma_z, min_cos;
  pdAnaUtils::GetBeamQualityCuts(static_cast<AnaEventPD*>(&event),
				 mean_x,mean_y,mean_z,
				 sigma_x,sigma_y,sigma_z,
				 min_cos);

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  
  // Main track must exist
  if(!box.MainTrack) return false;
  
  // get beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if(!beamPart)return false;

  //computations
  double normalized_x = 0;
  double normalized_y = 0;
  double normalized_z = 0;
  double cos = 0;

  normalized_x = abs(box.MainTrack->PositionStart[0]-mean_x)/sigma_x;
  normalized_y = abs(box.MainTrack->PositionStart[1]-mean_y)/sigma_y;
  normalized_z = abs(box.MainTrack->PositionStart[2]-mean_z)/sigma_z;
  for(int i = 0; i < 3; i++)cos = cos + box.MainTrack->DirectionStart[i]*beamPart->DirectionEnd[i];

  //check cut
  if(normalized_x < 3 && normalized_y < 3 && normalized_z < 3 && cos > min_cos)return true;
  else return false;
}
