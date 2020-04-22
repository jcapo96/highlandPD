#include "pandoraPreselection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pionAnalysisUtils.hxx"

//********************************************************************
pandoraPreselection::pandoraPreselection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void pandoraPreselection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")
  AddStep(StepBase::kAction, "find Pandora track",         new FindPandoraTrackAction());  
  AddStep(StepBase::kCut,    "candidate exists",           new CandidateExistsCut());
  
  // No preselection for the moment
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool FindPandoraTrackAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
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
    if(part->isPandora){
      box.MainTrack = part;
      break;
    }
  }
  return true;
}

//**************************************************
bool CandidateExistsCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  else return true;
}

//**************************************************
void pandoraPreselection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
}
