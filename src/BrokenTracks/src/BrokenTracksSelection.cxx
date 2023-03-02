#include "BrokenTracksSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
BrokenTracksSelection::BrokenTracksSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxdEdx) {
//********************************************************************

}

//********************************************************************
void BrokenTracksSelection::DefineSteps(){
//********************************************************************

  AddStep(StepBase::kAction,   "find Pandora track",         new FindBeamTrackAction(),      false);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "candidate exists",           new BeamTrackExistsCut(),       false);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "beam pdg",                   new BeamMuonCut(),              false);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "beam quality",               new BeamTrackQualityCut(),      false);
  AddStep(StepBase::kCut,      "denominator",                new TrackReachesAPABorderCut(), false);
  AddStep(StepBase::kCut,      "almost numerator",           new TrackEndsInAPABorderCut(),  false);
  AddStep(StepBase::kCut,      "numerator",                  new TrackIsBrokenCut(),         false);

  SetBranchAlias(0,"trunk");
}

//**************************************************
bool BeamMuonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut beamCut1(13);
  BeamPDGCut beamCut2(211);

  if(beamCut1.Apply(event,boxB) || beamCut2.Apply(event,boxB))
    return true;
  else return false;
}

//**************************************************
bool BeamTrackQualityCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 
  return box.MainTrack->isBeamPart;
}

//**************************************************
bool TrackReachesAPABorderCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  return (box.MainTrack->PositionEnd[2] >= 220.);
}

//**************************************************
bool TrackEndsInAPABorderCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 
  return (box.MainTrack->PositionEnd[2] >= 220. && box.MainTrack->PositionEnd[2] <= 234.);
}

//**************************************************
bool TrackIsBrokenCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 

  bool IsBroken = false;
  double coscut = 0.8;
  double discut = 13.;

  AnaParticlePD* part = box.MainTrack;
  
  int ndau = part->Daughters.size();
  for(int idau = 0; idau < ndau; idau++){
    AnaParticlePD* dau = (AnaParticlePD*)part->Daughters[idau];
    double cos = 0, dis = 0;
    for(int i = 0; i < 3; i++){
      cos += part->DirectionEnd[i]*dau->DirectionStart[i];
      dis += pow(dau->PositionStart[i]-part->PositionEnd[i],2);
    }
    cos = abs(cos);
    dis = sqrt(dis);
    if(cos > coscut && dis < discut){
      IsBroken = true;
      break;
    }
  }

  return IsBroken;
}

//**************************************************
void BrokenTracksSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();
  
  boxUtils::FillCandidateAndDaughters(event);
}
