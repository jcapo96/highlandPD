#include "secondaryKaonXSSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryKaonXSSelection::secondaryKaonXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************
  
  //Read parameters before base initialization because DefineSteps needs parameters to be defined
  _BeamPDGFilter = ND::params().GetParameterI("secondaryKaonAnalysis.XSSelection.BeamPDGFilter");

}

//********************************************************************
void secondaryKaonXSSelection::DefineSteps(){
//********************************************************************

  //copy steps from pdBaseSelection
  AddStep(StepBase::kAction,   "find beam track",            new FindBeamTrackAction()   );// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "Pandora track exists",       new BeamTrackExistsCut()    );// in pdBaseAnalysis/src/pandoraPreselection  

  //select events using beam instrumentation
  //AddStep(StepBase::kCut,      "beam quality cut"    ,       new BeamQualityCut()    );
  if(_BeamPDGFilter==211 )AddStep(StepBase::kCut, "beam pdg filter", new BeamFilterForXSCut());
  if(_BeamPDGFilter==2212)AddStep(StepBase::kCut, "beam pdg filter", new BeamPDGCut(2212)    );
  if(_BeamPDGFilter==321 )AddStep(StepBase::kCut, "beam pdg filter", new BeamPDGCut(321)     );

  //kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsForXSAction()   );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),             true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(0.64,5.84),       true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(0.210,0.237), true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(0.62,1.0),         true);
    AddStep(i, StepBase::kCut, "kaon-muon distance cut",     new KaonMuonDistanceCut(0,16),    true);
    AddStep(i, StepBase::kCut, "kaon chi2 cut",              new KaonChi2Cut(),                true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    std::stringstream ssi;
    ssi << i;
    SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  }

  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool BeamFilterForXSCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut* beamCut1 = new BeamPDGCut(211);
  BeamPDGCut* beamCut2 = new BeamPDGCut(13);

  if(beamCut1->Apply(event,boxB) || beamCut2->Apply(event,boxB))return true;
  else return false;
}

//**************************************************
bool BeamQualityCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;
  
  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  
  // Main track must exist
  if(!box.MainTrack) return false;
  
  // get beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if(!beamPart)return false;
  
  double data_mean_x  = -27.3;
  double data_sigma_x = 4.76;
  double data_mean_y  = 424.2;
  double data_sigma_y = 4.75;
  double data_mean_z  = 30.6;
  double data_sigma_z = 1.20;

  double mc_mean_x  = -28.6;
  double mc_sigma_x = 4.15;
  double mc_mean_y  = 422.6;
  double mc_sigma_y = 3.94;
  double mc_mean_z  = 29.7;
  double mc_sigma_z = 0.54;

  double mincos = 0.93;

  double normalized_x = 0;
  double normalized_y = 0;
  double normalized_z = 0;
  double cos = 0;

  if(event.GetIsMC()){
    normalized_x = (box.MainTrack->PositionStart[0]-mc_mean_x)/mc_sigma_x;
    normalized_y = (box.MainTrack->PositionStart[1]-mc_mean_y)/mc_sigma_y;
    normalized_z = (box.MainTrack->PositionStart[2]-mc_mean_z)/mc_sigma_z;
  }
  else{
    normalized_x = (box.MainTrack->PositionStart[0]-data_mean_x)/data_sigma_x;
    normalized_y = (box.MainTrack->PositionStart[1]-data_mean_y)/data_sigma_y;
    normalized_z = (box.MainTrack->PositionStart[2]-data_mean_z)/data_sigma_z;
  }

  for(int i = 0; i < 3; i++)cos = cos + box.MainTrack->DirectionStart[i]*beamPart->DirectionEnd[i];

  if(normalized_x < 3 && normalized_y < 3 && normalized_z < 3 && cos > mincos)return true;
  else return false;
}

//**************************************************
bool GetKaonsForXSAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // Get the beam track
  AnaParticlePD* beamTrack = box.MainTrack;
    
  //look over daughters and look for candidates
  for(int i = 0; i < (int)beamTrack->Daughters.size(); i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(beamTrack->Daughters[i]);
    if(part->DaughtersIDs.size() == 1)box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool KaonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its chi2
  AnaParticlePD* kaon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]);
  std::pair<double,int> result = pdAnaUtils::Chi2PID(*kaon,321);
  if(result.first/result.second < 40)return true;
  else return false;
}

//**************************************************
void secondaryKaonXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}

