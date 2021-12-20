#include "secondaryKaonFromBeamSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

//********************************************************************
secondaryKaonFromBeamSelection::secondaryKaonFromBeamSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void secondaryKaonFromBeamSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselecion
  AddStep(StepBase::kAction,   "find Pandora track",         new FindPandoraTrackAction());// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "candidate exists",           new CandidateExistsCut()    );// in pdBaseAnalysis/src/pandoraPreselection  
  //select events using beam instrumentation
  AddStep(StepBase::kCut,      "beam pdg filter",            new BeamFilterCut()         );
  //secondary kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsFromBeamAction());
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),      true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(),         true);
    AddStep(i, StepBase::kCut, "kaon daughter CNN cut",      new MuonCNNCut(),          true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(),     true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(),          true);
    AddStep(i, StepBase::kCut, "kaon-muon angle cut",        new KaonMuonAngleCut(),    true);
    AddStep(i, StepBase::kCut, "kaon-muon distance cut",     new KaonMuonDistanceCut(), true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    std::stringstream ssi;
    ssi << i;
    SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  }  

  // This number tells the selection to stop when a given accum level (accumulated cut level) is not passed. 
  // For example, 2 as argument would mean that selection is stopped for the current toy experiment when the 
  // third cut (numbering starts at 0) is not passed. This is a way of saving time.
  // -1 means no preselection
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool BeamFilterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)boxB;
  
  //get the beam
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  //for MC
  if(event.GetIsMC()){
    if(beam->BeamParticle){
      if(beam->BeamParticle->TrueObject){
        //cast the beam particle
        AnaTrueParticlePD* beamTruePart = static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);
        if(abs(beamTruePart->PDG)==13 || abs(beamTruePart->PDG)==211 || abs(beamTruePart->PDG)==321 || abs(beamTruePart->PDG)==2212)return true;
      }
    }
  }
  //for data
  else{
    bool hadron = false;
    for(size_t i = 0; i < beam->PDGs.size(); i++){
      if(beam->PDGs[i] == 211 || beam->PDGs[i] == 321 || beam->PDGs[i] == 2212){
        hadron = true;
        break;
      }
    }
    if(hadron)return true;
  }
  
  return false;
}

//**************************************************
bool GetKaonsFromBeamAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB); 
    
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part->isPandora)continue; //skip beam particle
    if(part->Daughters.size() == 1 && part->ParentID != -1 && part->BeamOrigin)box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
void secondaryKaonFromBeamSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}

