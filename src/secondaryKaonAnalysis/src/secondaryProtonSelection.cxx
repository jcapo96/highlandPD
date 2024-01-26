#include "secondaryProtonSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "secondaryKaonAnalysis.hxx"
#include "EventBoxKaon.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryProtonSelection::secondaryProtonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void secondaryProtonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselection
  AddStep(StepBase::kAction, "FIND PANDORA TRACK"   , new FindBeamTrackAction(), true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut   , "BEAM TRACK EXISTS"    , new BeamTrackExistsCut() , true);// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kAction, "GET VECTOR OF PROTONS", new GetProtonsAction()   , true);
  AddStep(StepBase::kCut   , "CANDIDATE EXISTS"     , new EventHasKaonCut()    , true);
  AddStep(StepBase::kCut   , "SECONDARY PROTONS CUT", new SecondaryProtonsCut(0,25), true);
  // //next cuts have to be applied to each branch
  // AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  // for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++)
  //   AddStep(i, StepBase::kCut, "Proton chi2 cut"    , new ProtonChi2Cut(0,25)  , true);

  // // This number tells the selection to stop when a given accum level (accumulated cut level) is not passed. 
  // // For example, 2 as argument would mean that selection is stopped for the current toy experiment when the 
  // // third cut (numbering starts at 0) is not passed. This is a way of saving time.
  // // -1 means no preselection
  // //Set the branch aliases to the different branches 
  // for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
  //   std::stringstream ssi;
  //   ssi << i;
  //   SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  // }
  SetBranchAlias(0,"trunk");
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool GetProtonsAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
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
    if(part->DaughtersIDs.size() == 0 && part->ParentID != -1 && part->Type==2)
      box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool SecondaryProtonsCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB); 
    
  // loop over candidates
  int npass = 0;
  int ipart = 0;
  while(ipart < (int)box.Candidates.size()){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box.Candidates[ipart]);
    double chi2 = part->Chi2Proton;
    double ndf  = part->Chi2ndf;
    if(chi2/ndf > _lower_cut && chi2/ndf < _upper_cut && chi2 > 0 && chi2>0){
      npass++;
      ipart++;
    }
    else 
      box.Candidates.erase(box.Candidates.begin()+ipart);
  }
    
  if(npass>0)return true;  
  else return false;
}


//**************************************************
void secondaryProtonSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************
  
  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 
  
  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();
  
  boxUtils::FillProtonCandidates(event);
}

