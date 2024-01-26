#include "secondaryKaonXSSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryKaonXSSelection::secondaryKaonXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************
  
}

//********************************************************************
void secondaryKaonXSSelection::DefineSteps(){
//********************************************************************

  //copy steps from pdBaseSelection
  AddStep(StepBase::kAction,   "find beam track",            new FindBeamTrackAction()   );// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "Pandora track exists",       new BeamTrackExistsCut()    );// in pdBaseAnalysis/src/pandoraPreselection  

  //select events using beam instrumentation
  //AddStep(StepBase::kCut,      "beam quality cut"    ,       new BeamQualityCut()    );
  AddStep(StepBase::kCut, "beam pdg filter", new BeamFilterForXSCut());

  //kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsForXSAction()   );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  AddStep(StepBase::kAction,   "get forced daughter",        new GetForcedDaughterAction()   );

  SetPreSelectionAccumLevel(-1);

  SetBranchAlias(0,"NO BRANCH :)");
}

//**************************************************
bool BeamFilterForXSCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut* beamCut1 = new BeamPDGCut(211);
  BeamPDGCut* beamCut2 = new BeamPDGCut(13);
  BeamPDGCut* beamCut3 = new BeamPDGCut(11);

  if(beamCut1->Apply(event,boxB) || beamCut2->Apply(event,boxB) || beamCut3->Apply(event,boxB))return true;
  else return false;
}

// //**************************************************
// bool BeamQualityCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
// //**************************************************
  
//   (void)event;
  
//   //cast the box
//   ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  
//   // Main track must exist
//   if(!box.MainTrack) return false;
  
//   // get beam particle
//   AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
//   AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
//   if(!beamPart)return false;
  
//   double data_mean_x  = -27.3;
//   double data_sigma_x = 4.76;
//   double data_mean_y  = 424.2;
//   double data_sigma_y = 4.75;
//   double data_mean_z  = 30.6;
//   double data_sigma_z = 1.20;

//   double mc_mean_x  = -28.6;
//   double mc_sigma_x = 4.15;
//   double mc_mean_y  = 422.6;
//   double mc_sigma_y = 3.94;
//   double mc_mean_z  = 29.7;
//   double mc_sigma_z = 0.54;

//   double mincos = 0.93;

//   double normalized_x = 0;
//   double normalized_y = 0;
//   double normalized_z = 0;
//   double cos = 0;

//   if(event.GetIsMC()){
//     normalized_x = (box.MainTrack->PositionStart[0]-mc_mean_x)/mc_sigma_x;
//     normalized_y = (box.MainTrack->PositionStart[1]-mc_mean_y)/mc_sigma_y;
//     normalized_z = (box.MainTrack->PositionStart[2]-mc_mean_z)/mc_sigma_z;
//   }
//   else{
//     normalized_x = (box.MainTrack->PositionStart[0]-data_mean_x)/data_sigma_x;
//     normalized_y = (box.MainTrack->PositionStart[1]-data_mean_y)/data_sigma_y;
//     normalized_z = (box.MainTrack->PositionStart[2]-data_mean_z)/data_sigma_z;
//   }

//   for(int i = 0; i < 3; i++)cos = cos + box.MainTrack->DirectionStart[i]*beamPart->DirectionEnd[i];

//   if(normalized_x < 3 && normalized_y < 3 && normalized_z < 3 && cos > mincos)return true;
//   else return false;
// }

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
    std::pair<double,int> pid = pdAnaUtils::Chi2PID(*part,321);
    double chi = pid.first/pid.second;
    double trunc = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,part->Hits[2]);
    if(chi>0 && chi<3.65 && trunc>2.4 && trunc<4.6)box.Candidates.push_back(part);
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
bool GetForcedDaughterAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  // loop over candidates
  for(int ican = 0; ican < (int)box.Candidates.size(); ican++){
    AnaParticlePD* candidate = box.Candidates[ican];
    // skip candidates with daughters
    if(!candidate->Daughters.empty())continue;
    double dis_min = 10000;
    int i_closest = -1;
    //loop over particles and look for the one starting closer to the end point of the candidate
    for(int ipart = 0; ipart < nParts; ipart++){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[ipart]);
      if(candidate->UniqueID == part->UniqueID)continue; //skip itself
      if(part->Type!=2)continue; //skip showers
      double dis = 0;
      for(int i = 0; i < 3; i++)dis += pow(part->PositionStart[i]-candidate->PositionEnd[i],2);
      dis = sqrt(dis);
      if(dis < dis_min){
	dis_min = dis;
	i_closest = ipart;
      }
    }
    if(i_closest == -1 || dis_min>5)continue;
    //associate the closest one as forced daughter
    AnaParticlePD* forced_dau = static_cast<AnaParticlePD*>(parts[i_closest]);
    candidate->forced_daughter = true;
    candidate->Daughters.push_back(forced_dau);
    
    //check true matching if possible
    AnaTrueParticlePD* truecandidate = static_cast<AnaTrueParticlePD*>(candidate->TrueObject);
    AnaTrueParticlePD* trueforced_dau = static_cast<AnaTrueParticlePD*>(forced_dau->TrueObject);
    if(!truecandidate || !trueforced_dau)continue;
    for(int i = 0; i < truecandidate->Daughters.size(); i++){
      if(truecandidate->Daughters[i] == trueforced_dau->ID){
	candidate->forced_daughter_matched = true;
	break;
      }
    }
    std::cout << "forced daughter found with distance " << dis_min << " " << candidate->forced_daughter_matched << std::endl;
  }
    
  return true;  
}


//**************************************************
void secondaryKaonXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonXS(event);
}

