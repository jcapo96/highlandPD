#include "TPad.h"
#include "TMultiGraph.h"
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
  AddStep(StepBase::kAction,   "find beam track",            new FindBeamTrackAction(), true   );// in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,      "Pandora track exists",       new BeamTrackExistsCut(), true    );// in pdBaseAnalysis/src/pandoraPreselection  

  //select events using beam instrumentation
  AddStep(StepBase::kCut,      "beam quality cut"    ,       new BeamQualityCut()    );
  //AddStep(StepBase::kCut,      "beam quality cut"    ,       new BQCForXS()    );
  //AddStep(StepBase::kCut, "beam pdg filter", new BeamFilterForXSCut(), true);

  //split the selection in branches, pions/protons/kaonsfor each possible candidate
  AddSplit(3);

  //pions, protons and kaons cuts
  AddStep(0,StepBase::kCut, "beam pdg filter", new BeamFilterForXSCut(), true);
  AddStep(1,StepBase::kCut, "beam pdg filter", new BeamPDGCut(2212), true);
  AddStep(2,StepBase::kCut, "beam pdg filter", new BeamPDGCut(321), true);

  //add cuts to branches
  for(int i = 0; i < 3; i++){
    AddStep(i,StepBase::kAction,   "get a vector of kaon candidates",      new GetKaonsForXSAction(),true   );
    AddStep(i,StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
    AddStep(i,StepBase::kAction,   "get forced daughter",        new GetForcedDaughterAction()   );
    //AddStep(i,StepBase::kCut,      "dau mom muon cut",           new DauMomMuonCut(), true);
    AddStep(i,StepBase::kCut,      "all cut",           new AllCut(), true);
    //AddStep(i,StepBase::kAction,   "test",           new DrawingTestAction(), true);
  }

  //set alias for branches
  SetBranchAlias(0,"beam pions"  ,0);
  SetBranchAlias(1,"beam protons",1);
  SetBranchAlias(2,"beam kaons"  ,2);

  SetPreSelectionAccumLevel(-1);

  //  SetBranchAlias(0,"NO BRANCH :)");
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

//**************************************************
bool BQCForXS::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // get particles
  AnaParticlePD* part = box.MainTrack;
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticlePD* beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if(!part || !beampart)return false;

  //get anaeventinfo
  AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(static_cast<AnaEventB*>(&event)->EventInfo);
  double mom = info->NominalBeamMom;

  double data_6_mean_x = 2.51;
  double data_6_mean_y = 1.90;
  double data_6_mean_z = 1.02;
  double data_7_mean_x = 3.55;   
  double data_7_mean_y = 1.68;   
  double data_7_mean_z = 1.76;   
  double data_6_sigma_x = 0.26;  
  double data_6_sigma_y = 1.1;   
  double data_6_sigma_z = 1.08; 
  double data_7_sigma_x = 0.29;  
  double data_7_sigma_y = 0.91;  
  double data_7_sigma_z = 1.51;                        
  double mc_6_mean_x = 0.99;    
  double mc_6_mean_y = 0.;       
  double mc_6_mean_z = -0.01;
  double mc_7_mean_x = -1.50; 
  double mc_7_mean_y = 0.71;     
  double mc_7_mean_z = 0.;    
  double mc_6_sigma_x = 0.086;   
  double mc_6_sigma_y = 0.11;    
  double mc_6_sigma_z = 0.20;  
  double mc_7_sigma_x = 0.20; 
  double mc_7_sigma_y = 0.22;    
  double mc_7_sigma_z = 0.19;

  double meanx,sigmax,meany,sigmay,meanz,sigmaz;
  if(event.GetIsMC()){
    if(mom == 6){
      meanx = mc_6_mean_x;
      meany = mc_6_mean_y;
      meanz = mc_6_mean_z;
      sigmax = mc_6_sigma_x;
      sigmay = mc_6_sigma_y;
      sigmaz = mc_6_sigma_z;
    }
    else if(mom == 7){
      meanx = mc_7_mean_x;
      meany = mc_7_mean_y;
      meanz = mc_7_mean_z;
      sigmax = mc_7_sigma_x;
      sigmay = mc_7_sigma_y;
      sigmaz = mc_7_sigma_z;
    }
  }
  else{
    if(mom == 6){
      meanx = data_6_mean_x;
      meany = data_6_mean_y;
      meanz = data_6_mean_z;
      sigmax = data_6_sigma_x;
      sigmay = data_6_sigma_y;
      sigmaz = data_6_sigma_z;
    }
    else if(mom == 7){
      meanx = data_7_mean_x;
      meany = data_7_mean_y;
      meanz = data_7_mean_z;
      sigmax = data_7_sigma_x;
      sigmay = data_7_sigma_y;
      sigmaz = data_7_sigma_z;
    }
  }

  double deltaX = abs(part->PositionStart[0]-beampart->PositionEnd[0]-meanx)/sigmax;
  double deltaY = abs(part->PositionStart[1]-beampart->PositionEnd[1]-meany)/sigmay;
  double deltaZ = abs(part->PositionStart[2]-beampart->PositionEnd[2]-meanz)/sigmaz;
  double coss = part->DirectionStart[0]*beampart->DirectionEnd[0]+part->DirectionStart[1]*beampart->DirectionEnd[1]+part->DirectionStart[2]*beampart->DirectionEnd[2];

  if(deltaX<3 && deltaY<3 && deltaZ<3 && coss>0.93)return true;
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
    if((int)part->Daughters.size()>1)continue;
    std::pair<double,int> pid = pdAnaUtils::Chi2PID(*part,321);
    if(pid.first<0)continue;
    double chi_kaon = pid.first/pid.second;
    pid = pdAnaUtils::Chi2PID(*part,2212);
    if(pid.first<0)continue;
    double chi_prot = pid.first/pid.second;
    if(chi_kaon>0 && chi_kaon<20 && chi_prot>10)
      box.Candidates.push_back(part);
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
      if(part->isPandora)continue; //skip beam particle
      if(part->Type!=2)continue; //skip showers maybe not needed?
      double dis = 0;
      for(int i = 0; i < 3; i++)dis += pow(part->PositionStart[i]-candidate->PositionEnd[i],2);
      dis = sqrt(dis);
      if(dis < dis_min){
	dis_min = dis;
	i_closest = ipart;
      }
    }
    if(i_closest == -1 || dis_min>8)continue;
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
  }
    
  return true;  
}

//**************************************************
bool DauMomMuonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //loop over candidates
  for(std::vector<AnaParticlePD*>::iterator it = box.Candidates.begin(); it != box.Candidates.end(); ){
    //check it has a single daughter
    if((*it)->Daughters.size()!=1){
      it = box.Candidates.erase(it);
      continue;
    }
    //check mom condition
    AnaParticlePD* dau = static_cast<AnaParticlePD*>((*it)->Daughters[0]);
    double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
    if(abs(mom-0.23)<0.03)
      it++;
    else
      it = box.Candidates.erase(it);
  }

  return !box.Candidates.empty();
}

//**************************************************
bool AllCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //loop over candidates
  for(std::vector<AnaParticlePD*>::iterator it = box.Candidates.begin(); it != box.Candidates.end(); ){
    //check it has a single daughter
    if((*it)->Daughters.size()!=1){
      it = box.Candidates.erase(it);
      continue;
    }
    //check conditions
    AnaParticlePD* dau = static_cast<AnaParticlePD*>((*it)->Daughters[0]);
    double candidates_truncated_dedx = (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,(*it)->Hits[2]);
    double candidates_dau_truncated_dedx = (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dau->Hits[2]);
    double truncated_dedx = (candidates_truncated_dedx+candidates_dau_truncated_dedx)/2;
    double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
    double dis = pdAnaUtils::ComputeDistanceMotherDaughter((*it),dau);
    bool forced = (*it)->forced_daughter;
    if(abs(mom-0.215)<0.045 && truncated_dedx<4.2 && ((forced && dis<2.5) || (!forced && dis<4)))
      it++;
    else
      it = box.Candidates.erase(it);
  }

  return !box.Candidates.empty();
}

//**************************************************
bool DrawingTestAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  //get anaeventinfo
  AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(static_cast<AnaEventB*>(&event)->EventInfo);
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // loop over candidates
  for(int ican = 0; ican < (int)box.Candidates.size(); ican++){
    AnaParticlePD* candidate = box.Candidates[ican];
    //get daughter
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(candidate->Daughters[0]);
    //skip if not forced
    if(!candidate->forced_daughter)continue;

    //title
    std::stringstream ssevent, ssmatch;
    ssevent << info->Event;
    if(candidate->forced_daughter_matched)
      ssmatch << 1;
    else
      ssmatch << 0;
    std::string title = ssevent.str()+"_"+ssmatch.str();
    
    TMultiGraph* mgg = new TMultiGraph();
    TGraph* tg1 = new TGraph();
    TGraph* tg2 = new TGraph();
    int counter = 0;
    for(int i = 0; i < candidate->Hits[2].size(); i++){
      tg1->SetPoint(counter,candidate->Hits[2][i].Position.Z(),candidate->Hits[2][i].Position.Y());
      counter++;
    }
    counter = 0;
    for(int i = 0; i < dau->Hits[2].size(); i++){
      tg2->SetPoint(counter,dau->Hits[2][i].Position.Z(),dau->Hits[2][i].Position.Y());
      counter++;
    }
    tg1->SetLineColor(1);
    tg1->SetMarkerColor(1);
    tg2->SetLineColor(2);
    tg2->SetMarkerColor(2);
    mgg->Add(tg1,"pl*");
    mgg->Add(tg2,"pl*");
    mgg->SetTitle(title.c_str());
    mgg->SetName(title.c_str());
    mgg->GetXaxis()->SetTitle("Z [cm]");
    mgg->GetYaxis()->SetTitle("Y [cm]");
    mgg->GetXaxis()->CenterTitle();
    mgg->GetYaxis()->CenterTitle();
    mgg->Draw("apl");
    gPad->Update();gPad->Print(("plots/"+title+".jpeg").c_str());
    delete tg1;
    delete tg2;
    delete mgg;
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

