#include "EventBoxKaon.hxx"
//#include "CutUtils.hxx"
#include "AnalysisUtils.hxx"
#include "pdDataClasses.hxx"

//********************************************************************
EventBoxKaon::EventBoxKaon():EventBoxB(){
//******************************************************************** 

}

//********************************************************************
EventBoxKaon::~EventBoxKaon(){
//********************************************************************

} 

//********************************************************************
void boxUtils::FillKaonCandidatesAndDaughters(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters], NMAXPARTICLES);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  int nCandidatesAndDaughters = 0;
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //skip beam particle
    if(part->isPandora)continue;
    //look for candidates
    if(part->DaughtersIDs.size() == 1){
      //check this particle hasn't been added before
      bool added = false;
      for(int j = 0; j < (int)NMAXPARTICLES; j++){
	if(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]] == part){
	  added = true;
	  break;
	}
      }
      //if it hasn't, add it
      if(!added){
	EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]++] = part;
	nCandidatesAndDaughters++;
      }
      //now go for the daughter
      AnaParticleB* dau = static_cast<AnaParticleB*>(part->Daughters[0]);
      if(!dau)continue;
      added = false;
      for(int j = 0; j < (int)NMAXPARTICLES; j++){
	if(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]] == dau){
	  added = true;
	  break;
	}
      }
      if(!added){
	EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]++] = dau;
	nCandidatesAndDaughters++;
      }
    }
  }
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxKaon::kCandidatesAndDaughters], 
                        EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters],
                        nCandidatesAndDaughters);
}

//********************************************************************
void boxUtils::FillTrueCandidatesAndDaughters(AnaEventB& event){
//********************************************************************
  
  (void)event;
  /*//  This method fills the EventBoxKaon with  the true candidate and its daughters for later use
  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  std::vector<AnaTrueParticle*> trueParts;
  int nCandidateAndDaughters = 0; 

  // Get all true parts in the event
  AnaTrueParticleB** parts = event.TrueParticles;
  Int_t nParts             = event.nTrueParticles;

  for(Int_t i = 0; i < nParts; i++){    
    AnaTrueParticlePD* part = static_cast<AnaTrueParticlePD*>(parts[i]);
    if (part->ProcessStart == AnaTrueParticle::primary){
      trueParts[nCandidateAndDaughters] = part;
      nCandidateAndDaughters++;
      const std::vector<Int_t>& daughtersIDs = part->Daughters;
      for (UInt_t j=0;j<daughtersIDs.size(); ++j){      
        AnaTrueParticleB* dau = anaUtils::GetTrueParticleByID(event, daughtersIDs[j]);
        if (!dau) continue;
        trueParts[nCandidateAndDaughters] = static_cast<AnaTrueParticle*>(dau);
        nCandidateAndDaughters++;
      }
      break;
    }
  }
      
  EventBox->nTrueObjectsInGroup[EventBoxKaon::kCandidateAndDaughters]=0;
  anaUtils::CreateArray(EventBox->TrueObjectsInGroup[EventBoxKaon::kCandidateAndDaughters], nCandidateAndDaughters);

  for (Int_t i=0;i<nCandidateAndDaughters; ++i){
    AnaTrueParticleB* part = trueParts[i];
    EventBox->TrueObjectsInGroup[EventBoxKaon::kCandidateAndDaughters][EventBox->nTrueObjectsInGroup[EventBoxKaon::kCandidateAndDaughters]++] = part;
  }

  anaUtils::ResizeArray(EventBox->TrueObjectsInGroup [EventBoxKaon::kCandidateAndDaughters], 
                        EventBox->nTrueObjectsInGroup[EventBoxKaon::kCandidateAndDaughters],
                        nCandidateAndDaughters);*/
}

//********************************************************************
void boxUtils::FillKaonXS(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxKaon::kKaonXS], 10*NMAXPARTICLES);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  int nCandidatesAndDaughters = 0;
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    EventBox->RecObjectsInGroup[EventBoxKaon::kKaonXS][EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS]++] = part;
    nCandidatesAndDaughters++;
  }
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxKaon::kKaonXS], 
                        EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS],
                        nCandidatesAndDaughters);
}
