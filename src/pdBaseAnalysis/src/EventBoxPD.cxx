#include "EventBoxPD.hxx"
#include "CutUtils.hxx"
#include "DataClasses.hxx"

//********************************************************************
EventBoxPD::EventBoxPD():EventBoxB(){
//******************************************************************** 

}

//********************************************************************
EventBoxPD::~EventBoxPD(){
//********************************************************************

} 

//********************************************************************
void boxUtils::FillLongTracks(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxPD];

  // TODO. We need a better method to get the tracks (returning a variable size array). Otherwise we have to guess how many tracks we have
  AnaTrackB* selTracks[NMAXPARTICLES*10];
  int nLong =anaUtils::GetAllTracksUsingDet(event,det, selTracks);

  nLong = std::min(nLong, (Int_t)NMAXPARTICLES);

  EventBox->nRecObjectsInGroup[EventBoxPD::kAllTracks]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxPD::kAllTracks], nLong);
    
  EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracks]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxPD::kLongTracks], nLong);

  EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracksInFV]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxPD::kLongTracksInFV], nLong);


  //loop over fgd tracks
  for (Int_t i=0;i<nLong; ++i){
    AnaParticleB* track = selTracks[i];
    EventBox->RecObjectsInGroup[EventBoxPD::kAllTracks][EventBox->nRecObjectsInGroup[EventBoxPD::kAllTracks]++] = track;
    if (track->NHits>20){      
      EventBox->RecObjectsInGroup[EventBoxPD::kLongTracks][EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracks]++] = track;
      if (cutUtils::FiducialCut(*track, det))
        EventBox->RecObjectsInGroup[EventBoxPD::kLongTracksInFV][EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracksInFV]++] = track;
    }
  }

  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxPD::kAllTracks], 
                        EventBox->nRecObjectsInGroup[EventBoxPD::kAllTracks],
                        nLong);
  
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxPD::kLongTracks], 
                        EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracks],
                        nLong);

  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxPD::kLongTracksInFV], 
                        EventBox->nRecObjectsInGroup[EventBoxPD::kLongTracksInFV],
                        nLong);
}


//********************************************************************
void boxUtils::FillCandidateAndDaughters(AnaEventB& event){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxPD];

  // TODO. We need a better method to get the tracks (returning a variable size array). Otherwise we have to guess how many tracks we have
  AnaParticle* selTracks[NMAXPARTICLES*10];
  int nCandidateAndDaughters = 0; 

  // Get all reconstructed parts in the event
  AnaParticleB** parts = event.Particles;
  Int_t nParts         = event.nParticles;

  for (Int_t i=0;i<nParts; ++i){    
    AnaParticle* part = static_cast<AnaParticle*>(parts[i]);
    if (part->Charge==-8888){
      selTracks[nCandidateAndDaughters] = part;
      nCandidateAndDaughters++;
      const std::vector<AnaRecObjectC*>& daughters = part->Daughters;
      for (UInt_t j=0;j<part->Daughters.size(); ++j){      
        selTracks[nCandidateAndDaughters] = static_cast<AnaParticle*>(daughters[j]);
        nCandidateAndDaughters++;
      }
      break;
    }
  }
      
  EventBox->nRecObjectsInGroup[EventBoxPD::kCandidateAndDaughters]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxPD::kCandidateAndDaughters], nCandidateAndDaughters);

  for (Int_t i=0;i<nCandidateAndDaughters; ++i){
    AnaParticleB* track = selTracks[i];
    EventBox->RecObjectsInGroup[EventBoxPD::kCandidateAndDaughters][EventBox->nRecObjectsInGroup[EventBoxPD::kCandidateAndDaughters]++] = track;
  }

  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxPD::kCandidateAndDaughters], 
                        EventBox->nRecObjectsInGroup[EventBoxPD::kCandidateAndDaughters],
                        nCandidateAndDaughters);

  
}
