#include "EventBoxdEdx.hxx"
//#include "CutUtils.hxx"
#include "AnalysisUtils.hxx"
#include "pdDataClasses.hxx"

//********************************************************************
EventBoxdEdx::EventBoxdEdx():EventBoxB(){
//******************************************************************** 

}

//********************************************************************
EventBoxdEdx::~EventBoxdEdx(){
//********************************************************************

} 

//********************************************************************
void boxUtils::FillAllTracks(AnaEventB& event){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxdEdx];

  EventBox->nRecObjectsInGroup[EventBoxdEdx::kAllTracks]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxdEdx::kAllTracks],30);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    EventBox->RecObjectsInGroup[EventBoxdEdx::kAllTracks][EventBox->nRecObjectsInGroup[EventBoxdEdx::kAllTracks]++] = part;
  }
  
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxdEdx::kAllTracks], 
                        EventBox->nRecObjectsInGroup[EventBoxdEdx::kAllTracks],
                        nParts);
}
