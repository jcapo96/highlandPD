#ifndef EventBoxPD_hxx
#define EventBoxPD_hxx

#include "BaseDataClasses.hxx"
#include "EventBoxUtils.hxx"

/*
  The EventBox is used to store objects that are use many times when processing a single event, and 
  are not directly available in the event model. For example if we need a subsample of particles, this can be stored 
  in this box  
 */ 


class EventBoxPD:public EventBoxB{
 public :

  enum RecObjectGroupEnum{
    kTracksUnassigned=0,
    kCandidateAndDaughters,  // Beam particle and its daughters tracks
    kAllTracks,       // All tracks
    kLongTracks,      // Tracks with more than 20 hits
    kLongTracksInFV   // Tracks with more than 20 hits in FV    
  };
  
  enum TrueObjectGroupEnum{
    kTrueParticlesUnassigned=0,
    kTrueParticlesChargedInTPCInBunch,
  };

  EventBoxPD();
  virtual ~EventBoxPD();
};


namespace boxUtils{

  /// Fill in the EventBox several arrays of tracks with Subdet2
  void FillLongTracks(AnaEventB& event, SubDetId::SubDetEnum det = SubDetId::kSubdet1);
  void FillCandidateAndDaughters(AnaEventB& event);
  void FillTrueCandidateAndDaughters(AnaEventB& event);
}

#endif
