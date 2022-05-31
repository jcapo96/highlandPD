#ifndef EventBoxKaon_hxx
#define EventBoxKaon_hxx

#include "BaseDataClasses.hxx"
//#include "EventBoxUtils.hxx"
#include "EventBoxPD.hxx"

/*
  The EventBox is used to store objects that are use many times when processing a single event, and 
  are not directly available in the event model. For example if we need a subsample of particles, this can be stored 
  in this box  
 */ 


class EventBoxKaon:public EventBoxB{
 public :

  enum RecObjectGroupEnum{
    kCandidatesAndDaughters = 0,
  };
  
  enum TrueObjectGroupEnum{
    kTrueCandidatesAndDaughters = 0,
  };

  EventBoxKaon();
  virtual ~EventBoxKaon();
};


namespace boxUtils{

  void FillKaonCandidatesAndDaughters(AnaEventB& event, SubDetId::SubDetEnum det = SubDetId::kSubdet1);
  void FillTrueCandidatesAndDaughters(AnaEventB& event);
}

#endif
