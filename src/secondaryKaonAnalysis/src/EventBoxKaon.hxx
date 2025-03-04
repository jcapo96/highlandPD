#ifndef EventBoxKaon_hxx
#define EventBoxKaon_hxx

#include "BaseDataClasses.hxx"
//#include "EventBoxUtils.hxx"
#include "EventBoxPD.hxx"

/*
  The EventBox is used to store objects that are use many times when processing a single event, and 
  are not directly available in the event model. For example if we need a subsample of particles, this can be stored 
  in this box.
  This derived class stores relevant groups of particle for the secondary kaon analysis
 */ 


class EventBoxKaon:public EventBoxB{
 public :

  enum RecObjectGroupEnum{
    kCandidatesAndDaughters = 0,
    kKaonXS
  };
  
  enum TrueObjectGroupEnum{
    kTrueCandidatesAndDaughters = 0,
  };

  EventBoxKaon();
  virtual ~EventBoxKaon();
};


namespace boxUtils{

  //relevant particles for the secondary kaon analysis (candidates and daughters)
  void FillKaonCandidatesAndDaughters(AnaEventB& event, SubDetId::SubDetEnum det = SubDetId::kSubdet1);
  void FillTrueCandidatesAndDaughters(AnaEventB& event);

  //relevant particles for the XS analysis
  void FillKaonXS(AnaEventB& event, SubDetId::SubDetEnum det = SubDetId::kSubdet1);

  //relevant particles for the secondary proton selection for dEdx correction studies
  void FillProtonCandidates(AnaEventB& event, SubDetId::SubDetEnum det = SubDetId::kSubdet1);
}

#endif
