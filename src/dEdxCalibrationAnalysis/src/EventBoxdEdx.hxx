#ifndef EventBoxdEdx_hxx
#define EventBoxdEdx_hxx

#include "BaseDataClasses.hxx"
//#include "EventBoxUtils.hxx"
#include "EventBoxPD.hxx"

/*
  The EventBox is used to store objects that are use many times when processing a single event, and 
  are not directly available in the event model. For example if we need a subsample of particles, this can be stored 
  in this box  
 */ 


class EventBoxdEdx:public EventBoxB{
 public :

  enum RecObjectGroupEnum{
    kAllTracksdEdx = 0,
  };
  
  enum TrueObjectGroupEnum{
    kAllTrueTracksdEdx = 0,
  };

  EventBoxdEdx();
  virtual ~EventBoxdEdx();
};


namespace boxUtils{

  void FillAllTracksdEdx(AnaEventB& event);
}

#endif
