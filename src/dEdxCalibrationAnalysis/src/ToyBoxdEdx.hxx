#ifndef ToyBoxdEdx_h
#define ToyBoxdEdx_h

#include "ToyBoxPD.hxx"
#include "pdDataClasses.hxx"


class ToyBoxdEdx:public ToyBoxPD{
public :
  
  ToyBoxdEdx();
  virtual ~ToyBoxdEdx(){}

  /// This method should be implemented by the derived class. If so it does nothing here 
  virtual void Reset();

  /// Reset this base class
  virtual void ResetBase();

public:

  /// Vector of tracks
  std::vector<AnaParticlePD*> Tracks;
};

#endif
