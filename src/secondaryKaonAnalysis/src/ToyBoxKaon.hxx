#ifndef ToyBoxKaon_h
#define ToyBoxKaon_h

#include "ToyBoxPD.hxx"
#include "pdDataClasses.hxx"


class ToyBoxKaon:public ToyBoxPD{
public :
  
  ToyBoxKaon();
  virtual ~ToyBoxKaon(){}

  /// This method should be implemented by the derived class. If so it does nothing here 
  virtual void Reset();

  /// Reset this base class
  virtual void ResetBase();

  void UpdateBestCandidateIndex(const int AccumLevel, const int Index);
  
public:

  /// Vector of candidates
  std::vector<AnaParticlePD*> Candidates;
  int BestCandidateIndex;
  int MaxAccumLevel; //not the real accum level, just a counter
};

#endif
