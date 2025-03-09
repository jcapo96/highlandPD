#ifndef ToyBoxNeutralKaon_h
#define ToyBoxNeutralKaon_h

#include "ToyBoxPD.hxx"
#include "pdDataClasses.hxx"


class ToyBoxNeutralKaon:public ToyBoxPD{
public :

  ToyBoxNeutralKaon();
  virtual ~ToyBoxNeutralKaon(){}

  /// This method should be implemented by the derived class. If so it does nothing here
  virtual void Reset();

  /// Reset this base class
  virtual void ResetBase();

  void UpdateBestCandidateIndex(const int AccumLevel, const int Index);

public:

    /// Vector of candidates
    // The vector of candidates is a vector of arrays of AnaParticlePD pointers
    // The first element of the array is the kaon candidate, the second is the pi1 candidate, and the third is the pi2 candidate
    std::vector<AnaParticlePD*> neutralKaonCandidates;
    int BestNeutralKaonCandidateIndex;
    int MaxAccumLevel; //not the real accum level, just a counter
};

#endif
