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
  void UpdateBestNeutralParticleCandidateIndex(const int AccumLevel, const int Index);

public:

    /// Maximum accumulation level (not the real accum level, just a counter)
    int MaxAccumLevel;

    /// Simple counter for beam daughters
    int nBeamDaughters;

    /// Counter for all particles with valid start positions
    int nAllParticles;

    /// Flag indicating if K0 exists in truth as beam daughter
    bool hasK0InTruth;

    /// Vector of neutral particle candidates
    std::vector<AnaNeutralParticlePD*> neutralParticleCandidates;

    /// Number of neutral particle candidates
    int nNeutralParticleCandidates;

    /// Index of best neutral particle candidate
    int BestNeutralParticleCandidateIndex;
};

#endif
