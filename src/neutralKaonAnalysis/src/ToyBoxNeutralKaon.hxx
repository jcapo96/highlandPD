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
  void UpdateBestTrueVertexCandidateIndex(const int AccumLevel, const int Index);
  void UpdateBestReconVertexCandidateIndex(const int AccumLevel, const int Index);

public:

    /// Maximum accumulation level (not the real accum level, just a counter)
    int MaxAccumLevel;

    /// Simple counter for beam daughters
    int nBeamDaughters;

    /// Counter for all particles with valid start positions
    int nAllParticles;

    /// Flag indicating if K0 exists in truth as beam daughter
    bool hasK0InTruth;

    /// Vector of true vertex candidates
    std::vector<AnaTrueVertexPD*> trueVertexCandidates;

    /// Vector of reconstructed vertex candidates
    std::vector<AnaVertexPD*> reconVertexCandidates;

    /// Number of true vertex candidates
    int nTrueVertexCandidates;

    /// Number of reconstructed vertex candidates
    int nReconVertexCandidates;

    /// Index of best true vertex candidate
    int BestTrueVertexCandidateIndex;

    /// Index of best reconstructed vertex candidate
    int BestReconVertexCandidateIndex;
};

#endif
