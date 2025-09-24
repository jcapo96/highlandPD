#include "ToyBoxNeutralKaon.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxNeutralKaon::ToyBoxNeutralKaon():ToyBoxPD(){
//********************************************************************

    MaxAccumLevel = -1;

    // Initialize new members for Preliminary K0 Selection
    nBeamDaughters = 0;
    nAllParticles = 0;
    hasK0InTruth = false;

    // Initialize vertex candidate vectors
    trueVertexCandidates.clear();
    reconVertexCandidates.clear();

    // Initialize vertex candidate counters and indices
    nTrueVertexCandidates = 0;
    nReconVertexCandidates = 0;
    BestTrueVertexCandidateIndex = -1;
    BestReconVertexCandidateIndex = -1;
}

//********************************************************************
void ToyBoxNeutralKaon::Reset(){
//********************************************************************

  if (!_ResetCheckDone){
    if( typeid(*this) !=  typeid(ToyBoxNeutralKaon)){
      std::cerr << "ERROR in ToyBoxKaon::Reset(). Either this mandatory method is not implemented "
                << "by the derived class '" << typeid(*this).name() << "' "
                << "or ToyBoxKaon::Reset() is called from the Reset method of the derived class. "
                << "Please correct any of these bugs. " << std::endl;

      exit(1);
    }
    _ResetCheckDone=true;
  }
}

//********************************************************************
void ToyBoxNeutralKaon::ResetBase(){
//********************************************************************

    ToyBoxPD::ResetBase();
    MaxAccumLevel = -1;

    // Reset new members for Preliminary K0 Selection
    nBeamDaughters = 0;
    nAllParticles = 0;
    hasK0InTruth = false;

    // Clear vertex candidate vectors
    trueVertexCandidates.clear();
    reconVertexCandidates.clear();

    // Reset vertex candidate counters and indices
    nTrueVertexCandidates = 0;
    nReconVertexCandidates = 0;
    BestTrueVertexCandidateIndex = -1;
    BestReconVertexCandidateIndex = -1;

    SoftReset();  // just reset internal stuff
}

//********************************************************************
void ToyBoxNeutralKaon::UpdateBestCandidateIndex(const int AccumLevel, const int Index){
//********************************************************************

  if(AccumLevel>MaxAccumLevel){
    MaxAccumLevel = AccumLevel;
  }
}

//********************************************************************
void ToyBoxNeutralKaon::UpdateBestTrueVertexCandidateIndex(const int AccumLevel, const int Index){
//********************************************************************

  if(AccumLevel>MaxAccumLevel){
    BestTrueVertexCandidateIndex = Index;
    MaxAccumLevel = AccumLevel;
  }
}

//********************************************************************
void ToyBoxNeutralKaon::UpdateBestReconVertexCandidateIndex(const int AccumLevel, const int Index){
//********************************************************************

  if(AccumLevel>MaxAccumLevel){
    BestReconVertexCandidateIndex = Index;
    MaxAccumLevel = AccumLevel;
  }
}
