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

    // Initialize neutral particle candidate vector
    neutralParticleCandidates.clear();

    // Initialize neutral particle candidate counter and index
    nNeutralParticleCandidates = 0;
    BestNeutralParticleCandidateIndex = -1;
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

    // Clear neutral particle candidate vector
    neutralParticleCandidates.clear();

    // Reset neutral particle candidate counter and index
    nNeutralParticleCandidates = 0;
    BestNeutralParticleCandidateIndex = -1;

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
void ToyBoxNeutralKaon::UpdateBestNeutralParticleCandidateIndex(const int AccumLevel, const int Index){
//********************************************************************

  if(AccumLevel>MaxAccumLevel){
    BestNeutralParticleCandidateIndex = Index;
    MaxAccumLevel = AccumLevel;
  }
}
