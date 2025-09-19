#include "ToyBoxNeutralKaon.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxNeutralKaon::ToyBoxNeutralKaon():ToyBoxPD(){
//********************************************************************

    neutralKaonCandidates.clear();
    BestNeutralKaonCandidateIndex = -1;
    MaxAccumLevel = -1;

    // Initialize new members for Preliminary K0 Selection
    nBeamDaughters = 0;
    hasK0InTruth = false;
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
    neutralKaonCandidates.clear();
    BestNeutralKaonCandidateIndex = -1;
    MaxAccumLevel = -1;

    // Reset new members for Preliminary K0 Selection
    nBeamDaughters = 0;
    hasK0InTruth = false;

    SoftReset();  // just reset internal stuff
}

//********************************************************************
void ToyBoxNeutralKaon::UpdateBestCandidateIndex(const int AccumLevel, const int Index){
//********************************************************************

  if(AccumLevel>MaxAccumLevel){
    BestNeutralKaonCandidateIndex =  Index;
    MaxAccumLevel = AccumLevel;
  }
}
