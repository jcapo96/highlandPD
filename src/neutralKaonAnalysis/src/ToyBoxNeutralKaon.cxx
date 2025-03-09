#include "ToyBoxNeutralKaon.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxNeutralKaon::ToyBoxNeutralKaon():ToyBoxPD(){
//********************************************************************

    neutralKaonCandidates.clear();
    BestNeutralKaonCandidateIndex = -1;
    MaxAccumLevel = -1;
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
    std::cout << "Size of Candidates before clear: " << neutralKaonCandidates.size() << std::endl;
    neutralKaonCandidates.clear();
    BestNeutralKaonCandidateIndex = -1;
    MaxAccumLevel = -1;
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
