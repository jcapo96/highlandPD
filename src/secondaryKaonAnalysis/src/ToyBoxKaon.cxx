#include "ToyBoxKaon.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxKaon::ToyBoxKaon():ToyBoxPD(){
//******************************************************************** 
  
  Candidates.clear();
} 

//********************************************************************
void ToyBoxKaon::Reset(){
//******************************************************************** 

  if (!_ResetCheckDone){
    if( typeid(*this) !=  typeid(ToyBoxKaon)){
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
void ToyBoxKaon::ResetBase(){
//******************************************************************** 

  ToyBoxPD::ResetBase();
  Candidates.clear();
  SoftReset();  // just reset internal stuff
}  
