#include "ToyBoxdEdx.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxdEdx::ToyBoxdEdx():ToyBoxPD(){
//******************************************************************** 
  
  Tracks.clear();
} 

//********************************************************************
void ToyBoxdEdx::Reset(){
//******************************************************************** 

  if (!_ResetCheckDone){
    if( typeid(*this) !=  typeid(ToyBoxdEdx)){
      std::cerr << "ERROR in ToyBoxdEdx::Reset(). Either this mandatory method is not implemented " 
                << "by the derived class '" << typeid(*this).name() << "' "
                << "or ToyBoxdEdx::Reset() is called from the Reset method of the derived class. "
                << "Please correct any of these bugs. " << std::endl;
      
      exit(1);
    }
    _ResetCheckDone=true;
  }  
}

//********************************************************************
void ToyBoxdEdx::ResetBase(){
//******************************************************************** 

  ToyBoxPD::ResetBase();
  Tracks.clear();
  SoftReset();  // just reset internal stuff
}
