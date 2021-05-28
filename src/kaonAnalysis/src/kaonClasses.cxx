#define kaonClasses_C

#include "kaonClasses.hxx"

const Double_t kDoubleUnassigned = -999.;
const Int_t    kIntUnassigned = -999;
const Float_t  kFloatUnassigned = -999.;

//********************************************************************
kaonAnaTrueVertex::kaonAnaTrueVertex():AnaTrueVertex(){
//********************************************************************

  Branch    = kIntUnassigned;
  DecayMode = kIntUnassigned;
  ChainMuon = kIntUnassigned;
}

//********************************************************************
kaonAnaTrueVertex::kaonAnaTrueVertex(const kaonAnaTrueVertex& vertex):AnaTrueVertex(vertex){
//********************************************************************

  Branch    = vertex.Branch;
  DecayMode = vertex.DecayMode;
  ChainMuon = vertex.ChainMuon;
}

//********************************************************************
void kaonAnaTrueVertex::Print() const{
//********************************************************************

  std::cout << "-------- kaonAnaTrueVertex --------- " << std::endl;

  AnaTrueVertex::Print();

  std::cout << "Branch:           " << Branch << std::endl;
  std::cout << "Decay Mode:       " << DecayMode << std::endl;
  std::cout << "ChainMuon:        " << ChainMuon << std::endl;
}

