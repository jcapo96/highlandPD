#include "CreateMiniTreePD.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"

//********************************************************************
CreateMiniTreePD::CreateMiniTreePD(int argc, char *argv[]):CreateMiniTree(argc, argv){
//********************************************************************

  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree",true));
  input().AddConverter("pduneana",         new PDSPAnalyzerTreeConverter());
  input().AddConverter("MichelRemovingTree", new MichelRemovingTreeConverter("michelremoving/Event"));
}

//********************************************************************
bool CreateMiniTreePD::Initialize(){
//********************************************************************

  // call the base class
  CreateMiniTree::Initialize();
  
  _saveHits  = (bool)ND::params().GetParameterI("pdIO.MiniTree.SaveHits");
  
  return true;
}

//********************************************************************
void CreateMiniTreePD::FilterParticleInfo(AnaParticleB& partB){
//********************************************************************

  AnaParticlePD& part = *static_cast<AnaParticlePD*>(&partB);

  if (!_saveHits){
    for (size_t i=0;i<3;i++)
      part.Hits[i].clear();
  }
}
