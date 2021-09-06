#include "CreateMiniTreePD.hxx"
#include "pdMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
CreateMiniTreePD::CreateMiniTreePD(int argc, char *argv[]):CreateMiniTree(argc, argv){
  //********************************************************************

  input().AddConverter("minitree",         new pdMiniTreeConverter());
  input().AddConverter("pduneana",         new PDSPAnalyzerTreeConverter());
}
