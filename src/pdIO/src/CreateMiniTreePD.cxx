#include "CreateMiniTreePD.hxx"
#include "pdMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
CreateMiniTreePD::CreateMiniTreePD(int argc, char *argv[]):CreateMiniTree(argc, argv){
  //********************************************************************

  input().AddConverter("pionana",          new hitPionTreeConverter());
  input().AddConverter("minitree",         new pdMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("pduneana",         new PDSPAnalyzerTreeConverter());
}
