#include "CreateMiniTreePionAna.hxx"
#include "pionTreeConverter.hxx"

//********************************************************************
CreateMiniTreePionAna::CreateMiniTreePionAna(int argc, char *argv[]):CreateMiniTree(argc, argv){
  //********************************************************************

  input().AddConverter("pionana",        new pionTreeConverter());
}
