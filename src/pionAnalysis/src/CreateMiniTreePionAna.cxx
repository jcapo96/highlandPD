#include "CreateMiniTreePionAna.hxx"
#include "pionTreeConverter.hxx"
//#include "LArSoftTreeConverter.hxx"

//********************************************************************
CreateMiniTreePionAna::CreateMiniTreePionAna(int argc, char *argv[]):CreateMiniTree(argc, argv){
  //********************************************************************

  input().AddConverter("pionana",        new pionTreeConverter());
  //  input().AddConverter("LArSoftTree",    new LArSoftTreeConverter());
}
