#ifndef CreateMiniTreePD_h
#define CreateMiniTreePD_h

#include "CreateMiniTree.hxx"
#include "pdDataClasses.hxx"

class CreateMiniTreePD: public CreateMiniTree {
 public:

  CreateMiniTreePD(int argc, char *argv[]);
  virtual ~CreateMiniTreePD(){}

protected:

  virtual bool Initialize();
  void FilterParticleInfo(AnaParticleB& part);
  
protected:
  Bool_t _saveHits;  
};

#endif
