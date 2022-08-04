#ifndef CreateMiniTreeCosmics_h
#define CreateMiniTreeCosmics_h

#include "CreateMiniTreePD.hxx"

class CreateMiniTreeCosmics: public CreateMiniTreePD {
 public:

  CreateMiniTreeCosmics(int argc, char *argv[]);
  virtual ~CreateMiniTreeCosmics(){}

protected:
  bool Initialize();
  bool CheckSaveParticle(const AnaParticleB& part);
  
protected:

  Float_t _cutLength;
  Float_t _cutStartZmin;
  Float_t _cutStartZmax;
};

#endif
