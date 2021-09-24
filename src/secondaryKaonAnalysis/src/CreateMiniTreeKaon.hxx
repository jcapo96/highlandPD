#ifndef CreateMiniTreeKaon_h
#define CreateMiniTreeKaon_h

#include "CreateMiniTreePD.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class CreateMiniTreeKaon: public CreateMiniTreePD {
 public:

  CreateMiniTreeKaon(int argc, char *argv[]);
  virtual ~CreateMiniTreeKaon(){}

protected:

  bool RecoCandidateExists();
  bool TruthCandidateExists();

  void DeleteUninterestingParticles();
  void DeleteUninterestingTrueParticles();

};

#endif
