#ifndef CreateMiniTreeKaon_h
#define CreateMiniTreeKaon_h

#include "CreateMiniTreePD.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

//Migue:
//This class generates minitrees with particles compatible 
//with the candidate definition (descendants of the beam 
//particle with a single daughter)

class CreateMiniTreeKaon: public CreateMiniTreePD {
 public:

  CreateMiniTreeKaon(int argc, char *argv[]);
  virtual ~CreateMiniTreeKaon(){}

protected:

  bool RecoCandidateExists();
  bool TruthCandidateExists();

  void DeleteUninterestingBunches() {return;} //no bunches to delete in protodune

  void DeleteUninterestingParticles();
  void DeleteUninterestingTrueParticles();

  void FilterParticleInfo(AnaParticleB& part) {return;} //no info to filter right now

};

#endif
