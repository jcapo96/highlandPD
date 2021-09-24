#ifndef kaonCreateMiniTree_h
#define kaonCreateMiniTree_h

#include "CreateMiniTree.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class kaonCreateMiniTree: public CreateMiniTree {
 public:

  kaonCreateMiniTree(int argc, char *argv[]);
  virtual ~kaonCreateMiniTree(){}

protected:

  bool RecoCandidateExists();
  bool TruthCandidateExists();

  void DeleteUninterestingParticles();
  void DeleteUninterestingTrueParticles();

};

#endif
