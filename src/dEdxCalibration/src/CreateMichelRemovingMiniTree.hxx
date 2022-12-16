#ifndef CreateMichelRemovingMiniTree_h
#define CreateMichelRemovingMiniTree_h

#include "CreateMiniTreePD.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class CreateMichelRemovingMiniTree: public CreateMiniTreePD {
 public:

  CreateMichelRemovingMiniTree(int argc, char *argv[]);
  virtual ~CreateMichelRemovingMiniTree(){}

protected:

  bool Process();

  bool SaveMiniTree();

  void ClearUninterestingHits(AnaParticlePD* part);
  void ClearUninterestingTracks();
  bool IsInterestingHit(AnaHitPD& hit);

public:

};

#endif
