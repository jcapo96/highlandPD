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

  bool Initialize();
  
  bool Process();

  bool SaveMiniTree();

  bool PartHasGoodAngles(AnaParticlePD* part);
  void ClearUninterestingHits(AnaParticlePD* part);
  void ClearUninterestingTracks();
  bool IsInterestingHit(AnaHitPD& hit);

  double _ThetaXZ_min;
  double _ThetaXZ_max;
  double _ThetaYZ_min;
  double _ThetaYZ_max;

public:

};

#endif
