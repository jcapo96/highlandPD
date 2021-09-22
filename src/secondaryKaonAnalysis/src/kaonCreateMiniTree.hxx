#ifndef kaonCreateMiniTree_h
#define kaonCreateMiniTree_h

#include "SimpleLoopBase.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class kaonCreateMiniTree: public SimpleLoopBase {
 public:

  kaonCreateMiniTree(int argc, char *argv[]);
  virtual ~kaonCreateMiniTree(){}

protected:

  //---- These are mandatory functions
  virtual bool Initialize     ();
  virtual bool InitializeSpill(){return true;}

  virtual void DefineOutputTree();

  virtual void FinalizeSpill(){}
  virtual void Finalize     ();

  virtual bool Process();
  //--------------------

  virtual void FillMiniTree();

  bool RecoCandidateExists();
  bool TruthCandidateExists();

  void DeleteUninterestingParticles();
  void DeleteUninterestingTrueParticles();

  enum miniTreeIndex{
    minitree = OutputManager::enumSpecialTreesLast+1
  };


protected:

  Int_t _totalRecoObjects;
  Int_t _savedRecoObjects;
  Int_t _totalTrueObjects;
  Int_t _savedTrueObjects;

  bool _saveCandidates;

  AnaSpill* _spill;
};

#endif
