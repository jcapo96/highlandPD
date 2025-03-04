#ifndef CreateMiniTreeCandidates_h
#define CreateMiniTreeCandidates_h

#include "CreateMiniTreePD.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

//Migue:
//This class generates minitrees only with 
//particles passing the secondarykaonselection

class CreateMiniTreeCandidates: public CreateMiniTreePD {
 public:

  CreateMiniTreeCandidates(int argc, char *argv[]);
  virtual ~CreateMiniTreeCandidates(){}

protected:

  bool Process();

  bool IsCandidate(AnaParticlePD* part);
  bool SaveMiniTree();

  void DeleteAllButCandidates();
};

#endif
