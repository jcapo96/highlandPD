#ifndef CreateMiniTreeEventCandidates_h
#define CreateMiniTreeEventCandidates_h

#include "CreateMiniTreePD.hxx"
#include "pdDataClasses.hxx"
#include "InputManager.hxx"
#include <set>

class CreateMiniTreeEventCandidates: public CreateMiniTreePD {
 public:

  CreateMiniTreeEventCandidates(int argc, char *argv[]);
  virtual ~CreateMiniTreeEventCandidates(){}

protected:

  bool Process();

  bool IsCandidate(AnaParticlePD* part);
  bool SaveMiniTree();

  void DeleteAllButCandidates();

public:

};

#endif
