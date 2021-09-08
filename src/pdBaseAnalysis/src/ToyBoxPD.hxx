#ifndef ToyBoxPD_h
#define ToyBoxPD_h

#include "ToyBoxB.hxx"
#include "pdDataClasses.hxx"


/* This is the container used to pass information from one step to another in the selection. It inherits from the base class in highland/src/psyche/psycheCore. 
*/


class ToyBoxPD:public ToyBoxB{
public :
  
  ToyBoxPD();
  virtual ~ToyBoxPD(){}

  /// This method should be implemented by the derived class. If so it does nothing here 
  virtual void Reset();

  /// Reset this base class
  virtual void ResetBase();
  
public:

  /// For storing the true vertex, for analyses with no reconstructed primary vertex
  AnaTrueVertexB* TrueVertex;

  /// The reconstructed Event Vertex
  AnaVertexB* Vertex;

  /// The MainTrack (the beam particle in general)
  AnaParticlePD* MainTrack;

  /// Vector of candidates
  //std::vector<AnaParticlePD*> Candidates;
};

#endif
