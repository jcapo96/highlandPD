#ifndef ToyBoxPD_h
#define ToyBoxPD_h

#include "ToyBoxB.hxx"
#include "pdDataClasses.hxx"

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

  /// The reconstructed EventVertex
  AnaVertexB* Vertex;

  /// The MainTrack (te beam particle in general)
  AnaParticlePD* MainTrack;

  // Distance from the MainTrack to its daughters
  std::vector<Float_t> DaughterDistanceToVertex;
};

#endif
