#ifndef kaonDataClasses_hxx
#define kaonDataClasses_hxx

#include "pdDataClasses.hxx"
#include "ParticleId.hxx"

class kaonAnaTrueVertex: public AnaTrueVertex{
public :

  kaonAnaTrueVertex();
  virtual ~kaonAnaTrueVertex(){}

  /// Clone this object.
  virtual kaonAnaTrueVertex* Clone() {
    return new kaonAnaTrueVertex(*this);
  }

  /// Dump the object to screen.
  void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  kaonAnaTrueVertex(const kaonAnaTrueVertex& vertex);

public:

  Int_t Branch;
  Int_t DecayMode;
  Int_t ChainMuon;
};


#endif
