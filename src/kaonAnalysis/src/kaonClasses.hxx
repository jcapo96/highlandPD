#ifndef kaonClasses_hxx
#define kaonClasses_hxx

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


class kaonCounters: public PDCounters{
  
public:
  
  kaonCounters(){
    ntrue_offspring_kaons = 0;
    ntrue_secondary_kaons = 0;
    ntrue_tertiary_kaons = 0;
    ntrue_overthree_kaons = 0;
  }
  virtual ~kaonCounters(){}
  
  Int_t ntrue_offspring_kaons;
  Int_t ntrue_secondary_kaons;
  Int_t ntrue_tertiary_kaons;
  Int_t ntrue_overthree_kaons;
};


#endif
