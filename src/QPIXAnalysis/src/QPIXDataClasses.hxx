#ifndef QPIXDataClasses_hxx
#define QPIXDataClasses_hxx

#include "pdDataClasses.hxx"
#include "ParticleId.hxx"

// Extension of AnaEventInfo
class AnaEventInfoQPIX: public AnaEventInfoPD{
public :

  AnaEventInfoQPIX();
  virtual ~AnaEventInfoQPIX();

  /// Clone this object.
  virtual AnaEventInfoQPIX* Clone() {
    AnaEventInfoQPIX* eventinfo = new AnaEventInfoQPIX(*this);
    return eventinfo;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaEventInfoQPIX(const AnaEventInfoQPIX& eventinfo);

public:

  enum SolarChainEnum {
    unknown_chain = 0,
    pp,
    Be7,
    pep,
    N13,
    O15,
    F17,
    B8,
    hep,
    dummy_chain
  };

  enum SolarNuReactionEnum {
    unknown_reaction = 0,
    ES,
    CC,
    dummy_reaction
  };

  SolarChainEnum SolarChain;
  SolarNuReactionEnum NuReaction;

  bool IsBackground;
  bool HasPotasium;
  bool HasGamma;
  int  DetectedPhotons[7];
  std::vector<int> PDSwf[7];
  
public:
  
  std::string ChainToString(const int index) const;
  std::string NuReactionToString(const int index) const;

};

#endif
