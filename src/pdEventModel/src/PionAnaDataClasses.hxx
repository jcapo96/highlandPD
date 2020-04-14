#ifndef PionAnaDataClasses_hxx
#define PionAnaDataClasses_hxx

#include "DataClasses.hxx"
#include "ParticleId.hxx"


/// AnaParticle
class AnaParticlePionAna: public AnaParticle{
public :

  AnaParticlePionAna();
  virtual ~AnaParticlePionAna();

  /// Clone this object.
  virtual AnaParticlePionAna* Clone() {
    return new AnaParticlePionAna(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaParticlePionAna(const AnaParticlePionAna& part);

public:

    enum PartTypeEnum {
      kUnknown=0, 
      kShower,
      kTrack
  };
    
  PartTypeEnum Type;

  bool isBeamPart;
};

/// AnaTrueParticle
class AnaTrueParticlePionAna: public AnaTrueParticle{
public :

  AnaTrueParticlePionAna();
  virtual ~AnaTrueParticlePionAna();

  /// Clone this object.
  virtual AnaTrueParticlePionAna* Clone() {
    return new AnaTrueParticlePionAna(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaTrueParticlePionAna(const AnaTrueParticlePionAna& truePart);

public:

  /// Vector Pi0 decays IDs
  std::vector<Int_t> Pi0_decay_ID; 

  Bool_t Matched;
  Int_t  Origin;
  
};


/// Representation of the beam information, including POT and quality.
class AnaBeamPionAna: public AnaBeam{
public :

  AnaBeamPionAna();
  virtual ~AnaBeamPionAna();

  /// Clone this object.
  virtual AnaBeamPionAna* Clone() {
    return new AnaBeamPionAna(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;
  
protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaBeamPionAna(const AnaBeamPionAna& beam);

public:

  /// Other relevant beam info
  int BeamTrigger;
  double TOF;
  int    CerenkovStatus[2];
  double CerenkovTime[2];
  double CerenkovPressure[2];
  double BeamTrackTime;
  double BeamMomentum;
  double BeamMomentumInTPC;

  int nFibers[3];
  size_t nMomenta;  
  size_t nTracks;  
  
  std::vector< int > PDGs;

};



class PionAnaCounters{

public:
  
  PionAnaCounters(){
    ntrue_beamdaughter_piplus=0;
    ntrue_beamdaughter_piminus=0;
    ntrue_beamdaughter_pi0=0;
    ntrue_beamdaughter_proton=0;
    ntrue_beamdaughter_neutron=0;
    ntrue_beamdaughter_nucleus=0;    
  }
  virtual ~PionAnaCounters(){}
  
  Int_t ntrue_beamdaughter_pi0;
  Int_t ntrue_beamdaughter_piplus;
  Int_t ntrue_beamdaughter_piminus;
  Int_t ntrue_beamdaughter_proton;
  Int_t ntrue_beamdaughter_neutron;
  Int_t ntrue_beamdaughter_nucleus;
};


#endif
