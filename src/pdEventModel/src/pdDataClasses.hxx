#ifndef pdDataClasses_hxx
#define pdDataClasses_hxx

#include "DataClasses.hxx"
#include "ParticleId.hxx"

const UInt_t NMAXHITSPERPLANE        = 300;
const UInt_t NMAXHITSPERPLANE_SELTRK = 2500;

class AnaWireID{
public:

  AnaWireID(){Wire=Plane=TPC=Cryostat=0;}
  AnaWireID(const AnaWireID& wireID);
  virtual ~AnaWireID(){}
  
  Int_t Wire;
  Int_t Plane;
  Int_t TPC;
  Int_t Cryostat;
};

class AnaHitPD{
public:
  AnaHitPD();
  virtual ~AnaHitPD(){}

  //Anselmo stuff for CNN calculation, not really used, to be removed at some point
  AnaHitPD(Int_t wire, Float_t integ, Float_t peakT, Float_t peakAmp, TVector3 pos){
    WireID.Plane =wire;
    Integral = integ;
    PeakTime = peakT;
    PeakAmplitude = peakAmp;
    Position = pos;
  }
    
  /// Dump the object to screen.
  virtual void Print() const;

  ///   Copy constructor
  AnaHitPD(const AnaHitPD& part);
  
public:

  // ---- Most of this information is needed to recompute the CNN. This is work in progress (is ABANDONED work)
  AnaWireID WireID;    

  Float_t Integral;
  Float_t PeakTime;
  Float_t PeakAmplitude;
  int TPCid;
  int PlaneID;

  UInt_t StartTick;
  UInt_t EndTick;
  
  UInt_t Channel;
  Int_t  View;

  //std::vector<Float_t> CNN; //!
  Float_t CNN[3];
  
  /// wave form associated to this hit
  std::vector<Float_t> Signal;  

  //------------------------------------------------------------

  /// Calorimetric information
  /// No SCE correction
  Float_t dQdx_NoSCE;
  Float_t dEdx_NoSCE;
  Float_t ResidualRange_NoSCE;
  Float_t Pitch_NoSCE;
  TVector3 Position_NoSCE;
  TVector3 Direction_NoSCE;


  /// SCE correction
  Float_t dQdx_SCE;
  Float_t dEdx_SCE;
  Float_t ResidualRange_SCE;

  /// SCE correction + e⁻ lifetime correction
  Float_t dQdx_elife;
  Float_t dEdx_elife;
  Float_t ResidualRange_elife;

  /// SCE correction + e⁻ lifetime correction + XYZT correction + dEdx calibration
  Float_t dQdx;
  Float_t dEdx;
  Float_t ResidualRange;
  Float_t Pitch;
  TVector3 Position;
  TVector3 Direction;

  /// deprecated but keep it for the moment to avoid errors
  Float_t dEdx_calib;

};

//trajectory point representation. Has position and direction information.
//this information is somehow redundant since each hit has a tjp associated 
//from which position and direction are obtained. This should be solved at
//some point. However, all tjp are needed to compute track lenght and sce
//systematics
class AnaTrajectoryPointPD{
public:
  AnaTrajectoryPointPD();
  virtual ~AnaTrajectoryPointPD(){}

  /// Dump the object to screen.
  virtual void Print() const;

  /// Copy constructor
  AnaTrajectoryPointPD(const AnaTrajectoryPointPD& tjp);
  
public:

  TVector3 Position;
  TVector3 Direction;

  TVector3 Position_NoSCE;
  TVector3 Direction_NoSCE;

  bool IsValid() {return (Position_NoSCE.X() != -999 && Position_NoSCE.Y() != -999 && Position_NoSCE.Z() != -999);}
};

/// AnaParticle
class AnaParticlePD: public AnaParticle{
public :

  AnaParticlePD();
  virtual ~AnaParticlePD();

  /// Clone this object.
  virtual AnaParticlePD* Clone() {
    return new AnaParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaParticlePD(const AnaParticlePD& part);

public:

    enum PartTypeEnum {
      kUnknown=0, 
      kShower,
      kTrack
  };

  /// MotherID
  Int_t ParentID;

  /// for kaon analysis
  bool IsCandidate;

  /// Track or Shower
  PartTypeEnum Type;

  /// Forced Track ID
  Int_t TrackID;

  /// Forced shower ID
  Int_t ShowerID;

  /// Is this the beam particle according to Pandora ?
  bool isPandora;

  /// Is this the beam particle or a descendant according to Pandora ?
  bool BeamOrigin;

  ///  Pandora beam particle that passes geometric cuts
  bool isBeamPart;

  /// Initial position and direction SCE corrected
  double PositionStartSCE[4];
  double PositionEndSCE[4];
  double DirectionStartSCE[3];
  double DirectionEndSCE[3];

  double ThetaXZ;
  double ThetaYZ;

  /// Vector of hits for each plane
  std::vector<AnaHitPD> Hits[3];

  /// Total Number of hits in each wire plane (The vector of Hits above might be a subsample)
  Int_t NHitsPerPlane[3];

  // Vector of trajectory points
  std::vector<AnaTrajectoryPointPD> TrjPoints;

  /// Libo truncated mean
  Float_t truncLibo_dEdx;
  
  // PID variables
  Float_t Chi2Proton;
  Float_t Chi2Muon;
  Float_t Chi2ndf;

  Float_t CNNscore[3];
  Float_t vtx_CNN_michelscore;
  Int_t   vtx_CNN_NHits;
  
  /// CALO variables
  Float_t CALO[3][10];
  
  /// Momentum by range for muon and proton hypotheses
  Float_t RangeMomentum[2];
  Float_t RangeMomentum_alt[2];

  /// Alternate length
  Float_t Length_alt;

  /// Generation of the particle (primary particle, daughter, granddaughter, etc.)
  Int_t Generation;

  /// test
  Double_t Distance_to_closest_particle;

  // ---- OBSOLETE PID VARIABLES ----------
  
  /// Particle ID hypothesis used in the fit (if any)
  Int_t FitPDG;
  
  /// PDG of the most probable particle hypothesis used at reconstruction level
  Int_t ReconPDG[3]; 

  /// PID variables
  Float_t PID[3][10];

  /// PIDA
  Float_t PIDA[3]; 

  /// Migue test
  bool forced_daughter;
  bool forced_daughter_matched;
};

/// AnaTrueParticle
class AnaTrueParticlePD: public AnaTrueParticle{
public :

  AnaTrueParticlePD();
  virtual ~AnaTrueParticlePD();

  /// Clone this object.
  virtual AnaTrueParticlePD* Clone() {
    return new AnaTrueParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaTrueParticlePD(const AnaTrueParticlePD& truePart);

public:

  /// Vector Pi0 decays IDs
  std::vector<Int_t> Pi0_decay_ID; 

  /// True-reco matching flag, need for 
  Bool_t Matched;

  /// Origin
  Int_t  Origin;

  /// Generation
  Int_t Generation;

  /// The particle length inside the TPC
  Float_t LengthInTPC;
  
  /// The true momentum at the TPC entrance
  Float_t MomentumInTPC;  
};

/// Extension of AnaEvent to include specific information of the ProtoDUNE beam line instrumentation
class AnaBeamPD: public AnaBeam{
public :

  AnaBeamPD();
  virtual ~AnaBeamPD();

  /// Clone this object.
  virtual AnaBeamPD* Clone() {
    return new AnaBeamPD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;
  
protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaBeamPD(const AnaBeamPD& beam);

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


/// ProtoDUNE counters

class PDCounters{

public:
  
  PDCounters(){
    ntrue_beamdaughter_piplus=0;
    ntrue_beamdaughter_piminus=0;
    ntrue_beamdaughter_pi0=0;
    ntrue_beamdaughter_proton=0;
    ntrue_beamdaughter_neutron=0;
    ntrue_beamdaughter_nucleus=0;    
  }
  virtual ~PDCounters(){}
  
  Int_t ntrue_beamdaughter_pi0;
  Int_t ntrue_beamdaughter_piplus;
  Int_t ntrue_beamdaughter_piminus;
  Int_t ntrue_beamdaughter_proton;
  Int_t ntrue_beamdaughter_neutron;
  Int_t ntrue_beamdaughter_nucleus;
};

class AnaWireCNN{
public:
  AnaWireCNN();
  virtual ~AnaWireCNN(){}
  AnaWireCNN(const AnaWireCNN& wire);

public:

  std::vector<float> adcs;
  Int_t wire;
  Int_t time;
};


// Extension of AnaBunch to include the APA wire wafeforms, needed to recompute the CNN 
class AnaBunchPD: public AnaBunch{
public :
  
  AnaBunchPD();
  virtual ~AnaBunchPD();
  
  /// Clone this object.
  virtual AnaBunchPD* Clone() {
    return new AnaBunchPD(*this);
  }
  
  /// Dump the object to screen.
  virtual void Print() const;
  
protected:
  
  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaBunchPD(const AnaBunchPD& bunch);
  
public:

  std::vector<AnaWireCNN> CNNwires;
};

class AnaSpillPD: public AnaSpill{
public :

  AnaSpillPD(){}
  virtual ~AnaSpillPD(){}

  /// Clone this object.
  virtual AnaSpillPD* Clone() {
    AnaSpillPD* spill = new AnaSpillPD(*this);
    spill->RedoLinks();
    spill->isClone=true;
    return spill;
  }

protected:

public:

};

// Extension of AnaEvent to include the APA wire wafeforms, needed to recompute the CNN
class AnaEventPD: public AnaEvent{
public :

  AnaEventPD();
  AnaEventPD(const AnaSpillPD& spill, const AnaBunchPD& bunch);
  virtual ~AnaEventPD();

  /// Clone this object.
  virtual AnaEventPD* Clone() {
    AnaEventPD* spill = new AnaEventPD(*this);
    spill->isClone=true;
    return spill;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaEventPD(const AnaEventPD& event);

public:

  std::vector<AnaWireCNN> CNNwires;
};

// Extension of AnaEventInfo
class AnaEventInfoPD: public AnaEventInfo{
public :

  AnaEventInfoPD();
  virtual ~AnaEventInfoPD();

  /// Clone this object.
  virtual AnaEventInfoPD* Clone() {
    AnaEventInfoPD* eventinfo = new AnaEventInfoPD(*this);
    return eventinfo;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaEventInfoPD(const AnaEventInfoPD& eventinfo);

public:

  Float_t NominalBeamMom;

  Bool_t EmptyEvent; //this should probably go in DataQuality
  Bool_t HasPandora; //this should probably go in DataQuality
};



#endif
