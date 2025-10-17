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

//** ------------------------------------------------------------ */
// True equivalent vertex class for ProtoDUNE
class AnaTrueEquivalentVertexPD{
  public :

    AnaTrueEquivalentVertexPD();
    virtual ~AnaTrueEquivalentVertexPD();

    /// Clone this object.
    virtual AnaTrueEquivalentVertexPD* Clone() {
      return new AnaTrueEquivalentVertexPD(*this);
    }

    /// Dump the object to screen.
    virtual void Print() const;

  protected:

    /// Copy constructor is protected, as Clone() should be used to copy this object.
    AnaTrueEquivalentVertexPD(const AnaTrueEquivalentVertexPD& vertex);

  public:

    /// Vector containing the true particles associated with this vertex
    std::vector<AnaTrueParticlePD*> TrueParticles;

    /// True Original distance between the two true particles
    Float_t OriginalDistance;

    /// Minimum distance between the two true particles
    Float_t MinimumDistance;

    /// Opening angle between the two true particles
    Float_t OpeningAngle;

    /// 3D coordinates of the vertex
    Float_t Position[3];

    /// 3D direction of the vertex
    Float_t Direction[3];

  /// Pandora-based vertex position (from DirectionStart/PositionStart)
  Float_t PositionPandora[3];

  /// Fitted vertex position (from geometric/Kalman fit)
  Float_t PositionFit[3];

  /// Fitted vertex direction (from geometric/TMinuit/Kalman fit)
  Float_t DirectionFit[3];

  /// Flag: 1 if Pandora calculation used simple average, 0 if used line intersection
  Int_t IsJustAverage;

  /// Degeneracy count before scoring (vertices within threshold)
  Int_t DegeneracyBeforeScoring;

  /// Degeneracy count after scoring (vertices within threshold)
  Int_t DegeneracyAfterScoring;

  /// Total unique particles in degenerate vertices
  Int_t NRecoParticles;

  /// Vector storing the 5 minimum distances to vertices within DegeneracyDistance threshold
  std::vector<Float_t> DegeneracyDistances;

  /// Vector storing the 5 minimum line-to-point distances for isolation calculation (Pandora-based)
  std::vector<Float_t> IsolationDistances;

  /// Vector storing the 5 minimum line-to-point distances using fitted tracks
  std::vector<Float_t> IsolationDistancesFit;

  /// Vector storing the 5 minimum point-to-point distances from vertex to particle PositionStart
  std::vector<Float_t> IsolationStartDistances;

  /// Total number of isolation particles that are true protons
  Int_t IsolationNProton;

  /// Total number of isolation particles that are true pions (both charges)
  Int_t IsolationNPion;

  /// Vector of flags (1=true proton, 0=not) for the 5 closest particles
  std::vector<Int_t> IsolationIsProton;

  /// Vector of chi2/ndf under proton hypothesis for the 5 closest particles
  std::vector<Float_t> IsolationChi2Proton;

  /// Vector of lengths for the 5 closest particles
  std::vector<Float_t> IsolationLength;

  /// Flag: 1 if any isolation particle is longer than both vertex particles, 0 otherwise
  Int_t IsolationIsLongest;
  };

//** ------------------------------------------------------------ */
// Extension of AnaVertexB for ProtoDUNE analysis
class AnaVertexPD: public AnaVertexB{
public :

  AnaVertexPD();
  virtual ~AnaVertexPD();

  /// Clone this object.
  virtual AnaVertexPD* Clone() {
    return new AnaVertexPD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaVertexPD(const AnaVertexPD& vertex);

public:

  /// Unique ID for this vertex within the event
  Int_t UniqueID;

  /// Number of particles associated with this vertex
  Int_t NParticles;

  /// Vector containing the particles associated with this vertex
  std::vector<AnaParticlePD*> Particles;

  /// 3D coordinates of the vertex
  Float_t Position[3];

  /// 3D momentum of the vertex
  Float_t Momentum[3];

  /// 3D direction of the vertex
  Float_t Direction[3];

  /// Energy of the vertex
  Float_t Energy;

  /// Opening angle of the vertex
  Float_t OpeningAngle;

  /// Angle with respect to the beam direction
  Float_t AngleWithBeam;

  /// Number of potential parents
  Int_t NPotentialParents;

  /// Generation of the vertex (1=primary, 2=secondary, etc.)
  Int_t Generation;

  /// Reaction/process that originated this vertex (interaction type)
  Int_t Process;

  /// Original distance between the particles in the vertex
  Float_t OriginalDistance;

  /// Fitted line parameters for daughter particles used in vertex reconstruction
  /// Each line is represented by 6 parameters: [x0, y0, z0, dx, dy, dz]
  /// Line equation: P(t) = (x0, y0, z0) + t * (dx, dy, dz)
  std::vector<std::vector<double>> FittedLineParams;

  /// Minimum distance between the two fitted lines (distance between closest points)
  Float_t MinimumDistance;

  /// Quality score from vertex fit (lower is better, from Chi2/minimization)
  Float_t Score;

  /// Parent ID of the vertex (both particles must have same parent)
  Int_t ParentID;

  /// Pandora-based vertex position (from DirectionStart/PositionStart)
  Float_t PositionPandora[3];

  /// Fitted vertex position (from geometric/Kalman fit)
  Float_t PositionFit[3];

  /// Fitted vertex direction (from geometric/TMinuit/Kalman fit)
  Float_t DirectionFit[3];

  /// Flag: 1 if Pandora calculation used simple average, 0 if used line intersection
  Int_t IsJustAverage;

  /// Degeneracy count before scoring (vertices within threshold)
  Int_t DegeneracyBeforeScoring;

  /// Degeneracy count after scoring (vertices within threshold)
  Int_t DegeneracyAfterScoring;

  /// Total unique particles in degenerate vertices
  Int_t NRecoParticles;

  /// Vector storing the 5 minimum distances to vertices within DegeneracyDistance threshold
  std::vector<Float_t> DegeneracyDistances;

  /// Vector storing the 5 minimum line-to-point distances for isolation calculation (Pandora-based)
  std::vector<Float_t> IsolationDistances;

  /// Vector storing the 5 minimum line-to-point distances using fitted tracks
  std::vector<Float_t> IsolationDistancesFit;

  /// Vector storing the 5 minimum point-to-point distances from vertex to particle PositionStart
  std::vector<Float_t> IsolationStartDistances;

  /// Total number of isolation particles that are true protons
  Int_t IsolationNProton;

  /// Total number of isolation particles that are true pions (both charges)
  Int_t IsolationNPion;

  /// Vector of flags (1=true proton, 0=not) for the 5 closest particles
  std::vector<Int_t> IsolationIsProton;

  /// Vector of chi2/ndf under proton hypothesis for the 5 closest particles
  std::vector<Float_t> IsolationChi2Proton;

  /// Vector of lengths for the 5 closest particles
  std::vector<Float_t> IsolationLength;

  /// Flag: 1 if any isolation particle is longer than both vertex particles, 0 otherwise
  Int_t IsolationIsLongest;

  AnaTrueEquivalentVertexPD* TrueEquivalentVertex;

  /// Ensure particles in this vertex have proper momentum assigned
  /// Uses the most robust momentum calculation method available
  void EnsureParticleMomentum();
};


//** ------------------------------------------------------------ */
// Forward declarations
class AnaTrueEquivalentNeutralParticlePD;

//** ------------------------------------------------------------ */
// Extension for neutral particle analysis in ProtoDUNE
class AnaNeutralParticlePD: public AnaParticleB{
public :

  AnaNeutralParticlePD();
  virtual ~AnaNeutralParticlePD();

  /// Clone this object.
  virtual AnaNeutralParticlePD* Clone() {
    return new AnaNeutralParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaNeutralParticlePD(const AnaNeutralParticlePD& neutralParticle);

public:

  /// Unique ID for this neutral particle within the event
  Int_t UniqueID;

  /// The vertex associated with this neutral particle
  AnaVertexPD* Vertex;

  /// The parent particle that decayed into this neutral particle
  AnaParticlePD* Parent;

  /// The impact parameter of the neutral particle
  Float_t ImpactParameter;

  /// Mass of the neutral particle (in GeV/c²)
  Float_t Mass;

  /// Momentum of the neutral particle (in GeV/c)
  Float_t Momentum;

  /// PDG code of the neutral particle
  Int_t PDG;

  /// Lifetime of the neutral particle (in ns)
  Float_t Lifetime;

  /// Decay length of the neutral particle (in cm)
  Float_t DecayLength;

  /// Number of reconstructed hits in the vertex
  Int_t NRecoHitsInVertex;

  /// Score for neutral particle compatibility (lower = more neutral-like)
  Double_t NeutralScore;

  /// Alignment between hits in cylinder and neutral particle direction (dot product)
  Double_t HitsAlignment;

  /// Number of hits in cylinder around neutral particle path
  Int_t NHitsInCylinder;

  /// Average perpendicular distance of hits to neutral particle path
  Double_t HitsAvgDistance;

  /// RMS of perpendicular distances of hits
  Double_t HitsRMSDistance;

  /// Longitudinal span fraction (span along path / total length)
  Double_t HitsLongitudinalSpan;

  /// Number of parent's daughters that are true protons with trajectory close to neutral start position
  Int_t NProtonInCreationVtx;

  /// Total number of parent's daughters near neutral creation vertex (all particle types)
  Int_t NParticlesInCreationVtx;

  /// Chi2/ndf under proton hypothesis for particles near creation vertex (up to 5)
  std::vector<Float_t> CreationVtxChi2Proton;

  /// Minimum distances from Pandora-based lines to neutral start (up to 5)
  std::vector<Float_t> CreationVtxDistances;

  /// True PDG codes of particles near creation vertex (up to 5)
  std::vector<Int_t> CreationVtxTruePDG;

  /// The reconstructed neutral particle associated with this neutral particle
  AnaParticlePD* RecoParticle;

  /// The true neutral particle associated with this reconstructed neutral particle
  AnaTrueEquivalentNeutralParticlePD* TrueEquivalentNeutralParticle;
};

//** ------------------------------------------------------------ */
// True neutral particle class for ProtoDUNE
class AnaTrueEquivalentNeutralParticlePD{
public:

  AnaTrueEquivalentNeutralParticlePD();
  virtual ~AnaTrueEquivalentNeutralParticlePD();

  /// Clone this object.
  virtual AnaTrueEquivalentNeutralParticlePD* Clone() {
    return new AnaTrueEquivalentNeutralParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaTrueEquivalentNeutralParticlePD(const AnaTrueEquivalentNeutralParticlePD& trueEquivalentNeutralPart);

public:

  /// The true vertex associated with this neutral particle
  AnaTrueEquivalentVertexPD* TrueEquivalentVertex;

  /// The true parent particle that decayed into this neutral particle
  AnaTrueParticlePD* TrueParent;

  /// 3D coordinates of the particle
  Float_t Position[3];

  /// 3D direction of the particle
  Float_t Direction[3];

  /// 3D end position of the particle
  Float_t PositionEnd[3];

  /// 3D end direction of the particle
  Float_t DirectionEnd[3];

  /// Length of the particle
  Float_t Length;

  /// Momentum of the particle
  Float_t Momentum;

  /// End momentum of the particle
  Float_t MomentumEnd;

  // Mass of the particle
  Float_t Mass;

  /// PDG code of the particle
  Int_t PDG;

  /// Generation of the particle
  Int_t Generation;

  /// Process of the particle
  Int_t Process;
};



#endif
