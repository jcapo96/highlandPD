#ifndef pdNeutralDataClasses_hxx
#define pdNeutralDataClasses_hxx

#include "DataClasses.hxx"
#include "pdDataClasses.hxx"

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

  };

//** ------------------------------------------------------------ */
// Forward declarations
class AnaTrueEquivalentNeutralParticlePD;

//** ------------------------------------------------------------ */
// Creation vertex class for neutral particle reconstruction
class AnaCreationVertexPD : public AnaVertexPD {
public:
  AnaCreationVertexPD();
  ~AnaCreationVertexPD();

  /// The beam particle
  AnaParticlePD* BeamParticle;

  /// The secondary particle paired with beam
  AnaParticlePD* SecondParticle;

  /// Chi2/ndf of secondary particle under proton hypothesis (lower = more proton-like)
  Float_t ProtonScore;

  /// Distance from beam end position to secondary particle start position
  Float_t DistanceScore;

  /// Closest point on beam fitted line (using creation vertex fit)
  Float_t ClosestPointBeam[3];

  /// Closest point on secondary particle fitted line
  Float_t ClosestPointSecond[3];

  /// Minimum distance between beam end line and secondary particle line
  Float_t MinDistanceScore;

  /// 3D position of creation vertex (midpoint of minimum distance between lines)
  Float_t Position[3];
};

//** ------------------------------------------------------------ */
// Annihilation vertex class for neutral particle reconstruction
class AnaAnnihilationVertexPD : public AnaVertexPD {
public:
  AnaAnnihilationVertexPD();
  ~AnaAnnihilationVertexPD();

  // Inherits all members from AnaVertexPD including Degeneracy
};

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

  /// The annihilation vertex associated with this neutral particle (V-topology)
  AnaAnnihilationVertexPD* AnnihilationVertex;

  /// The creation vertex associated with this neutral particle (beam + secondary)
  AnaCreationVertexPD* CreationVertex;

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

