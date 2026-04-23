#ifndef pdNeutralDataClasses_hxx
#define pdNeutralDataClasses_hxx

#include "pdDataClasses.hxx"

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

	/// Minimum distance between beam end line and secondary particle line
	Float_t MinDistanceScore;

	/// 3D position of creation vertex (midpoint of minimum distance between lines)
	Float_t Position[3];

	/// Closest points between beam and second-particle lines
	Float_t ClosestPointBeam[3];
	Float_t ClosestPointSecond[3];
};

//** ------------------------------------------------------------ */
// Annihilation vertex class for neutral particle reconstruction
class AnaAnnihilationVertexPD {
public:
	AnaAnnihilationVertexPD();
	~AnaAnnihilationVertexPD();

	Int_t UniqueID;
	Int_t NParticles;
	std::vector<AnaParticlePD*> Particles;
	Float_t PositionPandora[3];
	Float_t PositionFit[3];
	Float_t ClosestPointPandora1[3];
	Float_t ClosestPointPandora2[3];
	Float_t ClosestPointFit1[3];
	Float_t ClosestPointFit2[3];
	Float_t MinimumDistancePandora;
	Float_t MinimumDistanceFit;
	Float_t OriginalDistance;
	/// Which endpoint pair minimized distance: 0=SS, 1=SE, 2=ES, 3=EE (-1 unset). Used for daughter reversal policy.
	Int_t PairingEndpointCombo;
	Int_t Degeneracy;
	Float_t Momentum;
	Float_t InvariantMass;
	Float_t MomentumPandora;
	Float_t InvariantMassPandora;
	Float_t MomentumFit;
	Float_t InvariantMassFit;
	// Daughter momentum method flags (see pdAnnihilationUtils DaughterMomentumMethod):
	// -1 unassigned, 0 pion-range stopping (legacy), 1 calorimetric (legacy), 2 failed, 3 free-range ML, 4 joint K0s grid
	Int_t Daughter1MomentumMethod;
	Int_t Daughter2MomentumMethod;
	Int_t Daughter1HasPreexistingMomentum;
	Int_t Daughter2HasPreexistingMomentum;
	Int_t Daughter1ExtensionAttempted;
	Int_t Daughter2ExtensionAttempted;
	Int_t Daughter1ExtensionValid;
	Int_t Daughter2ExtensionValid;
	Float_t Daughter1ExtensionChi2Ndf;
	Float_t Daughter2ExtensionChi2Ndf;
	Int_t Daughter1ExtensionNValidHits;
	Int_t Daughter2ExtensionNValidHits;
	Float_t Daughter1ExtensionDedxBias;  // Gaussian mean of (measured - Bethe-Bloch) dEdx [MeV/cm]
	Float_t Daughter2ExtensionDedxBias;  // Gaussian mean of (measured - Bethe-Bloch) dEdx [MeV/cm]
	Float_t Daughter1ExtensionDedxSigma; // Gaussian sigma of (measured - Bethe-Bloch) dEdx [MeV/cm]
	Float_t Daughter2ExtensionDedxSigma; // Gaussian sigma of (measured - Bethe-Bloch) dEdx [MeV/cm]
	Int_t Daughter1ExtensionDedxFitOk;
	Int_t Daughter2ExtensionDedxFitOk;
	/// 0: independent free-range (or joint disabled / joint failed); 1: joint K0s grid+momentum applied.
	Int_t JointK0sMomentumUsed;
	/// Per-daughter momentum summaries [GeV/c] for diagnostics/tree export.
	Float_t Daughter1MomentumTLE;
	Float_t Daughter2MomentumTLE;
	Float_t Daughter1MomentumMCS;
	Float_t Daughter2MomentumMCS;
	Float_t Daughter1MomentumJoint;
	Float_t Daughter2MomentumJoint;
	Float_t JointK0sBestScore;
	Float_t JointK0sInvMassAtBest;
	/// Joint mass-penalty diagnostics (GeV for sigmas; -999 unset). Filled when joint fit is attempted with event σ_m.
	Float_t JointK0sSigmaP1GeV;
	Float_t JointK0sSigmaP2GeV;
	Float_t JointK0sSigmaMEventGeV;
	Float_t JointK0sDmDp1;
	Float_t JointK0sDmDp2;
	/// Post-fit: |C_mass|/max(|logL1+logL2|,eps) at joint optimum (−999 unset).
	Float_t JointK0sMomentumConstraintRatioR;
	/// Post-fit: χ²_TLE(after)−χ²_TLE(before); extension χ² from pdMomReconstruction template (−999 unset).
	Float_t JointK0sMomentumDedxChi2Degradation;
	/// Debug only: 0 unset; 1 data-driven (R<0.3 and Δχ²≤0); 2 constraint-dominated (R≥1); 3 calorimetric tension (Δχ²>5).
	Int_t JointK0sDebugClass;
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

	/// Reconstructed neutral length using annihilation Pandora position.
	Float_t LengthPandora;

	/// Reconstructed neutral length using annihilation Fit position.
	Float_t LengthFit;

	/// Angle (rad) between neutral axis (creation→annihilation Pandora) and Σp̂ vertex momentum (Pandora dirs).
	Float_t AlignmentPandora;

	/// Angle (rad) between neutral axis (creation→annihilation Fit) and Σp̂ vertex momentum (fit dirs).
	Float_t AlignmentFit;

};

#endif

