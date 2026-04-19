#define pdNeutralDataClasses_C

#include "pdNeutralDataClasses.hxx"

namespace {

const Int_t   kIntUnassigned = -999;
const Float_t kFloatUnassigned = -999.;

}

//********************************************************************
AnaCreationVertexPD::AnaCreationVertexPD() : AnaVertexPD(){
//********************************************************************
	BeamParticle = NULL;
	SecondParticle = NULL;
	ProtonScore = -999.0;
	DistanceScore = -999.0;
	MinDistanceScore = -999.0;
	Position[0] = -999.0;
	Position[1] = -999.0;
	Position[2] = -999.0;
	ClosestPointBeam[0] = kFloatUnassigned;
	ClosestPointBeam[1] = kFloatUnassigned;
	ClosestPointBeam[2] = kFloatUnassigned;
	ClosestPointSecond[0] = kFloatUnassigned;
	ClosestPointSecond[1] = kFloatUnassigned;
	ClosestPointSecond[2] = kFloatUnassigned;
}

//********************************************************************
AnaCreationVertexPD::~AnaCreationVertexPD(){
//********************************************************************
}

//********************************************************************
AnaAnnihilationVertexPD::AnaAnnihilationVertexPD(){
//********************************************************************
	UniqueID = kIntUnassigned;
	NParticles = 0;
	Particles.clear();
	PositionPandora[0] = kFloatUnassigned;
	PositionPandora[1] = kFloatUnassigned;
	PositionPandora[2] = kFloatUnassigned;
	PositionFit[0] = kFloatUnassigned;
	PositionFit[1] = kFloatUnassigned;
	PositionFit[2] = kFloatUnassigned;
	ClosestPointPandora1[0] = kFloatUnassigned;
	ClosestPointPandora1[1] = kFloatUnassigned;
	ClosestPointPandora1[2] = kFloatUnassigned;
	ClosestPointPandora2[0] = kFloatUnassigned;
	ClosestPointPandora2[1] = kFloatUnassigned;
	ClosestPointPandora2[2] = kFloatUnassigned;
	ClosestPointFit1[0] = kFloatUnassigned;
	ClosestPointFit1[1] = kFloatUnassigned;
	ClosestPointFit1[2] = kFloatUnassigned;
	ClosestPointFit2[0] = kFloatUnassigned;
	ClosestPointFit2[1] = kFloatUnassigned;
	ClosestPointFit2[2] = kFloatUnassigned;
	MinimumDistancePandora = kFloatUnassigned;
	MinimumDistanceFit = kFloatUnassigned;
	OriginalDistance = kFloatUnassigned;
	PairingEndpointCombo = -1;
	Degeneracy = kIntUnassigned;
	Momentum = kFloatUnassigned;
	InvariantMass = kFloatUnassigned;
	MomentumPandora = kFloatUnassigned;
	InvariantMassPandora = kFloatUnassigned;
	MomentumFit = kFloatUnassigned;
	InvariantMassFit = kFloatUnassigned;
	Daughter1MomentumMethod = -1;
	Daughter2MomentumMethod = -1;
	Daughter1HasPreexistingMomentum = -1;
	Daughter2HasPreexistingMomentum = -1;
	Daughter1ExtensionAttempted = -1;
	Daughter2ExtensionAttempted = -1;
	Daughter1ExtensionValid = -1;
	Daughter2ExtensionValid = -1;
	Daughter1ExtensionChi2Ndf = kFloatUnassigned;
	Daughter2ExtensionChi2Ndf = kFloatUnassigned;
	Daughter1ExtensionNValidHits = -1;
	Daughter2ExtensionNValidHits = -1;
	Daughter1ExtensionDedxBias = kFloatUnassigned;
	Daughter2ExtensionDedxBias = kFloatUnassigned;
	Daughter1ExtensionDedxSigma = kFloatUnassigned;
	Daughter2ExtensionDedxSigma = kFloatUnassigned;
	Daughter1ExtensionDedxFitOk = -1;
	Daughter2ExtensionDedxFitOk = -1;
	JointK0sMomentumUsed = 0;
	JointK0sBestScore = kFloatUnassigned;
	JointK0sInvMassAtBest = kFloatUnassigned;
	JointK0sSigmaP1GeV = kFloatUnassigned;
	JointK0sSigmaP2GeV = kFloatUnassigned;
	JointK0sSigmaMEventGeV = kFloatUnassigned;
	JointK0sDmDp1 = kFloatUnassigned;
	JointK0sDmDp2 = kFloatUnassigned;
	JointK0sMomentumConstraintRatioR = kFloatUnassigned;
	JointK0sMomentumDedxChi2Degradation = kFloatUnassigned;
	JointK0sDebugClass = 0;
}

//********************************************************************
AnaAnnihilationVertexPD::~AnaAnnihilationVertexPD(){
//********************************************************************
}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(): AnaParticleB(){
//********************************************************************

	UniqueID = kIntUnassigned;
	AnnihilationVertex = NULL;
	CreationVertex = NULL;
	Parent = NULL;
	LengthPandora = kFloatUnassigned;
	LengthFit = kFloatUnassigned;
	AlignmentPandora = kFloatUnassigned;
	AlignmentFit = kFloatUnassigned;
}

//********************************************************************
AnaNeutralParticlePD::~AnaNeutralParticlePD(){
//********************************************************************

}

//********************************************************************
AnaNeutralParticlePD::AnaNeutralParticlePD(const AnaNeutralParticlePD& neutralParticle): AnaParticleB(neutralParticle){
//********************************************************************

	UniqueID = neutralParticle.UniqueID;
	AnnihilationVertex = neutralParticle.AnnihilationVertex;
	CreationVertex = neutralParticle.CreationVertex;
	Parent = neutralParticle.Parent;
	LengthPandora = neutralParticle.LengthPandora;
	LengthFit = neutralParticle.LengthFit;
	AlignmentPandora = neutralParticle.AlignmentPandora;
	AlignmentFit = neutralParticle.AlignmentFit;
}

//********************************************************************
void AnaNeutralParticlePD::Print() const{
//********************************************************************

	std::cout << "-------- AnaNeutralParticlePD --------- " << std::endl;

	AnaParticleB::Print();

	std::cout << "UniqueID:              " << UniqueID << std::endl;
	std::cout << "AnnihilationVertex:    " << (AnnihilationVertex ? "Yes" : "No") << std::endl;
	std::cout << "CreationVertex:        " << (CreationVertex ? "Yes" : "No") << std::endl;
	std::cout << "Parent:                " << (Parent ? "Yes" : "No") << std::endl;
	std::cout << "LengthPandora:         " << LengthPandora << " cm" << std::endl;
	std::cout << "LengthFit:             " << LengthFit << " cm" << std::endl;
	std::cout << "AlignmentPandora:      " << AlignmentPandora << " rad" << std::endl;
	std::cout << "AlignmentFit:          " << AlignmentFit << " rad" << std::endl;
}
