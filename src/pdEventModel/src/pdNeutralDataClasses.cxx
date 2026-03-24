#define pdNeutralDataClasses_C

#include "pdNeutralDataClasses.hxx"
#include "AnalysisUtils.hxx"

// define a constant value for uninitialised parameters
const Double_t kDoubleUnassigned = -999.;
const Int_t    kIntUnassigned = -999;
const Float_t  kFloatUnassigned = -999.;

//********************************************************************
AnaTrueEquivalentVertexPD::AnaTrueEquivalentVertexPD(){
//********************************************************************

  TrueParticles.clear();
  OriginalDistance = kFloatUnassigned;
  MinimumDistance = kFloatUnassigned;
  OpeningAngle = kFloatUnassigned;
  Position[0] = kFloatUnassigned;
  Position[1] = kFloatUnassigned;
  Position[2] = kFloatUnassigned;
  Direction[0] = kFloatUnassigned;
  Direction[1] = kFloatUnassigned;
  Direction[2] = kFloatUnassigned;
  DirectionFit[0] = kFloatUnassigned;
  DirectionFit[1] = kFloatUnassigned;
  DirectionFit[2] = kFloatUnassigned;
  PositionPandora[0] = kFloatUnassigned;
  PositionPandora[1] = kFloatUnassigned;
  PositionPandora[2] = kFloatUnassigned;
  PositionFit[0] = kFloatUnassigned;
  PositionFit[1] = kFloatUnassigned;
  PositionFit[2] = kFloatUnassigned;
  DirectionFit[0] = kFloatUnassigned;
  DirectionFit[1] = kFloatUnassigned;
  DirectionFit[2] = kFloatUnassigned;
}

//********************************************************************
AnaTrueEquivalentVertexPD::~AnaTrueEquivalentVertexPD(){
//********************************************************************

}

//********************************************************************
AnaTrueEquivalentVertexPD::AnaTrueEquivalentVertexPD(const AnaTrueEquivalentVertexPD& vertex){
//********************************************************************

  TrueParticles = vertex.TrueParticles;
  OriginalDistance = vertex.OriginalDistance;
  MinimumDistance = vertex.MinimumDistance;
  OpeningAngle = vertex.OpeningAngle;
  Position[0] = vertex.Position[0];
  Position[1] = vertex.Position[1];
  Position[2] = vertex.Position[2];
  Direction[0] = vertex.Direction[0];
  Direction[1] = vertex.Direction[1];
  Direction[2] = vertex.Direction[2];
  DirectionFit[0] = vertex.DirectionFit[0];
  DirectionFit[1] = vertex.DirectionFit[1];
  DirectionFit[2] = vertex.DirectionFit[2];
  PositionPandora[0] = vertex.PositionPandora[0];
  PositionPandora[1] = vertex.PositionPandora[1];
  PositionPandora[2] = vertex.PositionPandora[2];
  PositionFit[0] = vertex.PositionFit[0];
  PositionFit[1] = vertex.PositionFit[1];
  PositionFit[2] = vertex.PositionFit[2];
  DirectionFit[0] = vertex.DirectionFit[0];
  DirectionFit[1] = vertex.DirectionFit[1];
  DirectionFit[2] = vertex.DirectionFit[2];
}

//********************************************************************
void AnaTrueEquivalentVertexPD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrueEquivalentVertexPD --------- " << std::endl;
  std::cout << "TrueParticles size:    " << TrueParticles.size() << std::endl;
  std::cout << "OriginalDistance:      " << OriginalDistance << " cm" << std::endl;
  std::cout << "MinimumDistance:       " << MinimumDistance << " cm" << std::endl;
  std::cout << "OpeningAngle:          " << OpeningAngle << " degrees" << std::endl;
  std::cout << "Position:              " << Position[0] << " " << Position[1] << " " << Position[2] << std::endl;
  std::cout << "Direction:             " << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
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
AnaAnnihilationVertexPD::AnaAnnihilationVertexPD() : AnaVertexPD(){
//********************************************************************
  // All initialization handled by base class AnaVertexPD
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
  TrueEquivalentNeutralParticle = NULL;
  ImpactParameter = kFloatUnassigned;
  Mass = kFloatUnassigned;
  Momentum = kFloatUnassigned;
  PDG = kIntUnassigned;
  Lifetime = kFloatUnassigned;
  DecayLength = kFloatUnassigned;
  NRecoHitsInVertex = kIntUnassigned;
  NeutralScore = kFloatUnassigned;
  HitsAlignment = kFloatUnassigned;
  NHitsInCylinder = kIntUnassigned;
  HitsAvgDistance = kFloatUnassigned;
  HitsRMSDistance = kFloatUnassigned;
  HitsLongitudinalSpan = kFloatUnassigned;
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
  TrueEquivalentNeutralParticle = neutralParticle.TrueEquivalentNeutralParticle;
  RecoParticle = neutralParticle.RecoParticle;
  ImpactParameter = neutralParticle.ImpactParameter;
  Mass = neutralParticle.Mass;
  Momentum = neutralParticle.Momentum;
  PDG = neutralParticle.PDG;
  Lifetime = neutralParticle.Lifetime;
  DecayLength = neutralParticle.DecayLength;
  NRecoHitsInVertex = neutralParticle.NRecoHitsInVertex;
  NeutralScore = neutralParticle.NeutralScore;
  HitsAlignment = neutralParticle.HitsAlignment;
  NHitsInCylinder = neutralParticle.NHitsInCylinder;
  HitsAvgDistance = neutralParticle.HitsAvgDistance;
  HitsRMSDistance = neutralParticle.HitsRMSDistance;
  HitsLongitudinalSpan = neutralParticle.HitsLongitudinalSpan;
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
  std::cout << "TrueEquivalentNeutralParticle:   " << (TrueEquivalentNeutralParticle ? "Yes" : "No") << std::endl;
  std::cout << "RecoParticle:          " << (RecoParticle ? "Yes" : "No") << std::endl;
  std::cout << "ImpactParameter:       " << ImpactParameter << " cm" << std::endl;
  std::cout << "Mass:                  " << Mass << " GeV/c²" << std::endl;
  std::cout << "Momentum:              " << Momentum << " GeV/c" << std::endl;
  std::cout << "PDG:                   " << PDG << std::endl;
  std::cout << "Lifetime:              " << Lifetime << " ns" << std::endl;
  std::cout << "DecayLength:           " << DecayLength << " cm" << std::endl;
  std::cout << "NRecoHitsInVertex:     " << NRecoHitsInVertex << std::endl;
  std::cout << "RecoParticle:          " << (RecoParticle ? "Yes" : "No") << std::endl;
}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::AnaTrueEquivalentNeutralParticlePD(){
//********************************************************************

  TrueEquivalentVertex = NULL;
  TrueParent = NULL;
  Position[0] = kFloatUnassigned;
  Position[1] = kFloatUnassigned;
  Position[2] = kFloatUnassigned;
  Direction[0] = kFloatUnassigned;
  Direction[1] = kFloatUnassigned;
  Direction[2] = kFloatUnassigned;
  PositionEnd[0] = kFloatUnassigned;
  PositionEnd[1] = kFloatUnassigned;
  PositionEnd[2] = kFloatUnassigned;
  DirectionEnd[0] = kFloatUnassigned;
  DirectionEnd[1] = kFloatUnassigned;
  DirectionEnd[2] = kFloatUnassigned;
  Length = kFloatUnassigned;
  Momentum = kFloatUnassigned;
  MomentumEnd = kFloatUnassigned;
  PDG = kIntUnassigned;
  Generation = kIntUnassigned;
  Process = kIntUnassigned;
  Mass = kFloatUnassigned;
}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::~AnaTrueEquivalentNeutralParticlePD(){
//********************************************************************

}

//********************************************************************
AnaTrueEquivalentNeutralParticlePD::AnaTrueEquivalentNeutralParticlePD(const AnaTrueEquivalentNeutralParticlePD& trueEquivalentNeutralPart){
//********************************************************************

  TrueEquivalentVertex = trueEquivalentNeutralPart.TrueEquivalentVertex;
  TrueParent = trueEquivalentNeutralPart.TrueParent;
  Position[0] = trueEquivalentNeutralPart.Position[0];
  Position[1] = trueEquivalentNeutralPart.Position[1];
  Position[2] = trueEquivalentNeutralPart.Position[2];
  Direction[0] = trueEquivalentNeutralPart.Direction[0];
  Direction[1] = trueEquivalentNeutralPart.Direction[1];
  Direction[2] = trueEquivalentNeutralPart.Direction[2];
  PositionEnd[0] = trueEquivalentNeutralPart.PositionEnd[0];
  PositionEnd[1] = trueEquivalentNeutralPart.PositionEnd[1];
  PositionEnd[2] = trueEquivalentNeutralPart.PositionEnd[2];
  DirectionEnd[0] = trueEquivalentNeutralPart.DirectionEnd[0];
  DirectionEnd[1] = trueEquivalentNeutralPart.DirectionEnd[1];
  DirectionEnd[2] = trueEquivalentNeutralPart.DirectionEnd[2];
  Length = trueEquivalentNeutralPart.Length;
  Momentum = trueEquivalentNeutralPart.Momentum;
  MomentumEnd = trueEquivalentNeutralPart.MomentumEnd;
  PDG = trueEquivalentNeutralPart.PDG;
  Generation = trueEquivalentNeutralPart.Generation;
  Process = trueEquivalentNeutralPart.Process;
  Mass = trueEquivalentNeutralPart.Mass;
}

//********************************************************************
void AnaTrueEquivalentNeutralParticlePD::Print() const{
//********************************************************************

  std::cout << "-------- AnaTrueNeutralParticlePD --------- " << std::endl;
  std::cout << "TrueEquivalentVertex:   " << (TrueEquivalentVertex ? "Set" : "NULL") << std::endl;
  std::cout << "TrueParent:             " << (TrueParent ? "Set" : "NULL") << std::endl;
  std::cout << "Position:               " << Position[0] << " " << Position[1] << " " << Position[2] << std::endl;
  std::cout << "Direction:              " << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
  std::cout << "PositionEnd:            " << PositionEnd[0] << " " << PositionEnd[1] << " " << PositionEnd[2] << std::endl;
  std::cout << "DirectionEnd:           " << DirectionEnd[0] << " " << DirectionEnd[1] << " " << DirectionEnd[2] << std::endl;
  std::cout << "Length:                 " << Length << " cm" << std::endl;
  std::cout << "Momentum:               " << Momentum << " GeV/c" << std::endl;
  std::cout << "PDG:                    " << PDG << std::endl;
  std::cout << "Generation:             " << Generation << std::endl;
  std::cout << "Process:                " << Process << std::endl;
  std::cout << "Mass:                   " << Mass << " GeV/c²" << std::endl;
}

