#include "SCEGeometricVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCEGeometricVariation::SCEGeometricVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  //Initialize SCE object to apply variations
  for(int i = 0; i < 100; i++)_sce[i] = NULL;
}

//********************************************************************
SCEGeometricVariation::~SCEGeometricVariation(){
//********************************************************************
  
  for(int i = 0; i < 100; i++)delete _sce[i];
}

//********************************************************************
void SCEGeometricVariation::Initialize(){
//********************************************************************
  
  std::cout << "SCEGeometricVariation::Initialize(). Initializing toy SCE" << std::endl;
  for(int i = 0; i < 100; i++){
    _sce[i] = new SpaceCharge();
    _sce[i]->Initialize();
  }
}

//********************************************************************
void SCEGeometricVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //initialize the SCE for this toy if it has not been initialized yet 
  int toy_index = toy.GetToyIndex();
  if(!_sce[toy_index]){
    _sce[toy_index] = new SpaceCharge();
    _sce[toy_index]->Initialize();
  }
  
  //Vary SCE map
  if(!_sce[toy_index]->IsVaried())
    VarySCEMap(toy);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
     AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
  
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    
    // Apply variation to TrajectoryPoints
    int ntps = part->TrjPoints.size();
    if(ntps > 0){
      _sce[toy_index]->ApplyTrjPointPositionCorrection(part);
      _sce[toy_index]->ApplyTrjPointDirectionCorrection(part);
     
      //recompute length
      part->Length = pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(part);

      // Apply variation to Position/Direction start/end
      pdAnaUtils::ComputeParticlePositionAndDirection(part);
    }
    
    //apply variation to hits
    //We don't have yet the hit-trjp relationship, this is redundant
    //but there is no other way so far
    if(!part->Hits[2].empty()){
      _sce[toy_index]->ApplyPositionCorrection(part);
      pdAnaUtils::ComputeResidualRange(part);
    }
    std::pair<double,int>result = pdAnaUtils::Chi2PID(*part,2212);
    part->Chi2Proton = result.first;
    result = pdAnaUtils::Chi2PID(*part,13);
    part->Chi2Muon = result.first;
    part->Chi2ndf = result.second;
  }
}

//********************************************************************
bool SCEGeometricVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  //_sce->ResetToNominal();
  
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over trajectorypoints
    for(int itrp = 0; itrp < (int)part->TrjPoints.size(); itrp++){
      part->TrjPoints[itrp].Position = original->TrjPoints[itrp].Position;
      part->TrjPoints[itrp].Direction = original->TrjPoints[itrp].Direction;
    }

    //set lenght and positions/direction
    part->Length = original->Length;
    part->PositionStart[0] = original->PositionStart[0];
    part->PositionStart[1] = original->PositionStart[1];
    part->PositionStart[2] = original->PositionStart[2];
    part->PositionEnd[0] = original->PositionEnd[0];
    part->PositionEnd[1] = original->PositionEnd[1];
    part->PositionEnd[2] = original->PositionEnd[2];
    part->DirectionStart[0] = original->DirectionStart[0];
    part->DirectionStart[1] = original->DirectionStart[1];
    part->DirectionStart[2] = original->DirectionStart[2];
    part->DirectionEnd[0] = original->DirectionEnd[0];
    part->DirectionEnd[1] = original->DirectionEnd[1];
    part->DirectionEnd[2] = original->DirectionEnd[2];

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].ResidualRange = original->Hits[2][ihit].ResidualRange;
      part->Hits[2][ihit].Position.SetX(original->Hits[2][ihit].Position.X());
      part->Hits[2][ihit].Position.SetY(original->Hits[2][ihit].Position.Y());
      part->Hits[2][ihit].Position.SetZ(original->Hits[2][ihit].Position.Z());
    }

    part->Chi2Proton = part->Chi2Proton;
    part->Chi2Muon   = part->Chi2Muon;
    part->Chi2ndf    = part->Chi2ndf;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//********************************************************************
void SCEGeometricVariation::VarySCEMap(const ToyExperiment& toy){
//********************************************************************
  
  int toy_index = toy.GetToyIndex();

  if(toy_index==0)std::cout << "SCEGeometricVariation::VarySCEMap(). Varying toy maps" << std::endl;
  if(toy_index==99)std::cout << "SCEGeometricVariation::VarySCEMap(). Toy maps varied" << std::endl;

  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin
  _sce[toy_index]->ApplyGlobalVariation(1 + sigma*toy.GetToyVariations(_index)->Variations[0]);    
}

//**************************************************
bool SCEGeometricVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
