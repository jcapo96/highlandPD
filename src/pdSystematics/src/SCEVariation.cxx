#include "SCEVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCEVariation::SCEVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  //check we have enough bins to run this systematic, otherwise exit
  if(GetNBins() > (int)NMAXSYSTSOURCEBINS){
    std::cout << "SCEVariation::SCEVariation()." << std::endl;
    std::cout << "the source files has more bins than allowed for the systematic propagation!" << std::endl;
    std::cout << "N MAX BINS allowed " << NMAXSYSTSOURCEBINS << std::endl;
    std::cout << "N BINS on source   " << GetNBins() << std::endl;
    std::cout << "if you need more bins, look at psyckeUtils/src/BinnedParams.hxx" << std::endl;
    std::exit(1);
  }
  
  //Initialize SCE object to apply variations
  for(int i = 0; i < 100; i++)_sce[i] = NULL;

  // Initialize Calorimetry object pointing to null SCE, it will be replaced afterwards
  _cal = new Calorimetry();
  _cal->Initialize(NULL);
}

//********************************************************************
SCEVariation::~SCEVariation(){
//********************************************************************

  for(int i = 0; i < 100; i++)delete _sce[i];
  delete _cal;
}

//********************************************************************
void SCEVariation::Initialize(){
//********************************************************************
  
  std::cout << "SCEVariation::Initialize(). Initializing toy SCE" << std::endl;
  for(int i = 0; i < 100; i++){
    _sce[i] = new SpaceCharge();
    _sce[i]->Initialize();
  }
}

//********************************************************************
void SCEVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);
  const AnaEventPD* eventpd = static_cast<const AnaEventPD*>(&event);
  
  //initialize the SCE for this toy if it has not been initialized yet 
  int toy_index = toy.GetToyIndex();
  _cal->SetSCE(_sce[toy_index],false);

  //Vary SCE map if it hasn't been varied yet
  if(!_sce[toy_index]->IsVaried())
    VarySCEMapGlobally(toy);

  //set appropiate lifetime depending on mc/data and run
  _cal->SetLifetime(*eventpd);
  
  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
     AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
  
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    if(part->Hits[2].empty())continue;
    
    //important, two things are varied: XYZ position of the hit and pitch
    //loop here to avoid looping more than once
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      _sce[toy_index]->ApplyPositionCorrection(part->Hits[2][ihit]);
      _cal->ApplySCECorrection(part->Hits[2][ihit]);
      _cal->ApplyLifetimeCorrection(part->Hits[2][ihit]);
    }
  }
}

//********************************************************************
bool SCEVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  //_sce->ResetToNominal();
  
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].dQdx = original->Hits[2][ihit].dQdx;
      part->Hits[2][ihit].dQdx_elife = original->Hits[2][ihit].dQdx_elife;
      part->Hits[2][ihit].Position.SetX(original->Hits[2][ihit].Position.X());
      part->Hits[2][ihit].Position.SetY(original->Hits[2][ihit].Position.Y());
      part->Hits[2][ihit].Position.SetZ(original->Hits[2][ihit].Position.Z());
    }
  }
  
  // Don't reset the spill to corrected
  return false;
}

//********************************************************************
void SCEVariation::VarySCEMapGlobally(const ToyExperiment& toy){
//********************************************************************
  
  int toy_index = toy.GetToyIndex();

  if(toy_index==0)std::cout << "SCEVariation::VarySCEMap(). Varying toy maps" << std::endl;
  if(toy_index==99)std::cout << "SCEVariation::VarySCEMap(). Toy maps varied" << std::endl;

  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin
  _sce[toy_index]->ApplyGlobalVariation(1 + sigma*toy.GetToyVariations(_index)->Variations[0]);
}

//********************************************************************
void SCEVariation::VarySCEMapLocally(const ToyExperiment& toy){
//********************************************************************
  
  int toy_index = toy.GetToyIndex();

  if(toy_index==0)std::cout << "SCEVariation::VarySCEMap(). Varying toy maps" << std::endl;
  if(toy_index==99)std::cout << "SCEVariation::VarySCEMap(). Toy maps Varied" << std::endl;

  //get map binning
  int nbinsx = _sce[toy_index]->GetNbinsX();
  int nbinsy = _sce[toy_index]->GetNbinsY();
  int nbinsz = _sce[toy_index]->GetNbinsZ();

  //loop over voxels
  Float_t xcenter,ycenter,zcenter;
  for(int ix = 0; ix < nbinsx; ix++){
    xcenter = _sce[toy_index]->GetBinCenterX(ix+1);
    for(int iy = 0; iy < nbinsy; iy++){
      ycenter = _sce[toy_index]->GetBinCenterY(iy+1);
      for(int iz = 0; iz < nbinsz; iz++){
	zcenter = _sce[toy_index]->GetBinCenterZ(iz+1);
	Int_t index = 0;
	Float_t sigma = 0;
	if(!GetBinSigmaValue(xcenter,ycenter,zcenter,sigma,index))
	  continue;
	_sce[toy_index]->ApplyVoxelVariation(ix+1,iy+1,iz+1,
					     1 + sigma*toy.GetToyVariations(_index)->Variations[index],
					     false); //do not reset splines now, only at the end
      }
    }
  }
  _sce[toy_index]->ResetSplines(); //now recompute splines;
  
}

//**************************************************
bool SCEVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
