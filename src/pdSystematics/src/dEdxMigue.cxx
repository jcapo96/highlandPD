#include "dEdxMigue.hxx"
#include "DataClasses.hxx"
#include "EventBoxPD.hxx"
#include <cassert>

//#define DEBUG

//********************************************************************
dEdxMigue::dEdxMigue(): EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","BeamPartIdEffWeight", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
  _cal = new Calorimetry();
}

//********************************************************************
void dEdxMigue::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  _cal->Initialize();

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      _cal->CalibrateHit(part->Hits[2][ihit]);
      std::cout << original->Hits[2][ihit].dEdx << " " << part->Hits[2][ihit].dEdx << std::endl;
    }
  }
}

//********************************************************************
bool dEdxMigue::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dEdxMigue::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
//**************************************************

  (void)event;

  
  return false;
}

//********************************************************************
Int_t dEdxMigue::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
//********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kLongTracks;
  }

  return ngroups;
}

