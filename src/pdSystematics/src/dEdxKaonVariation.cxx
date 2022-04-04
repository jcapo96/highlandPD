#include "dEdxKaonVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
dEdxKaonVariation::dEdxKaonVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","dEdxKaon", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void dEdxKaonVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for(Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    // vary the calibrated dEdx for each hit
    for(Int_t iplane = 2; iplane < 3; iplane++){
      for(UInt_t ihit = 0; ihit < part->Hits[iplane].size(); ihit++){

	// We need the errors part of the data file but as well the relative uncertainty for sigma
	Float_t mean_corr, mean_var;
	Int_t mean_index;

	// Get the systematics parameters for this particle type
	if(!GetBinValues(part->Hits[iplane][ihit].ResidualRange, mean_corr, mean_var, mean_index))return;

	//std::cout << part->Hits[iplane][ihit].ResidualRange << " " << mean_corr << " " << mean_var << " " << mean_index << std::endl;
	
	// Apply the variation to the event model quantities
	part->Hits[iplane][ihit].dEdx_calib = original->Hits[iplane][ihit].dEdx_calib * (mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);
      }
      
      // Recompute the top level derived quantities
      std::pair<double, int> chi2pid = pdAnaUtils::Chi2PID(*part, 2212);
      part->Chi2Proton  = chi2pid.first;
      part->Chi2ndf     = chi2pid.second;
    }
  }
}

//********************************************************************
bool dEdxKaonVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for(Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    for(Int_t iplane = 2; iplane < 3; iplane++){
      for(Int_t ihit = 0; ihit < part->Hits[iplane].size(); ihit++){
        part->Hits[iplane][ihit].dEdx_calib = original->Hits[iplane][ihit].dEdx_calib;
      }
    }

    // reset the top level quantities
    part->Chi2Proton      = original->Chi2Proton;
    part->Chi2ndf         = original->Chi2ndf;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dEdxKaonVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;

  //get true particle
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part.TrueObject);
  if(!truePart)return false;

  //systematic only for true kaons
  if(truePart->PDG == 321)return true;
  else return false;
}
