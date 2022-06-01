#include "ResidualRangeVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
ResidualRangeVariation::ResidualRangeVariation():EventVariationBase(100),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","ResidualRange", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  //This systematic uses a Uniform random throw
  _PDF = kUniform;
}

//********************************************************************
void ResidualRangeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
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

    //Get the systematic source values
    Float_t width;
    GetSigmaValueForBin(itrk, width);
    //std::cout << "BIN " << width << std::endl;

    // compute the maximum variation as a function of the end angle of the track with respect to the y axis
    double vy = abs(part->DirectionEnd[1]);
    // if end angle has a bad value (-0.9999 or 0.9999) try with start angle
    if(vy > 0.999)vy = abs(part->DirectionStart[1]);
    // if this doesn't work either, skip the track
    if(vy > 0.999)continue;

    //get the maximum possible variation given the angle with respect to the y axis
    double max_var = width / sqrt(1 - pow(vy,2));
    //std::cout << vy << " " << max_var << std::endl; 
    //std::cout << "variation " << toy.GetToyVariations(_index)->Variations[itrk] << std::endl;

    //reescale the variation to go from -max_var/2 to max_var/2
    double var = toy.GetToyVariations(_index)->Variations[itrk] * max_var - max_var/2;
    //std::cout << "reescaled variation " << var std::endl;

    // vary the residual range for each hit
    for(Int_t iplane = 2; iplane < 3; iplane++){
      for(Int_t ihit = 0; ihit < (int)part->Hits[iplane].size(); ihit++){
	part->Hits[iplane][ihit].ResidualRange = original->Hits[iplane][ihit].ResidualRange + var;
      }
    }
    // Recompute the top level derived quantities
    std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    part->Chi2Proton  = chi2pidprot.first;
    part->Chi2ndf     = chi2pidprot.second;
    std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    part->Chi2Muon    = chi2pidmuon.first;
  }
}

//********************************************************************
bool ResidualRangeVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for(Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    for(Int_t iplane = 2; iplane < 3; iplane++){
      for(Int_t ihit = 0; ihit < (int)part->Hits[iplane].size(); ihit++){
        part->Hits[iplane][ihit].ResidualRange = original->Hits[iplane][ihit].ResidualRange;
      }
    }

    // reset the top level quantities
    part->Chi2Proton = original->Chi2Proton;
    part->Chi2ndf    = original->Chi2ndf;
    part->Chi2Muon   = original->Chi2Muon;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool ResidualRangeVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
