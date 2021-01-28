#include "DerivedQuantitiesVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "Parameters.hxx"

//********************************************************************
DerivedQuantitiesVariation::DerivedQuantitiesVariation(): EventVariationBase(){
//********************************************************************

  _recombEnabled = (bool)ND::params().GetParameterI("pionAnalysis.Systematics.EnableRecombination");
}

//********************************************************************
void DerivedQuantitiesVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  
  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)         continue;
        
    std::vector<double> dedx_vector;

    // Loop over hits in the collection plane (plane 2)
    for (UInt_t j=0;j<part->Hits[2].size();j++){
      AnaHitPD& hit = part->Hits[2][j];

      if (!_recombEnabled){
        // recompute the dEdx from the varied dQdx
        // only when the recombination variation is not enabled, since there dEdx is recomputed
        hit.dEdx      = fCalo.dEdx_from_dQdx(hit.dQdx);
        hit.dEdx_corr = fCalo.dEdx_from_dQdx(hit.dQdx_corr);
      }
        
      // Fill the dedx_vector for the collection plane to recompute truncLibo_dEdx    
      dedx_vector.push_back(hit.dEdx_corr);
    }

    // Recompute the top level derived quantities
    std::pair<double, int> chi2pid = pdAnaUtils::Chi2PID(*part, NULL);
    part->Chi2Proton  = chi2pid.first;
    part->Chi2ndf     = chi2pid.second;
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dedx_vector);
  }
}

//********************************************************************
bool DerivedQuantitiesVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    part->Chi2Proton      = original->Chi2Proton;
    part->Chi2ndf         = original->Chi2ndf;
    part->truncLibo_dEdx  = original->truncLibo_dEdx;
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool DerivedQuantitiesVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t DerivedQuantitiesVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

