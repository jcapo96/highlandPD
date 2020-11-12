#include "LifetimeVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
LifetimeVariation::LifetimeVariation(): EventVariationBase(),BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","Lifetime", BinnedParams::k1D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
void LifetimeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  double calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //
  //  double norm_factor[3] = {1.0078, 1.0082, 0.9947};
  double norm_factor[3] = {1.0078, 1.0082, 0.9946};
  double fEventT0=500;

  
  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

    if (!original)         continue;
    
    // We need the errors part of the data file but as well the relative uncertainty for sigma
    Float_t mean_corr, mean_var;
    Int_t mean_index;
    
    // Get the systematics parameters
    if (!GetBinValues(0.5, mean_corr,  mean_var,  mean_index))  return;

    fCalo.SetElectronLifetime(mean_corr +  mean_var*toy.GetToyVariations(_index)->Variations[mean_index]);


    
    //    TVector3 dir = anaUtils::ArrayToTVector3(part->DirectionStart);
    //    Int_t plane = 2;
    //    Float_t pitch3D = pdAnaUtils::Compute3DWirePitch(plane, dir);


    std::vector<double> dedx_vector;
    
    // Apply the variation to the evnt model quantities
    for (Int_t i=2;i<3;i++){
      for (Int_t j=0;j<part->Hits[i].size();j++){

        AnaHitPD& hit = part->Hits[i][j];

                
        // dq/dx should be in e/cm
        // dE/dx is returned in MeV/cm
        
        double dQdx_e      = part->Hits[i][j].dQdx/calib_factor[2]*norm_factor[2];
        part->Hits[i][j].dEdx_corr = part->Hits[i][j].dEdx = fCalo.dEdx_from_dQdx_e(dQdx_e,
                                                                                    hit.PeakTime(),
                                                                                    fEventT0);
        dedx_vector.push_back(part->Hits[i][j].dEdx_corr);
      }
    }

    // Recompute the derived quantities. TODO: This should be done somewhere else, in a more transparent way
    std::pair<double, int> chi2pid = pdAnaUtils::Chi2PID(*part, NULL);
    part->Chi2Proton  = chi2pid.first;
    part->Chi2ndf     = chi2pid.second;
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dedx_vector);

  }
}

//********************************************************************
bool LifetimeVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************


  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=0;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        part->Hits[i][j].dEdx      = original->Hits[i][j].dEdx;
        part->Hits[i][j].dEdx_corr = original->Hits[i][j].dEdx_corr;
      }
    }
    part->Chi2Proton      = original->Chi2Proton;
    part->Chi2ndf         = original->Chi2ndf;
    part->truncLibo_dEdx  = original->truncLibo_dEdx;
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool LifetimeVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t LifetimeVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

