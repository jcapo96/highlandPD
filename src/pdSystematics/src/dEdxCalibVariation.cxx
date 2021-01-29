#include "dEdxCalibVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include "HitVariationBase.hxx"

/*
  This is the dEdx calibration variation, applying different variations at the hit level 
  (this avoid redoing the it loop in several places). Those has to be added as daughters 
  of this variation in the analysis algorithm. Some options are:
  - Lifetime
  - dQdxCalib
  - Recombination

  This variation also recomputes all derived quantities depending on dQdx:
  - dEdx
  - Chi2Proton
  - Truncated dEdx
*/


//********************************************************************
void dEdxCalibVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Check whether dEdx has to be recomputed or not. Cannot be done in constructor since it depends on the curent configuration 
  bool recompute_dEdx = true;
  for (UInt_t v=0;v<_daughters.size();v++){
    if (_daughters[v]->IsEnabled()){
      if (static_cast<HitVariationBase*>(_daughters[v])->dEdx_recomputed) recompute_dEdx = false;
    }
  }
    
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
    
    // Apply the variation to the event model quantities
    // Only plane 2 for the moment
    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){
        AnaHitPD& hit = part->Hits[i][j];                

        // Apply daughter hit variations        
        for (UInt_t v=0;v<_daughters.size();v++){
          if (_daughters[v]->IsEnabled()){
            static_cast<HitVariationBase*>(_daughters[v])->Apply(toy,hit);
          }
        }

        if (recompute_dEdx){          
          // recompute the dEdx from the varied dQdx
          // only when the recombination variation is not enabled, since there dEdx is recomputed
          hit.dEdx      = fCalo.dEdx_from_dQdx(hit.dQdx);
          //          hit.dEdx_calib = fCalo.dEdx_from_dQdx(hit.dQdx_corr);
        }

        // Fill the dedx_vector for the collection plane to recompute truncLibo_dEdx    
        if (i==2){
          dedx_vector.push_back(hit.dEdx);  // TODO: To be substituted by dEdx_corr      
        }
      }
    }
    // Recompute the top level derived quantities
    std::pair<double, int> chi2pid = pdAnaUtils::Chi2PID(*part, NULL);
    part->Chi2Proton  = chi2pid.first;
    part->Chi2ndf     = chi2pid.second;
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dedx_vector);
  }
}

//********************************************************************
bool dEdxCalibVariation::UndoSystematic(AnaEventC& event){
  //********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    for (Int_t i=2;i<3;i++){
      for (UInt_t j=0;j<part->Hits[i].size();j++){

        // Undo daughter hit variations
        for (UInt_t v=0;v<_daughters.size();v++){
          if (_daughters[v]->IsEnabled()){
            static_cast<HitVariationBase*>(_daughters[v])->UndoSystematic(original->Hits[i][j],part->Hits[i][j]);
          }
        }

        // Reset the dEdx
        part->Hits[i][j].dEdx       = original->Hits[i][j].dEdx;
        part->Hits[i][j].dEdx_calib = original->Hits[i][j].dEdx_calib;
      }
    }
    // reset the top level quantities
    part->Chi2Proton      = original->Chi2Proton;
    part->Chi2ndf         = original->Chi2ndf;
    part->truncLibo_dEdx  = original->truncLibo_dEdx;
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dEdxCalibVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
  //**************************************************

  (void)event;
  return true;
}
