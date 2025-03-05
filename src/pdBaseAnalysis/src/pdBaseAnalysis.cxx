#include "pdBaseAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
pdBaseAnalysis::pdBaseAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

//********************************************************************
bool pdBaseAnalysis::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  return true;
}

//********************************************************************
void pdBaseAnalysis::DefineInputConverters(){
//********************************************************************

}

//********************************************************************
void pdBaseAnalysis::DefineSelections(){
//********************************************************************

}

//********************************************************************
void pdBaseAnalysis::DefineCorrections(){
//********************************************************************

}

//********************************************************************
void pdBaseAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

}

//********************************************************************
void pdBaseAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

}

//********************************************************************
void pdBaseAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
   
}

//********************************************************************
void pdBaseAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

}

//********************************************************************
void pdBaseAnalysis::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 
  
}

//********************************************************************
void pdBaseAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************
  
}

//********************************************************************
void pdBaseAnalysis::FillTruthTree(){
//********************************************************************

  // Fill the "truth" tree

  if(!output().HasTree(OutputManager::truth))return;

  output().SetCurrentTree(OutputManager::truth);

  // Loop over all true particles
  std::vector<AnaTrueParticleB*> trueParts = GetSpill().TrueParticles;
  for(std::vector<AnaTrueParticleB*>::iterator it = trueParts.begin(); it!=trueParts.end(); it++) {
    AnaTrueParticlePD& part = *static_cast<AnaTrueParticlePD*>(*it);

    // Check if this particle needs to be saved in the truth tree
    if (!CheckFillTruthTreePD(&part)) continue;

    // Initialize the truth tree. We must do that for each saved particle since several entries may correspond to the same spill
    output().InitializeTree(OutputManager::truth);

    // accumulated cut levels to compute efficiencies. This is taken directly from the AnaTrueObject
    Int_t accumLevel=0;  
    for(std::vector<SelectionBase*>::iterator itf = sel().GetSelections().begin(); itf != sel().GetSelections().end(); itf++){
      if(!(*itf)->IsEnabled())continue;

      Int_t isel = (*itf)->GetEnabledIndex();

      for(UInt_t ibranch=0;ibranch<(*itf)->GetNBranches();ibranch++){
        if(part.AccumLevel.size()>0) accumLevel = part.AccumLevel[isel][ibranch];
        if(sel().GetNEnabledSelections()>1){
          if(sel().GetNMaxBranches()>1)
            output().FillMatrixVar(accum_level, accumLevel, isel, ibranch);
          else
            output().FillVectorVar(accum_level, accumLevel, isel);
        }
        else{
          if (sel().GetNMaxBranches()>1)
            output().FillVectorVar(accum_level, accumLevel, ibranch);
          else
            output().FillVar(accum_level, accumLevel);
        }
      }
    }

    // Reset the categories for the current track
    cat().ResetCurrentCategories();

    // Call the derive classes functions, that also fills the categories
    FillTruthTree(part);

    // Fill the truth tree provided the true codes for color drawing
    std::map< std::string, TrackCategoryDefinition* >::iterator its;
    Int_t categ_index = AnalysisAlgorithm::firstCategory;
    for(its=cat().GetCategories().begin();its!=cat().GetCategories().end();its++, categ_index++ ){
      std::string categ_name = its->first;
      TrackCategoryDefinition& categ = *(its->second);
      if(categ.IsMultiType()){
        for(unsigned int i=0;i<categ.GetNTypes();i++)
           output().FillVectorVar(categ_index, (int)cat().CheckCategoryType(categ_name,i),i);
      }
      else if(!categ.IsObject()) output().FillVar(categ_index, cat().GetCode(categ_name));
    }

    // Fill the tree
    output().GetTree(OutputManager::truth)->Fill();
  }
}

// //********************************************************************
// bool pdBaseAnalysis::CheckFillTruthTreePD(const AnaTrueParticlePD* part){
// //********************************************************************

//   return true;
// }

//********************************************************************
void pdBaseAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************

  // Event variables
  output().FillVar(run,    GetSpill().EventInfo->Run);
  output().FillVar(subrun, GetSpill().EventInfo->SubRun);
  output().FillVar(evt,    GetSpill().EventInfo->Event);

}

//********************************************************************
void pdBaseAnalysis::FillCategories(){
//********************************************************************

}

//********************************************************************
bool pdBaseAnalysis::FinalizeConfiguration(){
//********************************************************************
  
  // This is called before FillMicroTrees()

  // Save the accum level for the true vertex associated to the recon one such that efficiencies can be computed from the truth tree
  if (conf().GetCurrentConfigurationIndex() != ConfigurationManager::default_conf)
    return true;

  std::vector<AnaTrueParticlePD*> trueParticles;
  trueParticles.reserve(NMAXTRUEVERTICES); //to be changed
  std::vector<AnaTrueParticlePD*>::const_iterator iter;


  AnaTrueParticlePD* trueParticle = NULL;
  if (GetParticle()) trueParticle = static_cast<AnaTrueParticlePD*>(GetParticle()->TrueObject);
  else             trueParticle = static_cast<AnaTrueParticlePD*>(GetTrueParticle());

  if (trueParticle) trueParticles.push_back(trueParticle);

  // If true particle does not exist (e.g. can happen that reco particle is not yet available at this step of the selection) then 
  // store the accum_level to all true vertices of the bunch -> i.e. this basically corresponds to the fact that event as a whole 
  // passed some cuts
  if (trueParticles.size() == 0){
    // Loop over all true vertices in the event and found those that belong to the bunch being processed
    std::vector<AnaTrueParticleB*> particles = GetSpill().TrueParticles;
    for (std::vector<AnaTrueParticleB*>::iterator it = particles.begin(); it!=particles.end(); it++) {
      if (!(*it)) continue;
      AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(*it);

      // not sure how this works in pd, for the moment we store all
      // // Check the bunch
      // if (GetBunch().Bunch != vtx->Bunch) continue;

      trueParticles.push_back(truepart); 
    }
  }

  // Loop over all true vertices of interest
  for (iter = trueParticles.begin(); iter != trueParticles.end(); iter++){ 
    AnaTrueParticlePD* truepart = *iter;
    if (!truepart) continue;

    // When the AccumLevel has not been already saved for this particle 
    if (truepart->AccumLevel.size() == 0)
      truepart->AccumLevel.resize(sel().GetNEnabledSelections());
    
    for (std::vector<SelectionBase*>::iterator it = sel().GetSelections().begin(); it != sel().GetSelections().end(); it++){

      if (!(*it)->IsEnabled()) continue;

      Int_t isel = (*it)->GetEnabledIndex();

      if (truepart->AccumLevel[isel].size() == 0)
        truepart->AccumLevel[isel].resize((*it)->GetNBranches());

      for (UInt_t ibranch=0;ibranch<(*it)->GetNBranches();ibranch++)
        truepart->AccumLevel[isel][ibranch]=(*it)->GetAccumCutLevel(ibranch);
    }
  }
  return true;
}
