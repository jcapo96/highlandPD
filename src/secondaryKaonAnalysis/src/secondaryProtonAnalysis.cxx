#include "secondaryKaonAnalysis.hxx"
#include "secondaryProtonAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"

#include "secondaryProtonSelection.hxx"
#include "kaonAnalysisUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

#include "dQdxCalVariation.hxx"
#include "dQdxXCalVariation.hxx"
#include "dQdxYZCalVariation.hxx"
#include "dQdxNormVariation.hxx"
#include "RecombinationVariation.hxx"
#include "ResidualRangeVariation.hxx"
#include "BeamPartIdEffWeight.hxx"
#include "BrokenTrackWeight.hxx"
#include "SCEVariation.hxx"
#include "SCEGeometricVariation.hxx"
#include "NominalBeamMomWeight.hxx"
#include "BeamCompositionWeight.hxx"
#include "ProtonBackgroundWeight.hxx"

#include "HitPitchSCECorrection.hxx"
#include "HitPositionSCECorrection.hxx"
#include "CalorimetryCalibration.hxx"
#include "ParticlePositionSCECorrection.hxx"
#include "DaughterFinderCorrection.hxx"
#include "dEdxMCCorrection.hxx"

//********************************************************************
secondaryProtonAnalysis::secondaryProtonAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

//********************************************************************
bool secondaryProtonAnalysis::Initialize(){
//********************************************************************
  
  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("secondaryKaonAnalysis.SecondaryProtonAnalysis.MinAccumLevelToSave"));
  
  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle
  
  // Add standard categories for the candidates
  anaUtils::AddStandardObjectCategories("candidate"   ,kaonTree::ncandidates,"ncandidates",1);

  return true;
}

//********************************************************************
void secondaryProtonAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void secondaryProtonAnalysis::DefineSelections(){
//********************************************************************
  
  sel().AddSelection("SecondaryProtonSelection", "secondary proton selection",    new secondaryProtonSelection(false));  
}

//********************************************************************
void secondaryProtonAnalysis::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
  corr().AddCorrection(0, "sce geometric correction", new ParticlePositionSCECorrection());
  corr().AddCorrection(1, "dEdx MC correction"      , new dEdxMCCorrection());
}

//********************************************************************
void secondaryProtonAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void secondaryProtonAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void secondaryProtonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
    
  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_EventInfo(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());
  //standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());

  // -------- Add candidates variables ----------------------
  kaonTree::AddKaonVariables_KaonCandidatesReco    (output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesHitsReco(output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesTrue    (output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);

  //TEST
  AddVarMaxSizeVF(output(), candidates_hitvector_dedx    , "all candidates all hits dedx"    , candidates_allhits, 3000);
  AddVarMaxSizeVF(output(), candidates_hitvector_resrange, "all candidates all hits resrange", candidates_allhits, 3000);
  AddVarMaxSizeVF(output(), candidates_hitvector_thetaYZ , "all candidates all hits thetaYZ" , candidates_allhits, 3000);
  AddVarMaxSizeVF(output(), candidates_hitvector_thetaXZ , "all candidates all hits thetaXZ" , candidates_allhits, 3000);
  AddVarMaxSizeVF(output(), candidates_hitvector_tracklength , "all candidates all hits tracklength" , candidates_allhits, 3000);
}

//********************************************************************
void secondaryProtonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void secondaryProtonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 
  
  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_EventInfo(        output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);    
  
  // ---------- kaon candidates variables --------------
  if(box().Candidates.size()>0){
    for(int i = 0; i < std::min((int)box().Candidates.size(),(int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES); i++){
      AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(GetBunch(), box().Candidates[i]->ParentID));
      kaonTree::FillKaonVariables_KaonCandidatesReco    (output(), box().Candidates[i] ,parent);
      //kaonTree::FillKaonVariables_KaonCandidatesHitsReco(output(), box().Candidates[i]);
      kaonTree::FillKaonVariables_KaonCandidatesTrue    (output(), box().Candidates[i]);
      output().IncrementCounter(kaonTree::ncandidates);
    } 
  }
  
  //TEST
  //loop over candidates and add their hits to the tree
  if(box().Candidates.size()>0){
    for(int ipart = 0; ipart < std::min((int)box().Candidates.size(),(int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES); ipart++){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(box().Candidates[ipart]);
      double thetaYZ = atan(abs(part->DirectionStart[1]/part->DirectionStart[2]))*180/TMath::Pi();
      double thetaXZ = atan(abs(part->DirectionStart[0]/part->DirectionStart[2]))*180/TMath::Pi();
      for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	output().FillVectorVar(candidates_hitvector_dedx,part->Hits[2][ihit].dEdx);
	output().FillVectorVar(candidates_hitvector_resrange,part->Hits[2][ihit].ResidualRange);
	output().FillVectorVar(candidates_hitvector_thetaYZ,(Float_t)thetaYZ);
	output().FillVectorVar(candidates_hitvector_thetaXZ,(Float_t)thetaXZ);
	output().FillVectorVar(candidates_hitvector_tracklength,(Float_t)part->Length);
	output().IncrementCounter(candidates_allhits);
      }
    } 
  }
}

//********************************************************************
void secondaryProtonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
void secondaryProtonAnalysis::FillTruthTree(){
//********************************************************************

  // Fill the "truth" tree

  if(!output().HasTree(OutputManager::truth))return;

  output().SetCurrentTree(OutputManager::truth);

  // Loop over all true particles
  std::vector<AnaTrueParticleB*> trueParts = GetSpill().TrueParticles;
  for(std::vector<AnaTrueParticleB*>::iterator it = trueParts.begin(); it!=trueParts.end(); it++) {
    AnaTrueParticlePD& part = *static_cast<AnaTrueParticlePD*>(*it);

    // Check if this particle needs to be saved in the truth tree
    if (!CheckFillTruthTree(&part)) continue;

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

//********************************************************************
bool secondaryProtonAnalysis::CheckFillTruthTree(const AnaTrueParticlePD* part){
//********************************************************************

  if(!part)return false;
  else if(abs(part->PDG)==2212){
    return true;
  }
  else return false;
}

//********************************************************************
void secondaryProtonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  //baseAnalysis::FillTruthTreeBase(vtx);

  // Event variables
  output().FillVar(run,    GetSpill().EventInfo->Run);
  output().FillVar(subrun, GetSpill().EventInfo->SubRun);
  output().FillVar(evt,    GetSpill().EventInfo->Event);
}

//********************************************************************
void secondaryProtonAnalysis::FillCategories(){
//********************************************************************

  // For the beam track
  if(box().MainTrack)
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils
  
  // for the candidates
  if(box().Candidates.size()>0)
    for(int i = 0; i < (int)box().Candidates.size(); i++)
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]),              "candidate",   1);
    
  // For the beam 
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart){
    anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils

  }
}

//********************************************************************
bool secondaryProtonAnalysis::FinalizeConfiguration(){
//********************************************************************
  
  // This is called before FillMicroTrees()

  const Int_t CandidatesAccumLevel = 1; //TODO: think about this
  
  // Save the accum level for the true particle associated to the recon one such that efficiencies can be computed from the truth tree
  if(conf().GetCurrentConfigurationIndex() != ConfigurationManager::default_conf)
    return true;

  std::vector<AnaTrueParticleB*> trueParts = GetSpill().TrueParticles;
  // Loop over all true particles
  for(std::vector<AnaTrueParticleB*>::iterator itp = trueParts.begin(); itp!=trueParts.end(); itp++) {
    AnaTrueParticleB* part = *itp;
    if(!part) continue;
    
    // When the AccumLevel has not been already saved for this vertex 
    if(part->AccumLevel.size() == 0)
      part->AccumLevel.resize(sel().GetNEnabledSelections());
    
    for(std::vector<SelectionBase*>::iterator it = sel().GetSelections().begin(); it != sel().GetSelections().end(); it++){
      
      if(!(*it)->IsEnabled())continue;
      
      Int_t isel = (*it)->GetEnabledIndex();

      if (part->AccumLevel[isel].size() == 0)
        part->AccumLevel[isel].resize((*it)->GetNBranches());
      
      bool found = false;
      int branch = 0;
      int max_accum_level = -1;

      for(int i = 0; i < (int)box().Candidates.size(); i++){
	if(!box().Candidates[i]->TrueObject)continue;
	if(box().Candidates[i]->TrueObject->ID == part->ID){
	  if((*it)->GetAccumCutLevel(i)>max_accum_level){
	    max_accum_level = (*it)->GetAccumCutLevel(i);
	    branch = i;
	  }
	  found = true;
	  //break;
	}
      }

      int min_accum_level = 100;
      for(int i = 0; i < (int)box().Candidates.size(); i++){
	if((*it)->GetAccumCutLevel(i)<min_accum_level)min_accum_level=(*it)->GetAccumCutLevel(i);
      }

      if(!found){
	part->AccumLevel[0][0] = std::min(min_accum_level, CandidatesAccumLevel);
	std::fill(part->AccumLevel[0].begin()+1, part->AccumLevel[0].end(), -1);
      }
      else{
	part->AccumLevel[0][0] = (*it)->GetAccumCutLevel(branch);
	std::fill(part->AccumLevel[0].begin()+1, part->AccumLevel[0].end(), -1);
      }
    }
  }
  return true;
}
