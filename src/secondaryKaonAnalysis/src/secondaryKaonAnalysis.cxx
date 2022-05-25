#include "secondaryKaonAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"

#include "secondaryKaonSelection.hxx"
#include "secondaryKaonXSSelection.hxx"
#include "kaonAnalysisUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

#include "dEdxKaonVariation.hxx"
#include "ResidualRangeVariation.hxx"
#include "BeamPartIdEffWeight.hxx"

//********************************************************************
secondaryKaonAnalysis::secondaryKaonAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

//********************************************************************
bool secondaryKaonAnalysis::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("secondaryKaonAnalysis.MinAccumLevelToSave"));

  // Parameters
  _xs_selection = ND::params().GetParameterI("secondaryKaonAnalysis.XSSelection");
  if(_xs_selection)_selection_name = "secondaryKaonXSSelection";
  else             _selection_name = "secondaryKaonSelection";

  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle
  anaUtils::AddStandardCategories("bestcandidate");  // This is for the best candidate of each event
  anaUtils::AddStandardCategories("bestcandidatedau");  // This is for the best candidate of each event
  
  // Add standard categories for the candidates
  anaUtils::AddStandardObjectCategories("candidate"   ,kaonTree::ncandidates,"ncandidates",1);
  anaUtils::AddStandardObjectCategories("candidatedau",kaonTree::ncandidates,"ncandidates",1);

  // Add our own categories (in secondaryKaonAnalysisUtils.cxx)
  kaonAnaUtils::AddCustomCategories();

  return true;
}

//********************************************************************
void secondaryKaonAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void secondaryKaonAnalysis::DefineSelections(){
//********************************************************************

  if(_xs_selection)sel().AddSelection(_selection_name.c_str(), "kaon XS selection", new secondaryKaonXSSelection(false));
  else             sel().AddSelection(_selection_name.c_str(), "kaon selection",    new secondaryKaonSelection(false));
  
}

//********************************************************************
void secondaryKaonAnalysis::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
}

//********************************************************************
void secondaryKaonAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  evar().AddEventVariation(0,"residual range variation",                new ResidualRangeVariation());
  eweight().AddEventWeight(0,"beam particle identification efficiency", new BeamPartIdEffWeight());
}

//********************************************************************
void secondaryKaonAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if(_enableAllSystConfig){
    if(ND::params().GetParameterI("secondaryKaonAnalysis.Systematics.EnableResidualRangeVar"))conf().EnableEventVariation(0,all_syst);
    if(ND::params().GetParameterI("secondaryKaonAnalysis.Systematics.EnableBeamPartIdEff"))conf().EnableEventWeight(0,all_syst);
  }
}

//********************************************************************
void secondaryKaonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
  standardPDTree::AddStandardVariables_CandidateTrue(output());
  standardPDTree::AddStandardVariables_CandidateReco(output());
  standardPDTree::AddStandardVariables_CandidateHitsReco(output());

  // -------- Add candidates variables ----------------------
  kaonTree::AddKaonVariables_KaonCandidatesReco    (output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesHitsReco(output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesTrue    (output(),secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);

  // -------- Add best candidate variables ----------------------
  kaonTree::AddKaonVariables_KaonBestCandidateReco    (output());
  kaonTree::AddKaonVariables_KaonBestCandidateHitsReco(output());
  kaonTree::AddKaonVariables_KaonBestCandidateTrue    (output());

  // -------- Add toy variables ---------------------------------
  AddToyVarVF(output(), kaonTree::bestcandidate_hit_resrange_toy, "bestcandidate hit residual range", NMAXHITSPERPLANE);
}

//********************************************************************
void secondaryKaonAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for ProtoDUNE analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_BeamTrue(output());
  
  // Add specific variables for the kaon candidates
  kaonTree::AddKaonVariables_TrueKaonCandidates(output());
}

//********************************************************************
void secondaryKaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_BeamReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_CandidateReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_CandidateTrue(    output(), box().MainTrack);    
  standardPDTree::FillStandardVariables_CandidateHitsReco(output(), box().MainTrack);
  
  // ---------- kaon candidates variables --------------
  if(box().Candidates.size()>0){
    for(int i = 0; i < std::min((int)box().Candidates.size(),(int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES); i++){
      AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(GetBunch(), box().Candidates[i]->ParentID));
      kaonTree::FillKaonVariables_KaonCandidatesReco    (output(), box().Candidates[i] ,parent);
      kaonTree::FillKaonVariables_KaonCandidatesHitsReco(output(), box().Candidates[i]);
      kaonTree::FillKaonVariables_KaonCandidatesTrue    (output(), box().Candidates[i]);
      output().IncrementCounter(kaonTree::ncandidates);
    } 
  }
  
  // ---------- best kaon candidate variables --------------
  //get the kaon selection and get the branch with a larger AccumLevel
  SelectionBase* ksel = sel().GetSelection(_selection_name.c_str()); //todo fix this
  int branchmax = 0;
  int max       = 0;  
  for(UInt_t ibranch = 0; ibranch < ksel->GetNBranches(); ibranch++){
    if(ksel->GetAccumCutLevel(ibranch) > max){
      max = ksel->GetAccumCutLevel(ibranch);
      branchmax = ibranch;
    }
  }

  if(box().Candidates.size()>0){
    AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(GetBunch(), box().Candidates[branchmax]->ParentID));
    kaonTree::FillKaonVariables_KaonBestCandidateReco    (output(), box().Candidates[branchmax], parent);
    kaonTree::FillKaonVariables_KaonBestCandidateHitsReco(output(), box().Candidates[branchmax]);
    kaonTree::FillKaonVariables_KaonBestCandidateTrue    (output(), box().Candidates[branchmax]);
  }  
}

//********************************************************************
void secondaryKaonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

  // ---------- best kaon candidate variables --------------
  //get the kaon selection and get the branch with a larger AccumLevel
  SelectionBase* ksel = sel().GetSelection(_selection_name.c_str()); //todo fix this
  int branchmax = 0;
  int max       = 0;  
  for(UInt_t ibranch = 0; ibranch < ksel->GetNBranches(); ibranch++){
    if(ksel->GetAccumCutLevel(ibranch) > max){
      max = ksel->GetAccumCutLevel(ibranch);
      branchmax = ibranch;
    }
  }
 
  if(box().Candidates.size() == 0)return;

  AnaParticlePD* best = box().Candidates[branchmax];
  if(best->Hits[2].empty()){
    for(int j = 0; j < NMAXHITSPERPLANE; j++){
      output().FillToyVectorVar(kaonTree::bestcandidate_hit_resrange_toy,(Float_t)-999.,j);
    }
  }
  else{
    int jmin = std::max<int>(0,(int)best->Hits[2].size()-NMAXHITSPERPLANE);
    int jmax = (int)best->Hits[2].size(); 
    
    std::cout << "TOY" << std::endl;
    for(int j = jmin; j < jmax; j++){
      //std::cout << j-jmin << " " << best->Hits[2][j].ResidualRange << std::endl;
      output().FillToyVectorVar(kaonTree::bestcandidate_hit_resrange_toy,best->Hits[2][j].ResidualRange,j-jmin);
    }
  }
}

//********************************************************************
void secondaryKaonAnalysis::FillTruthTree(){
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
bool secondaryKaonAnalysis::CheckFillTruthTree(const AnaTrueParticlePD* part){
//********************************************************************

  if(!part)return false;
  else if(abs(part->PDG)==321){
    return true;
  }
  else return false;
}

//********************************************************************
void secondaryKaonAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  //baseAnalysis::FillTruthTreeBase(vtx);

  // Event variables
  output().FillVar(run,    GetSpill().EventInfo->Run);
  output().FillVar(subrun, GetSpill().EventInfo->SubRun);
  output().FillVar(evt,    GetSpill().EventInfo->Event);
  
  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_BeamTrue(output(), GetSpill().Beam);

  // Fill true kaons candidate info
  //const kaonAnaTrueVertex& kvtx = static_cast<const kaonAnaTrueVertex&>(vtx);
  kaonTree::FillKaonVariables_TrueKaonCandidates(output(), &part);
}

//********************************************************************
void secondaryKaonAnalysis::FillCategories(){
//********************************************************************

  // For the beam track
  if(box().MainTrack)
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils

  // for the candidates
  if(box().Candidates.size()>0){
    for(int i = 0; i < (int)box().Candidates.size(); i++){
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]),              "candidate",   1);
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]->Daughters[0]),"candidatedau",1);
      kaonAnaUtils::FillCandidateDaughterMuonCategory(box().Candidates[i], static_cast<AnaParticlePD*>(box().Candidates[i]->Daughters[0]));
    }
  }
  
  //get the kaon selection and get the branch with a larger AccumLevel
  SelectionBase* ksel = sel().GetSelection(_selection_name.c_str());
  int branchmax = 0;
  int max       = 0;  
  for(UInt_t ibranch = 0; ibranch < ksel->GetNBranches(); ibranch++){
    if(ksel->GetAccumCutLevel(ibranch) > max){
      max = ksel->GetAccumCutLevel(ibranch);
      branchmax = ibranch;
    }
  }

  // ---------- best kaon candidate variables --------------
  if(box().Candidates.size()>0){
    anaUtils::FillCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[branchmax]), "bestcandidate");
    anaUtils::FillCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[branchmax]->Daughters[0]), "bestcandidatedau");
  }


  // For the beam 
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart) anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
}

//********************************************************************
bool secondaryKaonAnalysis::FinalizeConfiguration(){
//********************************************************************
  
  // This is called before FillMicroTrees()

  const Int_t CandidatesAccumLevel = 2; //TODO: think about this
  
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
