#include "kaonAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"

#include "dEdxCorrection.hxx"
#include "BeamMomCorrection.hxx"
#include "BeamMomSmearingCorrection.hxx"
#include "CryoWallBeamMomCorrection.hxx"
#include "dEdxDataCorrection.hxx"
#include "BrokenTrackDataCorrection.hxx"

#include "LengthVariation.hxx"
#include "BeamCompositionWeight.hxx"
#include "LifetimeVariation.hxx"
#include "dQdxCalibVariation.hxx"
#include "dEdxCalibVariation.hxx"
#include "RecombinationVariation.hxx"
#include "TrackEffWeight.hxx"

#include "kaonSelection.hxx"
#include "kaonAnalysisUtils.hxx"
#include "pdDataClasses.hxx"
#include "kaonClasses.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

/*
  A highland Analysis inherits in ultimate instance from AnalysisAlgorithm. 
  In this particular case an intermediate class baseAnalysis is used, 
  which does all the work that is common for ProtoDUNE analysis 
  Then kaonAnalysis inherits from baseAnalysis, which inherits from AnalysisAlgorithm.
  
  The class that does the loops (over events, over toys, etc) is AnalysisLoop (under highland/src/highland/highlandTools/src). 
  There you can see the analysis flow, which is as follows

  LOOP over Spills{
    InitializeSpill
    LOOP over Bunches in the Spill{
      InitializeBunch
      LOOP over Configurations{
        InitializeConfiguration
        LOOP over Toy experiments for each configuration{
          InitializeToy
          LOOP over Enabled Selections{
            InitializeSelection
            Apply The Selection
            FinalizeSelection
          }
          FinalizeToy
          FillMicroTrees (Fill toy-dependent variables in micro-trees)
        }
        FillToyVarsInMicroTrees (Fill toy-independent variables)
        FillCategories (Fill categories for color code drawing)
        FinalizeConfiguration
      }
      FinalizeBunch
    }
    FillTruthTree (Fill the truth tree)
    FinalizeSpill
  }

  The Initialize.. and Finalize... methods can be implemented by the user to do more complicated analyses, but are not mandatory


  These is the list of mandatory methods to configure the analysis (call at initialization level):
  
  - Initialize:           This is all initializations that cannot be done in the constructor should be put 
                          (the ones that require access to a parameter in the parameters file)
  - DefineSelections:     Add to the SelectionManager the selections  we have defined in other files
  - DefineSystematics:    Add to the SystematicManager the systematics to be run (defined in other files)
  - DefineCorrections:    Add to the CorrectionManager the corrections to be run (defined in other files)
  - DefineConfigurations: Add to the ConfigurationManager the configurations to be considered
  - DefineMicroTrees:     Define the micro-tree variables
  - DefineTruthTree:      Define the variables of the truth tree

  These are the mandatory methods that are run for each Toy Experiment 

  - FillToyVarsInMicroTrees: Fill the micro-trees toy-dependent   variables (run for each toy experiment)

  These are the methods that are run for Event (Bunch)

  - FillMicroTrees:          Fill the micro-trees toy-independent variables. (MANDATORY) 
  - FillCategories:          Fill the micro-tree variables used for color code drawing (OPTIONAL)

  These are the mandatory methods that are run for each Spill

  - FillTruthTree:           Fill the truth tree variables
  - CheckFillTruthTree:      Check whether to include or not  a given true vertex in the truth tree

  These are the methods that are run at initialization level

  - FillConfigTree:           Fill the user defined variables in the single-entry config tree
  
*/


//********************************************************************
kaonAnalysis::kaonAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version (not yet properly done)
  //  ND::versioning().AddPackage("kaonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("KAONANALYSISROOT")));
}

//********************************************************************
bool kaonAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */
  
  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("kaonAnalysis.MinAccumLevelToSave"));

  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle
  anaUtils::AddStandardCategories("bestcandidate");  // This is for the best candidate of each event
  anaUtils::AddStandardCategories("bestcandidatedau");  // This is for the best candidate of each event
  anaUtils::AddStandardObjectCategories("daughter",standardPDTree::seltrk_ndau,"seltrk_ndau",1);  // This is for the daughters of the selected track
  anaUtils::AddStandardObjectCategories("all",standardPDTree::ntracks,"ntracks",1);  // This is for all the tracks in the event

  // Add standard categories for the candidates
  anaUtils::AddStandardObjectCategories("candidate"   ,kaonTree::candidates,"candidates",1);
  anaUtils::AddStandardObjectCategories("candidatedau",kaonTree::candidates,"candidates",1);

  // Add our own categories (in kaonAnalysisUtils.cxx)
  kaonAnaUtils::AddCustomCategories();

  return true;
}

//********************************************************************
bool kaonAnalysis::InitializeSpill(){
//********************************************************************
  
  /*this method is not mandatory, but filling the truth tree requires at least to have a true vertex
    we are using this method to create and add a dummy vertex to the spill
  */
  if(!baseAnalysis::InitializeSpill())return false;

  // Create a new truevertex for each true kaon
  for(int i = 1; i < (int)GetSpill().TrueParticles.size(); i++){ //skip first particle because it is the primary
    if(GetSpill().TrueParticles[i]->PDG==321){
      kaonAnaTrueVertex* vtx = new kaonAnaTrueVertex();
      //fill the vertex with the kaon
      vtx->TrueParticlesVect.clear();
      vtx->TrueParticlesVect.push_back(GetSpill().TrueParticles[i]);
      //set the bunch value
      vtx->Bunch = 1;
      //get the end of the kaon
      std::pair<Int_t, Int_t> result(0,-999);
      if     (GetSpill().TrueParticles[i]->ProcessEnd==2){
	result = kaonAnaUtils::GetKaonDecayMode(GetSpill().TrueParticles[i], GetSpill().TrueParticles);
	vtx->DecayMode = result.first;
      }
      else if(GetSpill().TrueParticles[i]->ProcessEnd==3){
	result = kaonAnaUtils::MuonFromKaonChain(GetSpill().TrueParticles[i], GetSpill().TrueParticles);
	vtx->ChainMuon = result.first;
      }
      //if it finish in a muon, store it also in the truevertex
      if(vtx->DecayMode == 1 || vtx->ChainMuon == 1)
	vtx->TrueParticlesVect.push_back(pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles,result.second));
      GetSpill().TrueVertices.push_back(static_cast<AnaTrueVertex*>(vtx));
    }
  }

  return true;
}

//********************************************************************
void kaonAnalysis::DefineInputConverters(){
//********************************************************************

  /*  InputConverters are the objects that read the input file (any format) and fill the internal highland event model.
      The user can write its own converter for a given input file type. 
   */
  
  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void kaonAnalysis::DefineSelections(){
//********************************************************************

  /* In this method the user will add to the SelectionManager (accessed by  sel() ) the selections to be run, 
     defined in other files. In general an analysis has a single selection, which could have multiple branches. 
     Each of the branches is in fact a different selection (different cuts and actions), but defining them as branches 
     we ussualy gain in speed and simplicity since the steps that are common are applied only once. 
     Sometimes the user does not want to expend time merging two selections in a single one (as branches), but preffers to run 
     both independently. This is in general done for quick checks (are the selections mutually exclusive ?, etc). 
     This is possible in highland2. An example on how to do that is shown below. 
  */

  // Add selections to the SelectionManager provided:
  // - Name of the selection 
  // - Title, the one dumped by the DrawingTools::DumpSelections() method. It is an explaination of the selection
  // - Pointer to the selection. The argument in the constructor (false) indicates the 
  //   step sequence is not broken when a cut is not passed. In this way we can save intermediate information for events 
  //   not passing the entire selection

  sel().AddSelection("kaonSelection", "kaon selection", new kaonSelection(false));     // true/false for forcing break  
}

//********************************************************************
void kaonAnalysis::DefineCorrections(){
//********************************************************************

  /* Corrections modify some aspect of the input data (real data or MC). 
     The entire analysis will proceed on the modified data
     In general Corrections are used to reduce data/MC differences
  */
  
  // Some corrections are defined in highland (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineCorrections();

  // Define additional corrections (kaonAnalysys/src/corrections)
}

//********************************************************************
void kaonAnalysis::DefineSystematics(){
//********************************************************************

  /* There are two types of systematic propagation methods: variations and weights, which are added to two different managers 
     The EventVariationManager (access with evar()) and the EventWeightManager (access with eweight()).     
   */
  
  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  //---- Define additional systematics (kaonAnalysys/src/systematics) -----
  eweight().AddEventWeight(kBeam, "beamComp", new BeamCompositionWeight());
}

//********************************************************************
void kaonAnalysis::DefineConfigurations(){
//********************************************************************

  /*  A "configuration" is defined by the systematics that are enabled, the number of toys being run and the random seed 
      used to generate the toys. Each configuration has a micro-tree associated in the output file (with the same name)
  */
  
  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if(_enableAllSystConfig){
    if(ND::params().GetParameterI("kaonAnalysis.Systematics.EnabledBeamComposition"))conf().EnableEventWeight(kBeam,all_syst);
  }
}

//********************************************************************
void kaonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  /*  We call Micro-trees to the standard analysis trees appearing in the output file. 
      There is always a Micro-Tree called "ana" which should contain the basic info to understand our selection. 
      The user can add extra Micro-Trees by adding configurations to the analysis (see DefineConfigurations method above).

      Here we give an example of different variables that can be added.
      Have a look at highland/src/highland2/highlandTools/src/OutputManager.hxx to see all available methods.
  */
  
  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  //standardPDTree::AddStandardVariables_CountersTrue(output());  
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
  standardPDTree::AddStandardVariables_CandidateTrue(output());
  standardPDTree::AddStandardVariables_CandidateReco(output());
  //standardPDTree::AddStandardVariables_CandidateHitsReco(output());
  standardPDTree::AddStandardVariables_CandidateDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateDaughtersReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateDaughtersHitsReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS, kaonAnalysisConstants::NMAXSAVEDHITSDAUGHTERS);
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),kaonAnalysisConstants::NMAXSAVEDPARTICLES);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),kaonAnalysisConstants::NMAXSAVEDPARTICLES);

  // -------- Add additional variables to the analysis tree ----------------------
  kaonTree::AddKaonVariables_CandidateDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  kaonTree::AddKaonVariables_CandidateGDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);
  kaonTree::AddKaonVariables_CandidateGDaughtersReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);
  kaonTree::AddKaonVariables_CandidateGGDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS, kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS);
  kaonTree::AddKaonVariables_CandidateGGDaughtersReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS, kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS);

  // -------- Add candidates variables ----------------------
  kaonTree::AddKaonVariables_KaonCandidatesReco    (output(),kaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesHitsReco(output(),kaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  kaonTree::AddKaonVariables_KaonCandidatesTrue    (output(),kaonAnalysisConstants::NMAXSAVEDCANDIDATES);

  // -------- Add best candidate variables ----------------------
  kaonTree::AddKaonVariables_KaonBestCandidateReco    (output());
  kaonTree::AddKaonVariables_KaonBestCandidateHitsReco(output());
  kaonTree::AddKaonVariables_KaonBestCandidateTrue    (output());
}

//********************************************************************
void kaonAnalysis::DefineTruthTree(){
//********************************************************************

  /*  The "truth" tree also appears in the output file. It contains all events in which we are interested in regardless on whether 
      the selection was passed or not. This is the tree that should be used to compute signal efficiencies
  */
  
  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for ProtoDUNE analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_CountersTrue(output());
  standardPDTree::AddStandardVariables_BeamTrue(    output());
  
  // Add specific variables for the kaon candidates
  kaonTree::AddKaonVariables_TrueKaonCandidates(output(), kaonAnalysisConstants::NMAXSAVEDTRUECANDIDATES);
}

//********************************************************************
void kaonAnalysis::InitializeBunch(){
//********************************************************************

  /* In ProtoDUNE a beam Bunch corresponds to a single event. At the beginning of the event the counters are computed
   */
}

//********************************************************************
void kaonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-independent variables (all except the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called once all toys have been run, what means that the value of all variables for the last toy will be saved. This is not a problem 
      for variables that are not expected to change from a toy to another.
  */

  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_BeamReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_CandidateReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_CandidateTrue(    output(), box().MainTrack);    
  //standardPDTree::FillStandardVariables_CandidateHitsReco(output(), box().MainTrack);
  
  // Get all reconstructed parts in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;
  
  // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDPARTICLES) recon parts in the event --------------
  // These are standard variables for the PD analysis
  for (Int_t i=0;i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDPARTICLES,nParts); ++i){    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //skip seltrk 
    if(part->isPandora)continue;
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    output().IncrementCounter(standardPDTree::ntracks);
  }

  // ---------- Additional candidate variables --------------
  if (box().MainTrack){        
 
    // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDDAUGHTERS) daughters in the candidate --------------
    Int_t ndau = (Int_t)box().MainTrack->Daughters.size();
    for (Int_t idau = 0; idau < std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,ndau); idau++){      
      AnaParticlePD* dau = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[idau]);
      // These are standard variables for the PD analysis
      standardPDTree::FillStandardVariables_CandidateDaughterReco(output(), dau);
      standardPDTree::FillStandardVariables_CandidateDaughterTrue(output(), dau);
      standardPDTree::FillStandardVariables_CandidateDaughterHitsReco(output(), dau, kaonAnalysisConstants::NMAXSAVEDHITSDAUGHTERS);
      // Additional variables for kaon Analysis
      kaonTree::FillKaonVariables_CandidateDaughterTrue(output(), box().MainTrack, dau);
       
      // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS) gdaughters in the candidate --------------
      Int_t ngdau = (Int_t)dau->DaughtersIDs.size();
      for (Int_t jgdau = 0; jgdau < std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS,ngdau); jgdau++){
	AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[jgdau]);
	kaonTree::FillKaonVariables_CandidateGDaughterReco(output(), gdau, jgdau);
	kaonTree::FillKaonVariables_CandidateGDaughterTrue(output(), gdau, jgdau);
	
	// ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS) gdaughters in the candidate --------------
	Int_t nggdau = (Int_t)gdau->DaughtersIDs.size();
	for (Int_t jggdau = 0; jggdau < (Int_t)std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS,nggdau); jggdau++){
	  AnaParticlePD* ggdau = static_cast<AnaParticlePD*>(gdau->Daughters[jggdau]);
	  kaonTree::FillKaonVariables_CandidateGGDaughterReco(output(), ggdau, jgdau, jggdau);
	  kaonTree::FillKaonVariables_CandidateGGDaughterTrue(output(), ggdau, jgdau, jggdau);
	}
      }
      output().IncrementCounter(standardPDTree::seltrk_ndau);
    } 
  }

  // ---------- kaon candidates variables --------------
  if(box().Candidates.size()>0){
    for(int i = 0; i < std::min((int)box().Candidates.size(),(int)kaonAnalysisConstants::NMAXSAVEDCANDIDATES); i++){
      kaonTree::FillKaonVariables_KaonCandidatesReco    (output(), box().Candidates[i]);
      kaonTree::FillKaonVariables_KaonCandidatesHitsReco(output(), box().Candidates[i]);
      kaonTree::FillKaonVariables_KaonCandidatesTrue    (output(), box().Candidates[i]);
      output().IncrementCounter(kaonTree::candidates);
    } 
  }

  //get the kaon selection and get the branch with a larger AccumLevel
  SelectionBase* ksel = sel().GetSelection("kaonSelection");
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
    kaonTree::FillKaonVariables_KaonBestCandidateReco    (output(), box().Candidates[branchmax]);
    kaonTree::FillKaonVariables_KaonBestCandidateHitsReco(output(), box().Candidates[branchmax]);
    kaonTree::FillKaonVariables_KaonBestCandidateTrue    (output(), box().Candidates[branchmax]);
  }

  //associate candidates with true kaons
  for(int i = 0; i < (int)box().Candidates.size(); i++){
    if(!box().Candidates[i]->TrueObject)continue;
    for(int j = 0; j < (int)GetSpill().TrueVertices.size(); j++){
      if(box().Candidates[i]->TrueObject->ID == GetSpill().TrueVertices[j]->TrueParticlesVect[0]->ID){
	static_cast<kaonAnaTrueVertex*>(GetSpill().TrueVertices[j])->Branch = i;
      }
    }
  }

  //since accum_level is for candidates, true kaon candidates have an assigned accum_level
  //despite having a candidate match. This gives problems with the efficiencies, because a 
  //truekaon can have an accum_level assigned from a non matching candidate.
  //For that reason we reset the accum_levels of true candidates without matching
  //and the accum_level of truekaon with matching for all non-matching branches
  for(int i = 0; i < (int)GetSpill().TrueVertices.size(); i++){
    kaonAnaTrueVertex* kvtx = static_cast<kaonAnaTrueVertex*>(GetSpill().TrueVertices[i]);
    if(kvtx->Branch==-999){
      kvtx->AccumLevel[0][0] = std::min(*min_element(kvtx->AccumLevel[0].begin(), kvtx->AccumLevel[0].end()),2);
      std::fill(kvtx->AccumLevel[0].begin()+1, kvtx->AccumLevel[0].end(), -2);
    }
    else {
      kvtx->AccumLevel[0][0] = kvtx->AccumLevel[0][kvtx->Branch];
      std::fill(kvtx->AccumLevel[0].begin()+1, kvtx->AccumLevel[0].end(), -2);
    }
  }
}

//********************************************************************
void kaonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-dependent variables (the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called at the end of each toy.

      There could be many variables that are toy-dependent. We don't need to save all of them as toy variables, but only the ones we are interested in plotting 
      for different toys. 

      TOY VARIABLES ARE VERY SPACE CONSUMING SO WE SHOULD MINIMISE ITS NUMBER !!!!
  */

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

  // Toy variables, since there is a systematic variation associated to them

}

//********************************************************************
bool kaonAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  /* To avoid unecessary events in the "truth" tree in this method we define the condition to include or not a given 
     true vertex in the tree. 
  */

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways for the moment
  return true;
}

//********************************************************************
void kaonAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);

  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_BeamTrue(    output(), GetSpill().Beam);

  // Fill true kaons candidate info
  const kaonAnaTrueVertex& kvtx = static_cast<const kaonAnaTrueVertex&>(vtx);
  kaonTree::FillKaonVariables_TrueKaonCandidates(output(), kvtx);
}

//********************************************************************
void kaonAnalysis::FillCategories(){
//********************************************************************

  /* This method fills the micro-tree variables related with event/object categories for color drawing. 
     Those variables are added automatically (no need to explicitely add them in DefineMicroTrees) to the 
     micro-trees, but must be filled by the analyser, provided the event and the relevant object 

     If this method is not implemented, the one from the base class (baseAnalysis::FillCategories()) will be called.      
  */

  // For the candidate
  if (box().MainTrack){
    anaUtils::FillCategories(          &GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils
    kaonAnaUtils::FillCustomCategories(&GetEvent(), box().MainTrack); 
    for(int i = 0; i < std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,(Int_t)box().MainTrack->Daughters.size()); i++){
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().MainTrack->Daughters[i]),"daughter",1);
    }
  }

  // for the candidates
  if(box().Candidates.size()>0){
    for(int i = 0; i < (int)box().Candidates.size(); i++){
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]),              "candidate",   1);
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]->Daughters[0]),"candidatedau",1);
      kaonAnaUtils::FillGDaughterKaonCategory("candidatedaumuon",&GetEvent(), box().Candidates[i], static_cast<AnaParticlePD*>(box().Candidates[i]->Daughters[0]),-1,-1);
    }
  }

  //get the kaon selection and get the branch with a larger AccumLevel
  SelectionBase* ksel = sel().GetSelection("kaonSelection");
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


  // For the beam track
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart) anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils

  // For all particles in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;
  for (Int_t i=0;i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDPARTICLES,nParts); ++i){    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //skip seltrk 
    if(part->isPandora)continue;
    anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(part), "all", 1);
  }
}

//********************************************************************
void kaonAnalysis::FinalizeSpill(){
//********************************************************************
  
  /*this method is not mandatory, but filling the truth tree requires at least to have a true vertex in the spill.
    Since the truevertices are not clonned, they are only deleted in the RawSpill. The spill we created was added
    to the corrected spill and it has to be deleted from there. Were are using this method to delete it
  */
  baseAnalysis::FinalizeSpill();
  for(int i = 0; i < (int)GetSpill().TrueVertices.size(); i++)
    delete GetSpill().TrueVertices[i];
  
  GetSpill().TrueVertices.clear();
}
