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

#include "HighlandMiniTreeConverter.hxx"

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
  ND::versioning().AddPackage("kaonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("KAONANALYSISROOT")));
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

  // Add our own categories (in kaonAnalysisUtils.cxx)
  kaonAnaUtils::AddCustomCategories();

  return true;
}

//********************************************************************
void kaonAnalysis::DefineInputConverters(){
//********************************************************************

  /*  InputConverters are the objects that read the input file (any format) and fill the internal highland event model.
      The user can write its own converter for a given input file type. 
   */
  
  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",       new HighlandMiniTreeConverter("highlandana/MiniTree"));
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
}

//********************************************************************
void kaonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  /*  We call Micro-trees to the standard analysis trees appearing in the output file. 
      There is always a Micro-Tree called "ana" which should contain the basic info to understand our selection. 
      The user can add extra Micro-Trees by adding configurations to the analysis (see DefineConfigurations method above).

      Here we give an example of different variables that can be added. Have a look at highland/src/highland2/highlandTools/src/OutputManager.hxx
      to see all available methods.
  */
  
  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_CountersTrue(output());  
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
  standardPDTree::AddStandardVariables_CandidateTrue(output());
  standardPDTree::AddStandardVariables_CandidateReco(output());
  standardPDTree::AddStandardVariables_CandidateHitsReco(output());
  standardPDTree::AddStandardVariables_CandidateDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateDaughtersReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateGDaughtersTrue(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateGDaughtersReco(output(),kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),kaonAnalysisConstants::NMAXSAVEDPARTICLES);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),kaonAnalysisConstants::NMAXSAVEDPARTICLES);

  // -------- Add additional variables to the analysis tree ----------------------
  
  AddVar3VF(       output(), seltrk_CNNscore,        "candidate reconstructed CNN score");
  AddVarI(         output(), seltrk_truedaukaons,    "candidate daughters that are kaons");
  AddVarI(         output(), seltrk_truedaukaon_nmu, "true grandaughters antimuons associated to a secondary kaon");
  AddToyVarF(      output(), seltrk_chi2_prot,       "candidate chi2 proton");
  AddToyVarF(      output(), seltrk_chi2_ndf,        "candidate chi2 ndf");
    
  seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMaxSize3MF(output(), seltrk_dau_CNNscore,    "candidate daughters reconstructed CNN score",seltrk_ndau,kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVF( output(), seltrk_dau_chi2_prot,   "candidate daughters chi2 proton",            seltrk_ndau,kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVF( output(), seltrk_dau_chi2_ndf,    "candidate daughters chi2 ndf",               seltrk_ndau,kaonAnalysisConstants::NMAXSAVEDDAUGHTERS);
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

  // Additional variables
  ntracks = standardPDTree::ntracks;
  AddVarI(  output(), ntracks,          "number of reconstructed tracks in the TPC");
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
  standardPDTree::FillStandardVariables_CandidateHitsReco(output(), box().MainTrack);
  
  // Get all reconstructed parts in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;
  
  // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDPARTICLES) recon parts in the event --------------
  // These are standard variables for the PD analysis
  for (Int_t i=0;i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDPARTICLES,nParts); ++i){    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    output().IncrementCounter(ntracks);
  }

  // ---------- Additional candidate variables --------------
  if (box().MainTrack){        
    output().FillVectorVarFromArray(seltrk_CNNscore,       box().MainTrack->CNNscore,3);

    // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDDAUGHTERS) daughters in the candidate --------------
    Int_t ndau = (Int_t)box().MainTrack->Daughters.size();
    for (Int_t i=0;i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,ndau); ++i){      
      AnaParticlePD* dau = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]);

      // These are standard variables for the PD analysis
      standardPDTree::FillStandardVariables_CandidateDaughterReco(output(), dau);
      standardPDTree::FillStandardVariables_CandidateDaughterTrue(output(), dau);

      output().FillMatrixVarFromArray(seltrk_dau_CNNscore,  dau->CNNscore,         3); 
      output().FillVectorVar(seltrk_dau_chi2_prot,          dau->Chi2Proton);
      output().FillVectorVar(seltrk_dau_chi2_ndf,           dau->Chi2ndf);

      // ---------- Save information about all (max kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS) gdaughters in the candidate --------------

      // Associate to each daughter all gdaughters (pointers) using the DaughterIDs
      /*dau->Daughters.clear();
      for (UInt_t j = 0; j < dau->DaughtersIDs.size(); j++){
	AnaParticleB* gdau = anaUtils::GetParticleByID(GetEvent(),dau->DaughtersIDs[j]);
	if (gdau) {dau->Daughters.push_back(gdau); gdau->Print();}      
	}*/
      

      std::cout << "MIGUE " << dau->Daughters.size() << std::endl;
      Int_t ngdau = (Int_t)dau->Daughters.size();
      for(int  j = 0; j < std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS,ngdau); j++){
      AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[i]);
      //standardPDTree::FillStandardVariables_CandidateGDaughterReco(output(), gdau, j);
      //standardPDTree::FillStandardVariables_CandidateGDaughterTrue(output(), gdau, j);
      gdau->Print();
      }

      output().IncrementCounter(seltrk_ndau);
    } 

    // ---------- Additional truth information -----------
    int nkaons = 0;
    int nmuons = 0;
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
    if(truePart){
      for(int i = 0; i < (Int_t)truePart->Daughters.size(); i++){
	AnaTrueParticlePD* trueDau = pdAnaUtils::GetTrueParticle(&GetEvent(),truePart->Daughters.at(i));
	if(trueDau){
	  if(trueDau->PDG==321){
	    nkaons++;
	    for(int j = 0; j < (Int_t)trueDau->Daughters.size(); j++){
	      AnaTrueParticlePD* trueGDau = pdAnaUtils::GetTrueParticle(&GetEvent(),trueDau->Daughters.at(j));
	      if(trueGDau){
		if(trueGDau->PDG==-13)nmuons++;
	      }
	    }
	  }
	}
      }
    }
    output().FillVar(seltrk_truedaukaons,    nkaons);
    output().FillVar(seltrk_truedaukaon_nmu, nmuons);
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

  // Additional variables. Number of reconstructed tracks
  output().FillVar(ntracks,  (Int_t)static_cast<AnaBunchB*>(GetSpill().Bunches[0])->Particles.size());  
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
  }

  // For the beam track
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart) anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
}

