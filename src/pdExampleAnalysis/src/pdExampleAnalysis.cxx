#include "pdExampleAnalysis.hxx"
#include "Parameters.hxx"
#include "pdExampleSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"

#include "pdAnalysisUtils.hxx"

#include "PDSPAnalyzerTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"

#include "baseToyMaker.hxx"

/*
  A highland Analysis inherits in ultimate instance from AnalysisAlgorithm. 
  In this particular case an intermediate class baseAnalysis is used, 
  which does all the work that is common for DUNE analysis 
  Then pdExampleAnalysis inherits from baseAnalysis, which inherits from AnalysisAlgorithm.
  
  The class that does the loops (over events, over toys, etc) is AnalysisLoop (under highlandTools). 
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

  - FillToyVarsInMicroTrees: Fill the micro-trees toy-dependent variables (run for each toy experiment)

  These are the methods that are run for Event (Bunch)

  - FillMicroTrees:          Fill the micro-trees toy-independent variables. (MANDATORY) 
  - FillCategories:          Fill the micro-tree variables used for color code drawing (OPTIONAL)

  These are the mandatory methods that are run for each Spill

  - FillTruthTree:           Fill the truth tree variables
  - CheckFillTruthTree:      Check whether to include or not a given true vertex in the truth tree

  These are the methods that are run at initialization level

  - FillConfigTree:           Fill the user defined variables in the single-entry config tree
  
*/


//********************************************************************
pdExampleAnalysis::pdExampleAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version
  //  ND::versioning().AddPackage("pdExampleAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
}

//********************************************************************
bool pdExampleAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */

  // Initialize the baseAnalysis
  if(!baseAnalysis::Initialize()) return false;

  // Save information about all tracks?
  _saveAllTracks = ND::params().GetParameterI("pdExampleAnalysis.SaveAllTracks");

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pdExampleAnalysis.MinAccumLevelToSave"));

  // Define categories for color drawing. Have a look at highland/src/highland2/highlandUtils/src/CategoriesUtils.hxx
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardObjectCategories("dau",standardPDTree::seltrk_ndau,"seltrk_ndau",1);

  return true;
}

//********************************************************************
void pdExampleAnalysis::DefineInputConverters(){
//********************************************************************
  
  /* In this method we add the to the InputManager (accessed by input() ) the InputConverters created 
     in separet files (see, for example, pdIO/PDSPAnalyzerTreeConverter.hxx)
     which define the allowed input file formats.
  */
  
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("pduneana",         new PDSPAnalyzerTreeConverter()); 

}

//********************************************************************
void pdExampleAnalysis::DefineSelections(){
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

  sel().AddSelection("pdExampleSelection", "protoDuneExample selection", new pdExampleSelection(false)); // true/false for forcing break
  
}

//********************************************************************
void pdExampleAnalysis::DefineCorrections(){
//********************************************************************

  /* Corrections modify some aspect of the input data (real data or MC). 
     The entire analysis will proceed on the modified data
  */

  // Some corrections are defined in baseAnalysis
  baseAnalysis::DefineCorrections();

}

//********************************************************************
void pdExampleAnalysis::DefineSystematics(){
//********************************************************************
  
  /* There are two types of systematic propagation methods: variations and weights, which are added to two different managers 
     The EventVariationManager (access with evar()) and the EventWeightManager (access with eweight()).     
  */
  
  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

}

//********************************************************************
void pdExampleAnalysis::DefineConfigurations(){
//********************************************************************

  /*  A "configuration" is defined by the systematics that are enabled, the number of toys being run and the random seed 
      used to generate the toys. Each configuration has a micro-tree associated in the output file (with the same name)
  */

  // Some configurations are defined in baseAnalysis
  baseAnalysis::DefineConfigurations();

  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if(_enableAllSystConfig){
    if(ND::params().GetParameterI("pdExampleAnalysis.Systematics.EnableHit"))
      conf().EnableEventVariation(SystId::kdEdx,         all_syst);
  }
}

//********************************************************************
void pdExampleAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  /*  We call Micro-trees to the standard analysis trees appearing in the output file. 
      There is always a Micro-Tree call "ana" which should contain the basic info to understand our selection. 
      The user can add extra Micro-Trees by adding configurations to the analysis (see DefineConfigurations method above).

      Here we give an example of different variables that can be added. Have a look at highland/src/highland2/highlandTools/src/OutputManager.hxx
      to see all available methods.
  */

  
  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis (these methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  // beam instrumentation info
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  
  // candidate track (beam track) info
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());
  standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
  				      
  // candidate's daughters info	      
  standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),pdExampleAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),pdExampleAnalysisConstants::NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersHitsReco(output(),pdExampleAnalysisConstants::NMAXSAVEDDAUGHTERS);
  
  // all tracks info
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),pdExampleAnalysisConstants::NMAXSAVEDPARTICLES);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),pdExampleAnalysisConstants::NMAXSAVEDPARTICLES);

  // Add additional variables manually, (should have been defined in pdExampleAnalysis.hxx)
  // Have a look at highland/src/highland2/highlandTools/src/OutputManager.hxx to see all available methods.
  AddVarF(output(), seltrk_csdarange_prot, "CSDA range proton hypothesis from beam momentum");

  // set counters for vector/matrix arrays
  seltrk_ndau = standardPDTree::seltrk_ndau;
  ntracks     = standardPDTree::ntracks;
  
}

//********************************************************************
void pdExampleAnalysis::DefineTruthTree(){
//********************************************************************

  /*  The "truth" tree also appears in the output file. It contains all interactions in which we are interested in regardless on whether 
      the selection was passed or not. This is the tree that should be used to compute signal efficiencies
  */

  // Variables from baseAnalysis (run, event, ...)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for PD analysis
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());

}

//********************************************************************
void pdExampleAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-independent variables (all except the ones added with AddToy...) defined in the method DefineMicroTrees.
      This method is called once all toys have been run, what means that the value of all variables for the last toy will be saved. 
      This is not a problem for variables that are not expected to change from a toy to another.
  */

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);    
  standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);

  // Get all reconstructed parts in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;

  // ---------- Save information about all (max NMAXSAVEDPARTICLES) recon parts in the event --------------
  Int_t maxtracks = std::min((Int_t)pdExampleAnalysisConstants::NMAXSAVEDPARTICLES,nParts);
  for(Int_t i = 0; i < maxtracks; i++){    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    output().IncrementCounter(ntracks);
  }
    
  // ---------- Additional candidate variables --------------
  if(box().MainTrack){        
    
    AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    if(beamPart){
      output().FillVar(seltrk_csdarange_prot, pdAnaUtils::ComputeCSDARange(1000*(beamPart->Momentum), 2212));
    }

    // ---------- Save information about all (max NMAXSAVEDDAUGHTERS) daughters in the candidate --------------
    Int_t ndau   = (Int_t)box().MainTrack->Daughters.size();
    Int_t maxdau = std::min((Int_t)pdExampleAnalysisConstants::NMAXSAVEDDAUGHTERS,ndau);
    for(Int_t i = 0; i < maxdau; i++){      
      AnaParticlePD* dau = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]);
      standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), dau);
      standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), dau);
      output().IncrementCounter(seltrk_ndau);
    }    
  }

}

//********************************************************************
void pdExampleAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-dependent variables (the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called at the end of each toy.

      There could be many variables that are toy-dependent. We don't need to save all of them as toy variables, 
      but only the ones we are interested in plotting for different toys. 

      TOY VARIABLES ARE VERY SPACE CONSUMING SO WE SHOULD MINIMISE ITS NUMBER !!!!
  */


  // Fill the common variables
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
  
}
//********************************************************************
bool pdExampleAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  /* To avoid unecessary events in the "truth" tree in this method we define the condition to include or not a given 
     true vertex in the tree. 
  */

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways for the moment
  return true;
}

//********************************************************************
void pdExampleAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables
  baseAnalysis::FillTruthTreeBase(vtx);

  // Fill standard variables for the PD analysis
  //  standardPDTree::FillStandardVariables_CountersTrue(output(), _globalCounters);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(    output(), GetSpill().Beam);

}

//********************************************************************
void pdExampleAnalysis::FillCategories(){
//********************************************************************

  /* This method fills the micro-tree variables related with track categories for color drawing. 
     Those variables are added automatically (no need to explicitely add them in DefineMicroTrees) to the 
     micro-trees, but must be filled by the analyser, provided the event and the relevant track 

     If this method is not implemented, the one from the base class (baseAnalysis::FillCategories()) will be called.      
  */

  
  // For the candidate
  if(box().MainTrack){
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,"");
    
    // Categories for daughters
    if(box().MainTrack->Daughters.size()>0){
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().MainTrack->Daughters[0]),"dau",1);
    }
  }
  
  // For the beam track
  AnaParticleB* beam = static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle;
  if(beam)anaUtils::FillCategories(&GetEvent(), beam,"beam");
}
