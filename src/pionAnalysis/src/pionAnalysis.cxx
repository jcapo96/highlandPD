#include "pionAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"

#include "dEdxVariation.hxx"
#include "dEdxCorrection.hxx"
#include "BeamMomCorrection.hxx"
#include "BeamMomSmearingCorrection.hxx"
#include "CryoWallBeamMomCorrection.hxx"
#include "dEdxDataCorrection.hxx"
#include "BrokenTrackDataCorrection.hxx"

#include "pionSelection.hxx"
#include "pionAnalysisUtils.hxx"
#include "pionTreeConverter.hxx"
#include "pdDataClasses.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "LArSoftTreeConverter.hxx"

// PionAnalysis from Jake and Francesca
// https://github.com/calcuttj/PionStudies.git

/*
  A highland Analysis inherits in ultimate instance from AnalysisAlgorithm. 
  In this particular case an intermediate class baseAnalysis is used, 
  which does all the work that is common for DUNE analysis 
  Then pionAnalysis inherits from baseAnalysis, which inherits from AnalysisAlgorithm.
  
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

const UInt_t NMAXSAVEDPARTICLES=700;
const UInt_t NMAXSAVEDDAUGHTERS=100;

PDCounters _globalCounters; 

//********************************************************************
pionAnalysis::pionAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version (not yet properly done)
  ND::versioning().AddPackage("pionAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("PIONANALYSISROOT")));
}

//********************************************************************
bool pionAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */
  
  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pionAnalysis.MinAccumLevelToSave"));

  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle

  // Add our own categories
  pionAnaUtils::AddCustomCategories();

  return true;
}

//********************************************************************
void pionAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("pionana",        new pionTreeConverter());
  input().AddConverter("minitree",       new HighlandMiniTreeConverter());
  input().AddConverter("LArSoftTree",    new LArSoftTreeConverter());
}

//********************************************************************
void pionAnalysis::DefineSelections(){
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

  sel().AddSelection("pionSelection", "pion selection", new pionSelection(false));     // true/false for forcing break  
}

//********************************************************************
void pionAnalysis::DefineCorrections(){
//********************************************************************

  /* Corrections modify some aspect of the input data (real data or MC). 
     The entire analysis will proceed on the modified data
     In general Corrections are used to reduce data/MC differences
  */
  
  // Some corrections are defined in highland (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineCorrections();

  // Define additional corrections (pionAnalysys/src/corrections)
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledEdx"))             corr().AddCorrection("dEdx",              new dEdxCorrection());
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledEdxData"))         corr().AddCorrection("dEdx data",         new dEdxDataCorrection());   
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledBeamMom"))         corr().AddCorrection("BeamMom",           new BeamMomCorrection());
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledBeamMomSmearing")) corr().AddCorrection("BeamMomSmearing",   new BeamMomSmearingCorrection());
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledCryoWallBeamMom")) corr().AddCorrection("CryoWallBeamMom",   new CryoWallBeamMomCorrection());
  if (ND::params().GetParameterI("pionAnalysis.Corrections.EnabledBrokenTrackData")) corr().AddCorrection("broken track data", new BrokenTrackDataCorrection());   
}

//********************************************************************
void pionAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  // Define additional systematics (pionAnalysys/src/systematics)
  evar().AddEventVariation(SystId::kdEdx,           "dEdx",           new dEdxVariation());
}

//********************************************************************
void pionAnalysis::DefineConfigurations(){
//********************************************************************

  /*  A "configuration" is defined by the systematics that are enabled, the number of toys being run and the random seed 
      used to generate the toys. Each configuration has a micro-tree associated in the output file (with the same name)
  */
  
  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if (_enableAllSystConfig){
    if (ND::params().GetParameterI("pionAnalysis.Variations.EnabledEdx"))        conf().EnableEventVariation(SystId::kdEdx,         all_syst);
  }
}

//********************************************************************
void pionAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  /*  We call Micro-trees to the standard analysis trees appearing in the output file. 
      There is always a Micro-Tree call "ana" which should contain the basic info to understand our selection. 
      The user can add extra Micro-Trees by adding configurations to the analysis (see DefineConfigurations method above).

      Here we give an example of different variables that can be added. Have a look at highland/src/highland2/highlandTools/src/OutputManager.hxx
      to see all available methods.
  */
  
  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for PD analysis
  standardPDTree::AddStandardVariables_CountersTrue(output());  
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
  standardPDTree::AddStandardVariables_CandidateTrue(output());
  standardPDTree::AddStandardVariables_CandidateReco(output());
  standardPDTree::AddStandardVariables_CandidateHitsReco(output());
  standardPDTree::AddStandardVariables_CandidateDaughtersTrue(output(),NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_CandidateDaughtersReco(output(),NMAXSAVEDDAUGHTERS);
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),NMAXSAVEDPARTICLES);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),NMAXSAVEDPARTICLES);

  // Additional variables
  AddVarF(  output(), seltrk_truebeta,    "candidate true beta");
  AddVarF(  output(), seltrk_truemom_tpc, "candidate true momentum at TPC entrance");
  AddVarF(  output(), seltrk_kinetic,     "candidate kinetic energy");
  AddVarI(  output(), seltrk_broken,      "candidate was broken");
  
  AddVar3VF(output(),seltrk_pos_z0,        "candidate start position extrapolated at z=0");
  AddVarF(output(),  seltrk_length_z0,     "candidate length assuming the particle starts at z=0");
  AddVarF(output(),  seltrk_mom_muon_z0,   "candidate reconstructed momentum at z=0 (muon)");
  AddVarF(output(),  seltrk_mom_prot_z0,   "candidate reconstructed momentum at z=0 (proton)");

  //  AddVarFixVF( output(), seltrk_dedx_exp,"candidate expected dEdx for different hypothesis",4);
  
  //new variables added for reproducing Jake/Fraceska plots
  AddVar3VF(       output(), seltrk_CNNscore,        "candidate reconstructed CNN score");
  AddToyVarF(      output(), seltrk_chi2_prot,       "candidate chi2 proton");
  AddToyVarF(      output(), seltrk_chi2_ndf,        "candidate chi2 ndf");

  seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMaxSize3MF(output(), seltrk_dau_CNNscore,    "candidate daughters reconstructed CNN score",seltrk_ndau,NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVF( output(), seltrk_dau_chi2_prot,   "candidate daughters chi2 proton",            seltrk_ndau,NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVF( output(), seltrk_dau_chi2_ndf,    "candidate daughters chi2 ndf",               seltrk_ndau,NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVF( output(), seltrk_dau_vtxdistance, "candidate daguhters distance to vtx",        seltrk_ndau,NMAXSAVEDDAUGHTERS);
  AddVarMaxSizeVI( output(), seltrk_dau_type,        "candidate daguhters type(track=0, shower=1)",seltrk_ndau,NMAXSAVEDDAUGHTERS);
}

//********************************************************************
void pionAnalysis::DefineTruthTree(){
//********************************************************************

  /*  The "truth" tree also appears in the output file. It contains all events in which we are interested in regardless on whether 
      the selection was passed or not. This is the tree that should be used to compute signal efficiencies
  */
  
  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for PD analysis
  standardPDTree::AddStandardVariables_CountersTrue(output());
  standardPDTree::AddStandardVariables_BeamTrue(    output());

  // Additional variables
  ntracks = standardPDTree::ntracks;
  AddVarI(  output(), ntracks,          "number of reconstructed tracks in the TPC");
}

//********************************************************************
void pionAnalysis::InitializeBunch(){
//********************************************************************

  /* In ProtoDUNE a beam Bunch corresponds to a single event. At the beginning of the event the counters are computed
   */
  
  pdAnaUtils::FillBeamDaughterCounters(GetEvent(), _globalCounters);
}

//********************************************************************
void pionAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-independent variables (all except the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called once all toys have been run, what means that the value of all variables for the last toy will be saved. This is not a problem 
      for variables that are not expected to change from a toy to another.
  */
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_CountersTrue(     output(), _globalCounters);  
  standardPDTree::FillStandardVariables_BeamReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_CandidateReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_CandidateTrue(    output(), box().MainTrack);    
  standardPDTree::FillStandardVariables_CandidateHitsReco(output(), box().MainTrack);
  
  // Get all reconstructed parts in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;
  
  // ---------- Save information about all (max NMAXSAVEDPARTICLES) recon parts in the event --------------
  // These are standard variables for the PD analysis
  for (Int_t i=0;i<std::min((Int_t)NMAXSAVEDPARTICLES,nParts); ++i){    
    AnaParticle* part = static_cast<AnaParticle*>(parts[i]);
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    output().IncrementCounter(ntracks);

    //    std::cout << "anselmo 2: " << static_cast<AnaParticle*>(part)->UniqueID << " " << static_cast<AnaParticle*>(part)->Length << " " << part->TrueObject << std::endl;
    if (part->TrueObject)
      static_cast<AnaTrueParticle*>(part->TrueObject)->Print();
  }

  // ---------- Additional candidate variables --------------
  if (box().MainTrack){        
    output().FillVar(seltrk_broken,                (Int_t)(box().MainTrack->NDOF==8888));
    output().FillVar(seltrk_kinetic,               pdAnaUtils::ComputeKineticEnergy(*box().MainTrack));              

    // Extrapolation to z=0
    Float_t posz[3]={0};
    output().FillVectorVarFromArray(seltrk_pos_z0, pdAnaUtils::ExtrapolateToZ(box().MainTrack,0,posz), 3); //default extrapolates to 0    

    Float_t diflength = sqrt(pow(posz[0]-box().MainTrack->PositionStart[0],2) + pow(posz[1]-box().MainTrack->PositionStart[1],2) + pow(posz[2]-box().MainTrack->PositionStart[2],2));
    output().FillVar(seltrk_length_z0,             box().MainTrack->Length + diflength);
    
    output().FillVar(seltrk_mom_muon_z0,           pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length + diflength, 13));
    output().FillVar(seltrk_mom_prot_z0,           pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length + diflength, 2212));
    
    //Jake/Franceska variable    
    output().FillVectorVarFromArray(seltrk_CNNscore,       box().MainTrack->CNNscore,3);

    // ---------- Save information about all (max NMAXSAVEDDAUGHTERS) daughters in the candidate --------------
    Int_t ndau = (Int_t)box().MainTrack->Daughters.size();
    for (Int_t i=0;i<std::min((Int_t)NMAXSAVEDDAUGHTERS,ndau); ++i){      
      AnaParticle* dau = static_cast<AnaParticle*>(box().MainTrack->Daughters[i]);

      // These are standard variables for the PD analysis
      standardPDTree::FillStandardVariables_CandidateDaughterReco(output(), dau);
      standardPDTree::FillStandardVariables_CandidateDaughterTrue(output(), dau);

      //Jake/Franceska variables
      output().FillMatrixVarFromArray(seltrk_dau_CNNscore,  dau->CNNscore,         3); 
      output().FillVectorVar(seltrk_dau_chi2_prot,          dau->Chi2Proton);
      output().FillVectorVar(seltrk_dau_chi2_ndf,           dau->Chi2ndf);
      output().FillVectorVar(seltrk_dau_vtxdistance,        box().DaughterDistanceToVertex[i]);
      if (dau->UniqueID!=-1 && dau->Chi2ndf>0)
        output().FillVectorVar(seltrk_dau_type,              0);
      else
        output().FillVectorVar(seltrk_dau_type,              1);
      
      output().IncrementCounter(seltrk_ndau);
    }    

    // ---------- Additional truth information -----------
    AnaTrueParticle* truePart = static_cast<AnaTrueParticle*>(box().MainTrack->TrueObject);
    if (truePart){
      output().FillVar(seltrk_truemom_tpc,              truePart->MomentumInTPC);
      if (anaUtils::GetParticleMass(ParticleId::GetParticle(truePart->PDG))>0)
        output().FillVar(seltrk_truebeta,          anaUtils::ComputeBetaGamma(truePart->Momentum, ParticleId::GetParticle(truePart->PDG)));
    }
  }
}

//********************************************************************
void pionAnalysis::FillToyVarsInMicroTrees(bool addBase){
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
  if (box().MainTrack){
    output().FillToyVar(seltrk_chi2_prot, box().MainTrack->Chi2Proton);
    output().FillToyVar(seltrk_chi2_ndf,  box().MainTrack->Chi2ndf);
  }  
}

//********************************************************************
bool pionAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  /* To avoid unecessary events in the "truth" tree in this method we define the condition to include or not a given 
     true vertex in the tree. 
  */

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways
  return true;
}

//********************************************************************
void pionAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_CountersTrue(output(), _globalCounters);
  standardPDTree::FillStandardVariables_BeamTrue(    output(), GetSpill().Beam);

  // Additional variables. Number of reconstructed tracks
  output().FillVar(ntracks,  (Int_t)static_cast<AnaBunchB*>(GetSpill().Bunches[0])->Particles.size());  
}

//********************************************************************
void pionAnalysis::FillCategories(){
//********************************************************************

  /* This method fills the micro-tree variables related with trackevent/object categories for color drawing. 
     Those variables are added automatically (no need to explicitely add them in DefineMicroTrees) to the 
     micro-trees, but must be filled by the analyser, provided the event and the relevant track 

     If this method is not implemented, the one from the base class (baseAnalysis::FillCategories()) will be called.      
  */

  // For the candidate
  if (box().MainTrack){
    anaUtils::FillCategories(          &GetEvent(), box().MainTrack,"");
    pionAnaUtils::FillCustomCategories(&GetEvent(), box().MainTrack, _globalCounters);
  }

  // For the beam track
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart) anaUtils::FillCategories(&GetEvent(), beamPart, "beam");
}

