#include "stoppingProtonAnalysis.hxx"
#include "Parameters.hxx"
//#include "stoppingKaonSelection.hxx"
#include "stoppingProtonSelection.hxx"
#include "stoppingMuonSelection.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"

#include "dEdxVariation.hxx"
#include "dEdxCorrection.hxx"
#include "BeamMomCorrection.hxx"
#include "BeamMomSmearingCorrection.hxx"
#include "CryoWallBeamMomCorrection.hxx"
#include "dEdxDataCorrection.hxx"
#include "BrokenTrackDataCorrection.hxx"

#include "pdAnalysisUtils.hxx"

#include "pionTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "LArSoftTreeConverter.hxx"

#include "baseToyMaker.hxx"

/*
  A highland Analysis inherits in ultimate instance from AnalysisAlgorithm. 
  In this particular case an intermediate class baseAnalysis is used, 
  which does all the work that is common for DUNE analysis 
  Then stoppingProtonAnalysis inherits from baseAnalysis, which inherits from AnalysisAlgorithm.
  
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


//********************************************************************
stoppingProtonAnalysis::stoppingProtonAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version
  ND::versioning().AddPackage("stoppingProtonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("STOPPINGPROTONANALYSISROOT")));
}

//********************************************************************
bool stoppingProtonAnalysis::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */

  
  // Initialize the baseAnalysis
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("stoppingProtonAnalysis.MinAccumLevelToSave"));

  // Define categories for color drawing
  anaUtils::AddStandardCategories();
  anaUtils::AddStandardCategories("beam");
  anaUtils::AddStandardCategories("dau");

  return true;
}

//********************************************************************
void stoppingProtonAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("pionana",        new pionTreeConverter());
  input().AddConverter("minitree",       new HighlandMiniTreeConverter());
  input().AddConverter("LArSoftTree",    new LArSoftTreeConverter());
}

//********************************************************************
void stoppingProtonAnalysis::DefineSelections(){
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


  sel().AddSelection("stoppingProtonSelection",  "protoDuneExample selection",     new stoppingProtonSelection(false));     // true/false for forcing break
  sel().AddSelection("stoppingMuonSelection",    "protoDuneExample selection",     new stoppingMuonSelection(false));
  
  //sel().AddSelection("stoppingKaonSelection",    "protoDuneExample selection",     new stoppingKaonSelection(false));     // true/false for forcing break
  
  
}

//********************************************************************
void stoppingProtonAnalysis::DefineCorrections(){
//********************************************************************

  /* Corrections modify some aspect of the input data (real data or MC). 
     The entire analysis will proceed on the modified data
  */
  
  // Some corrections are defined in baseAnalysis
  baseAnalysis::DefineCorrections();

  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledEdx"))             corr().AddCorrection("dEdx",              new dEdxCorrection());
  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledEdxData"))         corr().AddCorrection("dEdx data",         new dEdxDataCorrection());   
  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledBeamMom"))         corr().AddCorrection("BeamMom",           new BeamMomCorrection());
  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledBeamMomSmearing")) corr().AddCorrection("BeamMomSmearing",   new BeamMomSmearingCorrection());
  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledCryoWallBeamMom")) corr().AddCorrection("CryoWallBeamMom",   new CryoWallBeamMomCorrection());
  if (ND::params().GetParameterI("stoppingProtonAnalysis.Corrections.EnabledBrokenTrackData")) corr().AddCorrection("broken track data", new BrokenTrackDataCorrection());   

}

//********************************************************************
void stoppingProtonAnalysis::DefineSystematics(){
//********************************************************************

  // Some corrections are defined in baseAnalysis
  baseAnalysis::DefineSystematics();

  evar().AddEventVariation(SystId::kdEdx,           "dEdx",           new dEdxVariation());
}

//********************************************************************
void stoppingProtonAnalysis::DefineConfigurations(){
//********************************************************************

  /*  A "configuration" is defined by the systematics that are enabled, the number of toys being run and the random seed 
      used to generate the toys. Each configuration has a micro-tree associated in the output file (with the same name)
  */

  
  // Some configurations are defined in baseAnalysis
  baseAnalysis::DefineConfigurations();

  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if (_enableAllSystConfig){
    if (ND::params().GetParameterI("stoppingProtonAnalysis.Variations.EnabledEdx"))        conf().EnableEventVariation(SystId::kdEdx,         all_syst);
  }
}

//********************************************************************
void stoppingProtonAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  /*  We call Micro-trees to the standard analysis trees appearing in the output file. 
      There is always a Micro-Tree call "default" which should contain the basic info to understand our selection. 
      The user can add extra Micro-Trees by adding configurations to the analysis (see DefineConfigurations method above).

      Here we give an example of different variables that can be added. Have a look at highlandTools/vXrY/src/OutputManager.hxx
      to see all available methods.
  */

  
  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)
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

  
  // is this a signal event ? 
  AddVarI(  output(), true_signal,       "is true signal ?");

  AddVarF(  output(), seltrk_truebeta,    "candidate true beta");

  seltrk_ndau = standardPDTree::seltrk_ndau;
  AddVarMF(output(), seltrk_dau_dedx_binned,  "daughters dEdx binned",        seltrk_ndau,-100,4);

  AddVarF(  output(), seltrk_truemom_tpc,  "candidate true momentum in the TPC");
  
  // --- candidate recon variables ---
  AddVarI(  output(), seltrk_broken,     "candidate was broken");
  AddVarF(  output(), seltrk_kinetic,    "kinetic energy");

  AddVar3VF(output(),seltrk_pos_z0,        "candidate start position extrapolated at z=0");
  AddVarF(output(), seltrk_length_z0,     "candidate length assuming the particle starts at z=0");
  AddVarF(output(), seltrk_length_sce,    "candidate length corrected by sce");
  AddVarF(output(), seltrk_mom_muon_z0,   "candidate reconstructed momentum at z=0 (muon)");
  AddVarF(output(), seltrk_mom_prot_z0,   "candidate reconstructed momentum at z=0 (proton)");
  
  AddVarF(output(), seltrk_csdarange_muon,          "candidate reconstructed csdarange (muon) from beam mom");
  AddVarF(output(), seltrk_csdarange_prot,          "candidate reconstructed csdarange (proton) from beam mom");
  AddVarF(output(), seltrk_csdarange_muon_raw,      "candidate reconstructed csdarange (muon) from beam mom with no BI correction");
  AddVarF(output(), seltrk_csdarange_prot_raw,      "candidate reconstructed csdarange (proton) from beam mom with no BI correction");
  AddVarF(output(), seltrk_csdarange_tpc_muon,      "candidate reconstructed csdarange (muon) from cryo mom");
  AddVarF(output(), seltrk_csdarange_tpc_prot,      "candidate reconstructed csdarange (proton) from cryo mom");
  
  //  AddVarF(output(), seltrk_theta,    "candidate reconstructed polar angle");

  AddVarFixMF(output(), seltrk_dedx_binned,     "candidate average dEdx for several resrange bins",3,4);
  
  AddVarF(   output(), seltrk_pida_raw,  "candidate PIDA (raw)");
  AddVarF(   output(), seltrk_pida_corr, "candidate PIDA (corrected)");
  AddToyVarF(output(), seltrk_pida,      "candidate PIDA (variated)");

  AddVarFixVF(  output(), seltrk_pida2,     "pida2", 3);
  //  AddVarFixVI(  output(), seltrk_pdg,       "recon pdg", 3);
}

//********************************************************************
void stoppingProtonAnalysis::DefineTruthTree(){
//********************************************************************

  /*  The "truth" tree also appears in the output file. It contains all interactions in which we are interested in regardless on whether 
      the selection was passed or not. This is the tree that should be used to compute signal efficiencies
  */

  
  // Variables from baseAnalysis (run, event, ...)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for PD analysis
  standardPDTree::AddStandardVariables_CountersTrue(output());
  standardPDTree::AddStandardVariables_BeamTrue(    output());

  // Additional variables
  ntracks = standardPDTree::ntracks;
  AddVarI(  output(), ntracks,          "number of reconstructed tracks in the TPC");
}

//********************************************************************
void stoppingProtonAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-independent variables (all except the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called once all toys has been run, what means that the value of all variables for the last toy will be saved. This is not a problem 
      for variables that are not expected to change from a toy to another.
  */

  // Variables from baseAnalysis (run, event, ...)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  // The true signal for protoDUNE requires the initial beam particle to reach the active volume
  if (GetEvent().nTrueParticles>0){
    if (GetEvent().TrueParticles[0]->PositionEnd[2]>0)     output().FillVar(true_signal, 1);
    else output().FillVar(true_signal, 0);
  }

  // Fill standard variables for the PD analysis
  //  standardPDTree::FillStandardVariables_CountersTrue(     output(), _globalCounters);  
  standardPDTree::FillStandardVariables_BeamReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_CandidateReco(    output(), box().MainTrack);
  if(box().MainTrack)output().FillVar(seltrk_length_sce,            box().MainTrack->corrected_Length);
  standardPDTree::FillStandardVariables_CandidateTrue(    output(), box().MainTrack);    
  standardPDTree::FillStandardVariables_CandidateHitsReco(output(), box().MainTrack);

  // Get all reconstructed parts in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nParts         = GetEvent().nParticles;
  
  // ---------- Save information about all (max NMAXSAVEDPARTICLES) recon parts in the event --------------
  // These are standard variables for the PD analysis
  for (Int_t i=0;i<std::min((Int_t)NMAXSAVEDPARTICLES,nParts); ++i){    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    //    output().FillVar(trk_pandora,(Int_t)part->isPandora);
    output().IncrementCounter(ntracks);
  }

  // ---------- Additional candidate variables --------------
  if (box().MainTrack){        
    output().FillVar(seltrk_broken,                (Int_t)(box().MainTrack->NDOF==8888));
    output().FillVar(seltrk_kinetic,               pdAnaUtils::ComputeKineticEnergy(*box().MainTrack));              

    // This is the PIDA computed from the hits, and not the one precomputed in the Analysis Tree
    if (box().MainTrack->Original){
      if (box().MainTrack->Original->Original){
        output().FillVar(seltrk_pida_corr,  pdAnaUtils::ComputePIDA(*static_cast<const AnaParticlePD*>(box().MainTrack->Original->Original)));      
      }
      if (box().MainTrack->Original->Original->Original){
        output().FillVar(seltrk_pida_raw,  pdAnaUtils::ComputePIDA(*static_cast<const AnaParticlePD*>(box().MainTrack->Original->Original->Original)));
      }
    }

    // Extrapolation to z=0
    Float_t posz[3]={0};
    output().FillVectorVarFromArray(seltrk_pos_z0, pdAnaUtils::ExtrapolateToZ(box().MainTrack,0,posz), 3); //default extrapolates to 0    

    Float_t diflength = sqrt(pow(posz[0]-box().MainTrack->PositionStart[0],2) + pow(posz[1]-box().MainTrack->PositionStart[1],2) + pow(posz[2]-box().MainTrack->PositionStart[2],2));
    output().FillVar(seltrk_length_z0,             box().MainTrack->Length + diflength);
    
    //    output().FillVar(seltrk_mom_muon_z0,           pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length + diflength, 13));
    //    output().FillVar(seltrk_mom_prot_z0,           pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length + diflength, 2212));

    output().FillVectorVarFromArray(seltrk_pida2,  box().MainTrack->PIDA,3);
    //    output().FillVectorVarFromArray(seltrk_pdg,    box().MainTrack->ReconPDG,3);

    AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
    if(beamPart){
      //Float_t mom = beam->BeamMomentumInTPC;
      //if (beamPart->TrueObject) mom =  static_cast<AnaTrueParticlePD*>(beamPart->TrueObject)->Momentum;
      output().FillVar(seltrk_csdarange_muon,         pdAnaUtils::ComputeCSDARange(1000*(beamPart->Momentum), 13));
      output().FillVar(seltrk_csdarange_prot,         pdAnaUtils::ComputeCSDARange(1000*(beamPart->Momentum), 2212));
      if (beamPart->TrueObject){
        output().FillVar(seltrk_csdarange_muon_raw,     pdAnaUtils::ComputeCSDARange(1000*(static_cast<AnaTrueParticlePD*>(beamPart->TrueObject)->Momentum), 13));
        output().FillVar(seltrk_csdarange_prot_raw,     pdAnaUtils::ComputeCSDARange(1000*(static_cast<AnaTrueParticlePD*>(beamPart->TrueObject)->Momentum), 2212));
      }
      else{
        output().FillVar(seltrk_csdarange_muon_raw,     pdAnaUtils::ComputeCSDARange(1000*(static_cast<const AnaParticlePD*>(beamPart->Original->Original)->Momentum), 13));
        output().FillVar(seltrk_csdarange_prot_raw,     pdAnaUtils::ComputeCSDARange(1000*(static_cast<const AnaParticlePD*>(beamPart->Original->Original)->Momentum), 2212));
      }
      output().FillVar(seltrk_csdarange_tpc_muon,     pdAnaUtils::ComputeCSDARange(1000*((Float_t)beam->BeamMomentumInTPC), 13));
      output().FillVar(seltrk_csdarange_tpc_prot,     pdAnaUtils::ComputeCSDARange(1000*((Float_t)beam->BeamMomentumInTPC), 2212));
    }

    //    output().FillVar(seltrk_mom_muon2,         pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length/10., 13));
    //    output().FillVar(seltrk_mom_prot2,         pdAnaUtils::ComputeRangeMomentum(box().MainTrack->Length/10., 2212));
    

    Float_t *binned_dedx[3];
    for(int i = 0; i < 3; i++)
      binned_dedx[i] = new Float_t[4];

    pdAnaUtils::ComputeBinnedDeDx(static_cast<const AnaParticlePD*>(box().MainTrack),20,4,binned_dedx);
    for (int i=0;i<3;i++){
      output().FillMatrixVarFromArray(seltrk_dedx_binned,  binned_dedx[i] , i, 4);
    }
    
    for(int i = 0; i < 3; i++)
      delete binned_dedx[i];

    // ---------- Save information about all (max NMAXSAVEDDAUGHTERS) daughters in the candidate --------------
    Int_t ndau = (Int_t)box().MainTrack->Daughters.size();
    for (Int_t i=0;i<std::min((Int_t)NMAXSAVEDDAUGHTERS,ndau); ++i){      
      AnaParticlePD* dau = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]);

      // These are standard variables for the PD analysis
      standardPDTree::FillStandardVariables_CandidateDaughterReco(output(), dau);
      standardPDTree::FillStandardVariables_CandidateDaughterTrue(output(), dau);
      
      /*
      if (dau->UniqueID!=-1 && dau->Chi2ndf>0)
        output().FillVectorVar(seltrk_dau_type,              0);
      else
        output().FillVectorVar(seltrk_dau_type,              1);
      */

      Float_t *dau_binned_dedx[3];
      for(int i = 0; i < 3; i++)
        dau_binned_dedx[i] = new Float_t[4];

      pdAnaUtils::ComputeBinnedDeDx(dau,20,4,dau_binned_dedx);
      output().FillMatrixVarFromArray(seltrk_dau_dedx_binned,  dau_binned_dedx[2], 4);

      for(int i = 0; i < 3; i++)
        delete dau_binned_dedx[i];

      output().IncrementCounter(seltrk_ndau);
    }    

    // ---------- Additional truth information -----------
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(box().MainTrack->TrueObject);
    if (truePart){
      output().FillVar(seltrk_truemom_tpc,              truePart->MomentumInTPC);
      if (anaUtils::GetParticleMass(ParticleId::GetParticle(truePart->PDG))>0)
        output().FillVar(seltrk_truebeta,          anaUtils::ComputeBetaGamma(truePart->Momentum, ParticleId::GetParticle(truePart->PDG)));
    }
  }

}

//********************************************************************
void stoppingProtonAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-dependent variables (the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called at the end of each toy.

      There could be many variables that are toy-dependent. We don't need to save all of them as toy variables, but only the ones we are interested in plotting 
      for different toys. 

      TOY VARIABLES ARE VERY SPACE CONSUMING SO WE SHOULD MINIMISE ITS NUMBER !!!!
  */


  // Fill the common variables
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

  // PIDA has to be a toy variable, since there is a systematic variation associated to it
  if (box().MainTrack){
    output().FillToyVar(seltrk_pida,  pdAnaUtils::ComputePIDA(*box().MainTrack));
  }
  
}
//********************************************************************
bool stoppingProtonAnalysis::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  /* To avoid unecessary events in the "truth" tree in this method we define the condition to include or not a given 
     true vertex in the tree. 
  */

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways
  return true;
}

//********************************************************************
void stoppingProtonAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables
  baseAnalysis::FillTruthTreeBase(vtx);

  // The true signal for protoDUNE requires the initial beam particle to reach the active volume
  //  output().FillVar(true_signal, 0);

  // Fill standard variables for the PD analysis
  //  standardPDTree::FillStandardVariables_CountersTrue(output(), _globalCounters);
  standardPDTree::FillStandardVariables_BeamTrue(    output(), GetSpill().Beam);

  // Additional variables. Number of reconstructed tracks
  output().FillVar(ntracks,  (Int_t)static_cast<AnaBunchB*>(GetSpill().Bunches[0])->Particles.size());  
}

//********************************************************************
void stoppingProtonAnalysis::FillCategories(){
//********************************************************************

  /* This method fills the micro-tree variables related with track categories for color drawing. 
     Those variables are added automatically (no need to explicitely add them in DefineMicroTrees) to the 
     micro-trees, but must be filled by the analyser, provided the event and the relevant track 

     If this method is not implemented, the one from the base class (baseAnalysis::FillCategories()) will be called.      
  */

  
  // For the candidate
  if (box().MainTrack){
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,"");

    // Categories for first daugther of the candidate
    if (box().MainTrack->Daughters.size()>0) anaUtils::FillCategories(&GetEvent(), static_cast<AnaParticleB*>(box().MainTrack->Daughters[0]),"dau");
  }
  
  // For the beam track
  AnaParticleB* beam = static_cast<AnaBeamPD*>(GetSpill().Beam)->BeamParticle;
  if (beam) anaUtils::FillCategories(&GetEvent(), beam,"beam");
}
