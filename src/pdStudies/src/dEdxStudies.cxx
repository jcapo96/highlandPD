#include "dEdxStudies.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"

#include "pdDataClasses.hxx"

#include "hitPionTreeConverter.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "LArSoftTreeConverter.hxx"
#include "pandoraPreselection.hxx"


const UInt_t NMAXSAVEDPARTICLES=700;


//********************************************************************
dEdxStudies::dEdxStudies(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version (not yet properly done)
  ND::versioning().AddPackage("pdStudies", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("PDSTUDIESROOT")));
}

//********************************************************************
bool dEdxStudies::Initialize(){
//********************************************************************

  /* In this method we Initialize everything that requires access to parameters in the parameters file. 
     This is because in order to the overwride parameters file 
     option (-p param.dat) to work, parameters cannot be accessed in the constructors. 
  */
  
  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("dEdxStudies.MinAccumLevelToSave"));

  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle


  // Initialize the MVA algorithm
  //  fMVA.Initialize();
  

  return true;
}

//********************************************************************
void dEdxStudies::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("pionana",        new hitPionTreeConverter());
  input().AddConverter("minitree",       new HighlandMiniTreeConverter("kaonana/MiniTree"));
  input().AddConverter("LArSoftTree",    new LArSoftTreeConverter());
}

//********************************************************************
void dEdxStudies::DefineSelections(){
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

  sel().AddSelection("pandoraPreselection", "pandora selection", new pandoraPreselection(false));     // true/false for forcing break  
}

//********************************************************************
void dEdxStudies::DefineCorrections(){
//********************************************************************

  /* Corrections modify some aspect of the input data (real data or MC). 
     The entire analysis will proceed on the modified data
     In general Corrections are used to reduce data/MC differences
  */
  
  // Some corrections are defined in highland (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineCorrections();

}

//********************************************************************
void dEdxStudies::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

}

//********************************************************************
void dEdxStudies::DefineConfigurations(){
//********************************************************************

  /*  A "configuration" is defined by the systematics that are enabled, the number of toys being run and the random seed 
      used to generate the toys. Each configuration has a micro-tree associated in the output file (with the same name)
  */
  
  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

}

//********************************************************************
void dEdxStudies::DefineMicroTrees(bool addBase){
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
  standardPDTree::AddStandardVariables_BeamReco(output());
  standardPDTree::AddStandardVariables_BeamTrue(output());
  standardPDTree::AddStandardVariables_CandidateTrue(output());
  standardPDTree::AddStandardVariables_CandidateReco(output());
  standardPDTree::AddStandardVariables_CandidateHitsReco(output());
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),NMAXSAVEDPARTICLES);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),NMAXSAVEDPARTICLES);

  AddVarF(    output(), seltrk_dedx_fromdqdx, "candidate average dEdx from dqdx");

  AddVarMaxSizeVF( output(), seltrk_hit_dedx_0,        "hit dedx 0",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_dedx_1,        "hit dedx 1",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_dedx_2,        "hit dedx 2",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_dedx_3,        "hit dedx 3",seltrk_nhits,500);

  AddVarMaxSizeVF( output(), seltrk_hit_dedx_0_cal,        "hit dedx 0",seltrk_nhits,500);
  
  AddVarMaxSizeVF( output(), seltrk_hit_dqdx_0,        "hit dqdx 0",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_dqdx_3,        "hit dqdx 3",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_dqdx_3_cal,        "hit dqdx 3 calibrated",seltrk_nhits,500);

  AddVarMaxSizeVF( output(), seltrk_hit_time,          "hit time",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_integral,      "hit integral",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_amplitude,     "hit amplitude",seltrk_nhits,500);
  AddVarMaxSizeVF( output(), seltrk_hit_pitch3D,       "hit pitch3D",seltrk_nhits,500);

  
}

//********************************************************************
void dEdxStudies::DefineTruthTree(){
//********************************************************************

  /*  The "truth" tree also appears in the output file. It contains all events in which we are interested in regardless on whether 
      the selection was passed or not. This is the tree that should be used to compute signal efficiencies
  */
  
  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  // Additional variables
  ntracks = standardPDTree::ntracks;
  AddVarI(  output(), ntracks,          "number of reconstructed tracks in the TPC");
}

//********************************************************************
void dEdxStudies::InitializeBunch(){
//********************************************************************

  /* In ProtoDUNE a beam Bunch corresponds to a single event. At the beginning of the event the counters are computed
   */
  
}

//********************************************************************
void dEdxStudies::FillMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-independent variables (all except the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called once all toys have been run, what means that the value of all variables for the last toy will be saved. This is not a problem 
      for variables that are not expected to change from a toy to another.
  */


  //  std::vector<hl::MVAPIDResult> result;
  //  fMVA.RunPID(GetEvent(),result);
  //  _cnn.produce(*(static_cast<AnaEventPD*>(&GetEvent())));
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // Fill standard variables for the PD analysis
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
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
    standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
    output().IncrementCounter(ntracks);
  }

  // ---------- Additional candidate variables --------------
  if (box().MainTrack){        
    FilldEdxInfo(box().MainTrack);
  }
}

//********************************************************************
void dEdxStudies::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  /*  In this method we fill all toy-dependent variables (the ones added with AddToy...) defined in the method DefineMicroTrees. 
      This method is called at the end of each toy.

      There could be many variables that are toy-dependent. We don't need to save all of them as toy variables, but only the ones we are interested in plotting 
      for different toys. 

      TOY VARIABLES ARE VERY SPACE CONSUMING SO WE SHOULD MINIMISE ITS NUMBER !!!!
  */

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
bool dEdxStudies::CheckFillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  /* To avoid unecessary events in the "truth" tree in this method we define the condition to include or not a given 
     true vertex in the tree. 
  */

  (void) vtx; // to avoid warning for unused vtx variable
  
  // fill it allways
  return true;
}

//********************************************************************
void dEdxStudies::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);

  // Fill standard variables for the PD analysis
  standardPDTree::FillStandardVariables_BeamTrue(    output(), GetSpill().Beam);

  // Additional variables. Number of reconstructed tracks
  output().FillVar(ntracks,  (Int_t)static_cast<AnaBunchB*>(GetSpill().Bunches[0])->Particles.size());  
}

//********************************************************************
void dEdxStudies::FillCategories(){
//********************************************************************

  /* This method fills the micro-tree variables related with trackevent/object categories for color drawing. 
     Those variables are added automatically (no need to explicitely add them in DefineMicroTrees) to the 
     micro-trees, but must be filled by the analyser, provided the event and the relevant track 

     If this method is not implemented, the one from the base class (baseAnalysis::FillCategories()) will be called.      
  */

  // For the candidate
  if (box().MainTrack){
    anaUtils::FillCategories(          &GetEvent(), box().MainTrack,"");
  }

  // For the beam track
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart) anaUtils::FillCategories(&GetEvent(), beamPart, "beam");
}

//********************************************************************
void dEdxStudies::FilldEdxInfo(AnaParticlePD* part){
//********************************************************************

  bool debug = true;
  
  /*
    1. Original dEdX from tree
    2. Compute dEdX from dQdx taken from tree. Use pdAnaUtils::ComputeDeDxFromDqDx
    3. Compute dEdX from dQdx taken from tree. Use fCalo.dEdx_from_dQdx_e. Has e- lifetime correction
    4. Compute dEdx from dQdx, computed from hit::Integral. Use fCalo.dEdx_from_dQdx_e. Has e- lifetime correction
    
    
    
    2 & 3 are equivalent except for a small difference in Wion. Both use 4.57e-3 calib factor
    4 uses 
  */
  
  
  double fEventT0=500;
  TVector3 dir = anaUtils::ArrayToTVector3(part->DirectionStart);
  Int_t plane = 2;
  Float_t pitch3D = pdAnaUtils::Compute3DWirePitch(plane, dir);
  
  for (UInt_t j=0;j<part->Hits[plane].size();j++){

    AnaHitPD& hit = part->Hits[plane][j];

    if (debug)
      std::cout << "(1)  -->  dEdxStudies: dQdx, dEdx = " << hit.dQdx << ", " << hit.dEdx << std::endl;


    if (debug)
      std::cout << "(2) dEdxStudies: pdAnaUtils::ComputeDeDxFromDqDx. dQdx = " << hit.dQdx << std::endl;

    float dedx_2 = pdAnaUtils::ComputeDeDxFromDqDx(hit.dQdx, plane, hit.Position.X(), hit.Position.Y(), hit.Position.Z());

    if (debug)
      std::cout << "(2)  -->  dEdxStudies: pdAnaUtils::ComputeDeDxFromDqDx. dEdx = " << dedx_2 << std::endl;

    if (debug)
      std::cout << "(3) dEdxStudies: fCalo.dEdx_from_dQdx. dQdx = " << hit.dQdx << std::endl;
    
    float dedx_3 = fCalo.dEdx_from_dQdx(hit.dQdx,
                                        hit.PeakTime,
                                        fEventT0);
    if (debug)
      std::cout << "(3)  --> dEdxStudies: fCalo.dEdx_from_dQdx. dEdx = " << dedx_3 << std::endl;

    if (debug)
      std::cout << "(4) dEdxStudies: fCalo.dEdx_AREA " << std::endl;
    
    float dedx_4 = fCalo.dEdx_AREA(hit, pitch3D, fEventT0);

    if (debug)
      std::cout << "ALL RESULTS: " << hit.dEdx << " " << dedx_2 <<  " " << dedx_3 <<  " " << dedx_4 << std::endl;

      output().FillVectorVar(seltrk_hit_dedx_0, hit.dEdx);
      output().FillVectorVar(seltrk_hit_dedx_1, dedx_2);
      output().FillVectorVar(seltrk_hit_dedx_2, dedx_3);
      output().FillVectorVar(seltrk_hit_dedx_3, dedx_4);

      output().FillVectorVar(seltrk_hit_dedx_0_cal, hit.dEdx_corr);


      output().FillVectorVar(seltrk_hit_dqdx_0, hit.dQdx);



      output().FillVectorVar(seltrk_hit_dqdx_3, hit.Integral/pitch3D);


      Float_t cal_dQdx_e = pdAnaUtils::ComputeCalibrateddQdX(hit.Integral/pitch3D, hit.Position);
      output().FillVectorVar(seltrk_hit_dqdx_3_cal, cal_dQdx_e);      

      output().FillVectorVar(seltrk_hit_time,     (Float_t)hit.PeakTime);
      output().FillVectorVar(seltrk_hit_integral, (Float_t)hit.Integral);
      output().FillVectorVar(seltrk_hit_amplitude,(Float_t)hit.PeakAmplitude);
      output().FillVectorVar(seltrk_hit_pitch3D,  (Float_t)pitch3D);
      
      output().IncrementCounter(seltrk_nhits);
      

  }
}
