#include "dEdxCalibration.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"

#include "dEdxTrackSelection.hxx"

#include "SCEVariation.hxx"

//********************************************************************
dEdxCalibration::dEdxCalibration(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

//********************************************************************
bool dEdxCalibration::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("dEdxCalibration.MinAccumLevelToSave"));

  return true;
}

//********************************************************************
void dEdxCalibration::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",           new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered",   new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("MichelRemovingTree", new MichelRemovingTreeConverter("michelremoving2/Event"));
}

//********************************************************************
void dEdxCalibration::DefineSelections(){
//********************************************************************

  sel().AddSelection("non stopping crossing tracks", "non stopping crossing tracks", new dEdxTrackSelection(false));  
}

//********************************************************************
void dEdxCalibration::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
}

//********************************************************************
void dEdxCalibration::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  //evar().AddEventVariation(0, "SCE variation", new SCEVariation());
}

//********************************************************************
void dEdxCalibration::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

  // // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  // if(_enableAllSystConfig){
  //   conf().EnableEventVariation(0,all_syst);
  // }
}

//********************************************************************
void dEdxCalibration::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
    
  AddVarMF(output(), track_hit_x   , "trk hit x"   , ntracks, -30, 3000);
  AddVarMF(output(), track_hit_y   , "trk hit y"   , ntracks, -30, 3000);
  AddVarMF(output(), track_hit_z   , "trk hit z"   , ntracks, -30, 3000);
  AddVarMF(output(), track_hit_dqdx, "trk hit dqdx", ntracks, -30, 3000);

  // -------- Add toy variables ---------------------------------
  // AddToyVarMF(output(), toy_track_hit_x   , "trk hit x"   , 5, 1000);
  // AddToyVarMF(output(), toy_track_hit_y   , "trk hit y"   , 5, 1000);
  // AddToyVarMF(output(), toy_track_hit_z   , "trk hit z"   , 5, 1000);
  // AddToyVarMF(output(), toy_track_hit_dqdx, "trk hit dqdx", 5, 1000);
}

//********************************************************************
void dEdxCalibration::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void dEdxCalibration::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 
  
  // fill trk info
  if(box().Tracks.size()>0){
    for(int i = 0; i < std::min((int)box().Tracks.size(),30); i++){
      AnaParticlePD* part = box().Tracks[i];
      for(int j = 0; j < (int)std::min((int)part->Hits[2].size(),3000); j++){
  	output().FillMatrixVar(track_hit_x,         (Float_t)part->Hits[2][j].Position.X(),   -1, j);
  	output().FillMatrixVar(track_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),   -1, j);
  	output().FillMatrixVar(track_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),   -1, j);
  	output().FillMatrixVar(track_hit_dqdx,      (Float_t)part->Hits[2][j].dQdx,           -1, j);
      }
      output().IncrementCounter(ntracks);
    } 
  }
}

//********************************************************************
void dEdxCalibration::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

  // int ntracks = box().Tracks.size();
  // if(ntracks > 0){
  //   for(int itrk = 0; itrk < std::min(ntracks,5); itrk++){
  //     AnaParticlePD* part = box().Tracks[itrk];
  //     int nhits = part->Hits[2].size();
  //     for(int ihit = 0; ihit < std::min(nhits,1000); ihit++){
  // 	output().FillToyMatrixVar(toy_track_hit_x   ,(Float_t) part->Hits[2][ihit].Position.X(),itrk,ihit);
  // 	output().FillToyMatrixVar(toy_track_hit_y   ,(Float_t) part->Hits[2][ihit].Position.Y(),itrk,ihit);
  // 	output().FillToyMatrixVar(toy_track_hit_z   ,(Float_t) part->Hits[2][ihit].Position.Z(),itrk,ihit);
  // 	output().FillToyMatrixVar(toy_track_hit_dqdx,(Float_t) part->Hits[2][ihit].dQdx_elife,itrk,ihit);
  //     }
  //   }
  // }
}

//********************************************************************
void dEdxCalibration::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void dEdxCalibration::FillCategories(){
//********************************************************************

}
