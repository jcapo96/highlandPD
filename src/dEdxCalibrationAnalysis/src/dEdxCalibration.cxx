//header
#include "dEdxCalibration.hxx"

//highland
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"
#include "dEdxUtils.hxx"

//converters
#include "HighlandMiniTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"

//selections
#include "dEdxTrackSelection.hxx"

//corrections
#include "HitPositionSCECorrection.hxx"
#include "HitPitchSCECorrection.hxx"
#include "SCECorrection.hxx"
#include "LifetimeCorrection.hxx"

//systematics
#include "SCEVariation.hxx"
#include "LifetimeVariation.hxx"

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

  _SaveAna = ND::params().GetParameterI("dEdxCalibration.SaveAna");
  _SaveToy = ND::params().GetParameterI("dEdxCalibration.SaveToy");

  _ApplySCEPositionCorrection = ND::params().GetParameterI("dEdxCalibration.ApplySCEPositionCorrection");
  _ApplySCEPitchCorrection = ND::params().GetParameterI("dEdxCalibration.ApplySCEPitchCorrection");
  _ApplyLifetimeCorrection = ND::params().GetParameterI("dEdxCalibration.ApplyLifetimeCorrection");

  _ApplySCESystematic = ND::params().GetParameterI("dEdxCalibration.ApplySCESystematic");
  _ApplyLifetimeSystematic = ND::params().GetParameterI("dEdxCalibration.ApplyLifetimeSystematic");

  TTree::SetMaxTreeSize(200000000000);

  return true;
}

//********************************************************************
void dEdxCalibration::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  _MichelRemovingTree = ND::params().GetParameterI("dEdxCalibration.MichelRemovingTree");
  std::stringstream treenumber;
  if(_MichelRemovingTree==0)treenumber.str("");
  else                      treenumber << _MichelRemovingTree;

  input().AddConverter("MichelRemovingTree", new MichelRemovingTreeConverter(("michelremoving"+treenumber.str()+"/Event").c_str()));
  input().AddConverter("minitreefiltered",   new HighlandMiniTreeConverter("MiniTree"));
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
  if(_ApplySCEPositionCorrection && _ApplySCEPitchCorrection)
    corr().AddCorrection(kSCE_corr, "SCE correction"  , new SCECorrection());
  else{
    if(_ApplySCEPositionCorrection)
      corr().AddCorrection(kSCEPosition_corr, "SCEpos"  , new HitPositionSCECorrection());
    if(_ApplySCEPitchCorrection)   
      corr().AddCorrection(kSCEPitch_corr, "SCEpitch", new HitPitchSCECorrection());
  }
  if(_ApplyLifetimeCorrection)   
    corr().AddCorrection(kLifetime_corr, "Lifetime", new LifetimeCorrection());
}

//********************************************************************
void dEdxCalibration::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  if(_ApplySCESystematic)     evar().AddEventVariation(kSCE_syst     , "SCE variation"     , new SCEVariation());
  if(_ApplyLifetimeSystematic)evar().AddEventVariation(kLifetime_syst, "lifetime variation", new LifetimeVariation());
}

//********************************************************************
void dEdxCalibration::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();

  //enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if(_enableAllSystConfig){
    if(_ApplySCESystematic)     conf().EnableEventVariation(kSCE_syst,all_syst);
    if(_ApplyLifetimeSystematic)conf().EnableEventVariation(kLifetime_syst,all_syst);
  }
}

//********************************************************************
void dEdxCalibration::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
  
  if(_SaveAna){
    AddVarMaxSizeVF(output(), hit_x   , "hit x"   , nhits, 3000);
    AddVarMaxSizeVF(output(), hit_y   , "hit y"   , nhits, 3000);
    AddVarMaxSizeVF(output(), hit_z   , "hit z"   , nhits, 3000);
    AddVarMaxSizeVF(output(), hit_dqdx, "hit dqdx", nhits, 3000);
  }
  
  // -------- Add toy variables ---------------------------------
  if(_SaveToy){
    AddToyVarVF(output(), toy_hit_x   , "toy hit x"   , 1000);
    AddToyVarVF(output(), toy_hit_y   , "toy hit y"   , 1000);
    AddToyVarVF(output(), toy_hit_z   , "toy hit z"   , 1000);
    AddToyVarVF(output(), toy_hit_dqdx, "toy hit dqdx", 1000);
  }
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
  int savedhits = 0;
  if(_SaveAna){
    if(box().Tracks.size()>0){
      for(int i = 0; i < (int)box().Tracks.size(); i++){
	AnaParticlePD* part = box().Tracks[i];
	for(int j = 0; j < (int)part->Hits[2].size(); j++){
	  if(savedhits>=3000)break;
	  output().FillVectorVar(hit_x,         (Float_t)part->Hits[2][j].Position.X());
	  output().FillVectorVar(hit_y,         (Float_t)part->Hits[2][j].Position.Y());
	  output().FillVectorVar(hit_z,         (Float_t)part->Hits[2][j].Position.Z());
	  output().FillVectorVar(hit_dqdx,      (Float_t)part->Hits[2][j].dQdx_SCE  );
	  output().IncrementCounter(nhits);
	  savedhits++;
	}
	if(savedhits>=3000)break;
      } 
    }
  }
}

//********************************************************************
void dEdxCalibration::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

  if(_SaveToy){
    int ntracks = box().Tracks.size();
    if(ntracks > 0){
      int hitcounter = 0;
      for(int itrk = 0; itrk < ntracks; itrk++){
	AnaParticlePD* part = box().Tracks[itrk];
	int nhits = part->Hits[2].size();
	for(int ihit = 0; ihit < nhits; ihit++){
	  if(hitcounter>=1000)break;
	  if(dEdxUtils::IsInterestingHit(part->Hits[2][ihit])){
	    output().FillToyVectorVar(toy_hit_x   ,(Float_t)part->Hits[2][ihit].Position.X(),hitcounter);
	    output().FillToyVectorVar(toy_hit_y   ,(Float_t)part->Hits[2][ihit].Position.Y(),hitcounter);
	    output().FillToyVectorVar(toy_hit_z   ,(Float_t)part->Hits[2][ihit].Position.Z(),hitcounter);
	    output().FillToyVectorVar(toy_hit_dqdx,(Float_t)part->Hits[2][ihit].dQdx_elife,hitcounter);
	    hitcounter++;
	  }
	  if(hitcounter>=1000)break;}
      }
    }
  }
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
