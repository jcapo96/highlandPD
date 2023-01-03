//header
#include "dQdxYZCalibration.hxx"

//system
#include <sys/time.h>
#include <iomanip>

//root
#include "Math/MinimizerOptions.h"

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
dQdxYZCalibration::dQdxYZCalibration(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  h_global_yz = NULL;
    
  for(int iz = 0; iz < nbinsz; iz++)
      for(int iy = 0; iy < nbinsy; iy++)
	h_local_yz[iz][iy] = NULL;  
}

//********************************************************************
bool dQdxYZCalibration::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("dEdxCalibration.YZ.MinAccumLevelToSave"));

  _SaveAna = ND::params().GetParameterI("dEdxCalibration.YZ.SaveAna");
  _SaveToy = ND::params().GetParameterI("dEdxCalibration.YZ.SaveToy");

  _ApplySCEPositionCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCEPositionCorrection");
  _ApplySCEPitchCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCEPitchCorrection");
  _ApplyLifetimeCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplyLifetimeCorrection");

  _ApplySCESystematic = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCESystematic");
  _ApplyLifetimeSystematic = ND::params().GetParameterI("dEdxCalibration.YZ.ApplyLifetimeSystematic");

  TTree::SetMaxTreeSize(200000000000);

  //initialize histograms
  h_global_yz = new TH2F("yz_global","yz_global",NTOYS,0,NTOYS,100,0,100);

  for(int iz = 0; iz < nbinsz; iz++){
    std::stringstream ssz;
    ssz << iz;
    for(int iy = 0; iy < nbinsy; iy++){
      std::stringstream ssy;
      ssy << iy;
      h_local_yz[iz][iy] = new TH2F(("yz_local_"+ssz.str()+"_"+ssy.str()+"").c_str(),
				    ("yz_local_"+ssz.str()+"_"+ssy.str()+"").c_str(),
				    NTOYS,0,NTOYS,100,0,100);
    }
  }

  return true;
}

//********************************************************************
void dQdxYZCalibration::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  _MichelRemovingTree = ND::params().GetParameterI("dEdxCalibration.YZ.MichelRemovingTree");
  std::stringstream treenumber;
  if(_MichelRemovingTree==0)treenumber.str("");
  else                      treenumber << _MichelRemovingTree;

  input().AddConverter("MichelRemovingTree", new MichelRemovingTreeConverter(("michelremoving"+treenumber.str()+"/Event").c_str()));
  input().AddConverter("minitreefiltered",   new HighlandMiniTreeConverter("MiniTree"));
}

//********************************************************************
void dQdxYZCalibration::DefineSelections(){
//********************************************************************

  sel().AddSelection("non stopping crossing tracks", "non stopping crossing tracks", new dEdxTrackSelection(false));  
}

//********************************************************************
void dQdxYZCalibration::DefineCorrections(){
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
void dQdxYZCalibration::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  if(_ApplySCESystematic)     evar().AddEventVariation(kSCE_syst     , "SCE variation"     , new SCEVariation());
  if(_ApplyLifetimeSystematic)evar().AddEventVariation(kLifetime_syst, "lifetime variation", new LifetimeVariation());
}

//********************************************************************
void dQdxYZCalibration::DefineConfigurations(){
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
void dQdxYZCalibration::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);
  
  if(_SaveAna){
    AddVarMF(output(), track_hit_x   , "trk hit x"   , ntracks, -30, 3000);
    AddVarMF(output(), track_hit_y   , "trk hit y"   , ntracks, -30, 3000);
    AddVarMF(output(), track_hit_z   , "trk hit z"   , ntracks, -30, 3000);
    AddVarMF(output(), track_hit_dqdx, "trk hit dqdx", ntracks, -30, 3000);
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
void dQdxYZCalibration::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void dQdxYZCalibration::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // fill trk info
  if(_SaveAna){
    if(box().Tracks.size()>0){
      for(int i = 0; i < std::min((int)box().Tracks.size(),30); i++){
	AnaParticlePD* part = box().Tracks[i];
	for(int j = 0; j < (int)std::min((int)part->Hits[2].size(),3000); j++){
	  output().FillMatrixVar(track_hit_x,         (Float_t)part->Hits[2][j].Position.X(),   -1, j);
	  output().FillMatrixVar(track_hit_y,         (Float_t)part->Hits[2][j].Position.Y(),   -1, j);
	  output().FillMatrixVar(track_hit_z,         (Float_t)part->Hits[2][j].Position.Z(),   -1, j);
	  output().FillMatrixVar(track_hit_dqdx,      (Float_t)part->Hits[2][j].dQdx_elife,     -1, j);
	}
	output().IncrementCounter(ntracks);
      } 
    }
  }
}

//********************************************************************
void dQdxYZCalibration::FillToyVarsInMicroTrees(bool addBase){
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

  //fill histograms
  int ntracks = box().Tracks.size();
  if(!(ntracks > 0))return; 
  int itoy = conf().GetToyIndex();

  int zbin,ybin;
  for(int itrk = 0; itrk < ntracks; itrk++){
    AnaParticlePD* part = box().Tracks[itrk];
    int nhits = part->Hits[2].size();
    for(int ihit = 0; ihit < nhits; ihit++){
      if(dEdxUtils::IsInterestingHit(part->Hits[2][ihit])){
	//fill global histogram
	h_global_yz->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife);
	//fill local histogram
	zbin = part->Hits[2][ihit].Position.Z()/STEP;
	ybin = part->Hits[2][ihit].Position.Y()/STEP;
	h_local_yz[zbin][ybin]->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife);
      }
    }
  }
}

//********************************************************************
void dQdxYZCalibration::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void dQdxYZCalibration::FillCategories(){
//********************************************************************

}

//********************************************************************
void dQdxYZCalibration::Finalize(){
//********************************************************************
  
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Computing error histograms and toy corrections" << std::endl;
  std::cout << "----------------------------------" << std::endl;

  //speed up fitting procedure
  ROOT::EnableImplicitMT();
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);
  
  //compute means and erros, and fill final histograms
  TH2F* error_map = new TH2F("error_map_yz","error_map_yz",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX);
  TH3F* toy_correction = new TH3F("toy_correction_yz","toy_correction_yz",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX,NTOYS,0,NTOYS);
  
  //global values
  std::cout << "Fitting global values" << std::endl;

  //initialize clock
  timeval tim;
  gettimeofday(&tim, NULL);
  double t0=tim.tv_sec+(tim.tv_usec/1000000.0);
    
  TH1F* hsyst = new TH1F("hsyst","hsyst",1000,0,100);
  double global_mpv, global_syst_error;
  std::vector<double> global_mpv_vector;
  global_mpv_vector.clear();
  
  TF1* f = new TF1("f",dEdxUtils::Langaus,0,100,4);
  
  TH1D* hproj = new TH1D("hproj","hproj",100,0,100);
  for(int itoy = 0; itoy < NTOYS; itoy++){
    hproj=h_global_yz->ProjectionY("hproj",itoy+1,itoy+2);
    double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
    dEdxUtils::SetFitParameters(f,hmax);
    hproj->Fit("f","QNOR");
    hsyst->Fill(f->GetParameter(1));
    global_mpv_vector.push_back(f->GetParameter(1));
    hproj->Reset();
  }
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  hsyst->Reset();
  delete h_global_yz;

  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Global MPV = " << global_mpv << " +/- " << global_syst_error << std::endl;
  std::cout << "Fit time = " << t1-t0 << std::endl;

  //local values
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Fitting local values" << std::endl;
  double local_mpv, local_syst_error;
  //loop over voxels
  for(int iz = 0; iz < nbinsz; iz++){
    for(int iy = 0; iy < nbinsy; iy++){
      //loop over toys
      for(int itoy = 0; itoy < NTOYS; itoy++){
	hproj=h_local_yz[iz][iy]->ProjectionY("hproj",itoy+1,itoy+2);
	double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
	dEdxUtils::SetFitParameters(f,hmax);
	hproj->Fit("f","QNOR");
	hsyst->Fill(f->GetParameter(1));
	toy_correction->SetBinContent(iz+1,iy+1,itoy+1,global_mpv_vector[itoy]/f->GetParameter(1));
	hproj->Reset();
      }
      local_mpv = hsyst->GetMean();
      local_syst_error = hsyst->GetRMS();
      hsyst->Reset();
      error_map->SetBinContent(iz+1,iy+1,sqrt(pow(global_syst_error/global_mpv,2)+pow(local_syst_error/local_mpv,2))*100);

      delete h_local_yz[iz][iy];
    }
    std::cout << (double)(iz+1)/nbinsz*100 << "% of local fits completed" << std::setprecision(4) << std::endl;
  }

  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Fit time = " << t2-t1 << std::endl;
  
  //write histograms and close
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Writing final histograms to root file and finishing" << std::endl;
  error_map->Write();
  toy_correction->Write();
}
