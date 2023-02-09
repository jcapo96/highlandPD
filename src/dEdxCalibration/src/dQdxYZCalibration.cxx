//header
#include "dQdxYZCalibration.hxx"

//system
#include <sys/time.h>
#include <iomanip>
#include <climits>

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

  _SelectedTracks = 0;
  _MaxTracks = INT_MAX; //arbitrarily large number

  h_global_yz     = NULL;
  h_global_yz_toy = NULL;
    
  for(int iz = 0; iz < nbinsz; iz++){
    for(int iy = 0; iy < nbinsy; iy++){
	h_local_yz[iz][iy]     = NULL;
	h_local_yz_toy[iz][iy] = NULL;
    }
  }
}

//********************************************************************
bool dQdxYZCalibration::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("dEdxCalibration.YZ.MinAccumLevelToSave"));

  _MaxTracks = ND::params().GetParameterI("dEdxCalibration.YZ.MaxTracks");

  _SaveAna = ND::params().GetParameterI("dEdxCalibration.YZ.SaveAna");
  _SaveToy = ND::params().GetParameterI("dEdxCalibration.YZ.SaveToy");

  _ApplySCEPositionCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCEPositionCorrection");
  _ApplySCEPitchCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCEPitchCorrection");
  _ApplyLifetimeCorrection = ND::params().GetParameterI("dEdxCalibration.YZ.ApplyLifetimeCorrection");

  _ApplySCESystematic = ND::params().GetParameterI("dEdxCalibration.YZ.ApplySCESystematic");
  _ApplyLifetimeSystematic = ND::params().GetParameterI("dEdxCalibration.YZ.ApplyLifetimeSystematic");

  //initialize histograms
  h_global_yz     = new TH1F("yz_global"    ,"yz_global"    ,              100,0,100);
  h_global_yz_toy = new TH2F("yz_global_toy","yz_global_toy",NTOYS,0,NTOYS,100,0,100);

  for(int iz = 0; iz < nbinsz; iz++){
    std::stringstream ssz;
    ssz << iz;
    for(int iy = 0; iy < nbinsy; iy++){
      std::stringstream ssy;
      ssy << iy;
      h_local_yz[iz][iy]     = new TH1F(("yz_local_"+ssz.str()+"_"+ssy.str()+"").c_str(),
					("yz_local_"+ssz.str()+"_"+ssy.str()+"").c_str(),
					100,0,100);
      h_local_yz_toy[iz][iy] = new TH2F(("yz_local_toy_"+ssz.str()+"_"+ssy.str()+"").c_str(),
					("yz_local_toy_"+ssz.str()+"_"+ssy.str()+"").c_str(),
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

  //nothing to save in the trees
}

//********************************************************************
void dQdxYZCalibration::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)

  //nothing to save in the trees
}

//********************************************************************
void dQdxYZCalibration::FillMicroTrees(bool addBase){
//********************************************************************
  
  if(_SelectedTracks > _MaxTracks)Finalize();

  //fill histograms withoug systematics
  int ntracks = box().Tracks.size();
  if(!(ntracks > 0))return;
  _SelectedTracks += ntracks;
  
  int zbin,ybin;
  for(int itrk = 0; itrk < ntracks; itrk++){
    AnaParticlePD* part = box().Tracks[itrk];
    int nhits = part->Hits[2].size();
    for(int ihit = 0; ihit < nhits; ihit++){
      if(!dEdxUtils::IsInterestingHit(part->Hits[2][ihit]))continue;
      h_global_yz->Fill(part->Hits[2][ihit].dQdx_elife);
      //fill local histogram
      zbin = part->Hits[2][ihit].Position.Z()/STEP;
      ybin = part->Hits[2][ihit].Position.Y()/STEP;
      h_local_yz[zbin][ybin]->Fill(part->Hits[2][ihit].dQdx_elife);
    }
  }
}

//********************************************************************
void dQdxYZCalibration::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  //skip 'ana' configuration, this is meant to work only with all_syst
  if(conf().GetCurrentConfigurationIndex() == 0)
    return;

  if(!_enableAllSystConfig)
    return;
  
  //fill histograms
  int ntracks = box().Tracks.size();
  if(!(ntracks > 0))return; 
  int itoy = conf().GetToyIndex();

  int zbin,ybin;
  for(int itrk = 0; itrk < ntracks; itrk++){
    AnaParticlePD* part = box().Tracks[itrk];
    int nhits = part->Hits[2].size();
    for(int ihit = 0; ihit < nhits; ihit++){
      if(!dEdxUtils::IsInterestingHit(part->Hits[2][ihit]))continue;
      //fill global histogram
      h_global_yz_toy->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife);
      //fill local histogram
      zbin = part->Hits[2][ihit].Position.Z()/STEP;
      ybin = part->Hits[2][ihit].Position.Y()/STEP;
      h_local_yz_toy[zbin][ybin]->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife);
    }
  }
}

//********************************************************************
void dQdxYZCalibration::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************
  
}

//********************************************************************
void dQdxYZCalibration::FillCategories(){
//********************************************************************

}

//********************************************************************
void dQdxYZCalibration::FinalizeAnalysisLoop(){
//********************************************************************
  
  baseAnalysis::FinalizeBunch();
  input().DeleteSpill();
  baseAnalysis::Finalize();
}

//********************************************************************
void dQdxYZCalibration::Finalize(){
//********************************************************************
  
  //speed up fitting procedure
  ROOT::EnableImplicitMT();
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);

  FillHistograms();
  if(_enableAllSystConfig)
    FillToyHistograms();

  output().CloseOutputFile();
  std::exit(1);
}

//********************************************************************
void dQdxYZCalibration::FillHistograms(){
//********************************************************************

  std::cout << "----------------------------------" << std::endl;
  std::cout << "Computing corrections"              << std::endl;
  std::cout << "----------------------------------" << std::endl;

  //compute means and erros, and fill final histograms
  TH2F* correction = new TH2F("correction_yz","correction_yz",nbinsz,ZMIN,ZMAX,nbinsy,YMIN,YMAX);
  
  //global values
  std::cout << "Fitting global value" << std::endl;

  //initialize clock
  timeval tim;
  gettimeofday(&tim, NULL);
  double t0=tim.tv_sec+(tim.tv_usec/1000000.0);
    
  TF1* f = new TF1("f",dEdxUtils::Langaus,0,100,4);
  double hmax = h_global_yz->GetBinCenter(h_global_yz->GetMaximumBin());
  dEdxUtils::SetFitParameters(f,hmax);
  h_global_yz->Fit("f","QNOR");
  double global_mpv       = f->GetParameter(1);
  double global_mpv_error = f->GetParError(1);
  double global_rel_error = global_mpv_error/global_mpv;
  delete h_global_yz;
  h_global_yz = NULL;

  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Global MPV = " << global_mpv << " +/- " << global_mpv_error << std::endl;
  std::cout << "Fit time = " << t1-t0 << std::endl;

  //local values
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Fitting local values" << std::endl;
  //loop over voxels
  for(int iz = 0; iz < nbinsz; iz++){
    for(int iy = 0; iy < nbinsy; iy++){
      double hmax = h_local_yz[iz][iy]->GetBinCenter(h_local_yz[iz][iy]->GetMaximumBin());
      dEdxUtils::SetFitParameters(f,hmax);
      h_local_yz[iz][iy]->Fit("f","QNOR");
      double local_mpv       = f->GetParameter(1);
      double local_mpv_error = f->GetParError(1);
      double local_rel_error = local_mpv_error/local_mpv;
      correction->SetBinContent(iz+1,iy+1,global_mpv/local_mpv);
      correction->SetBinError(iz+1,iy+1,global_mpv/local_mpv*sqrt(pow(global_rel_error,2)+pow(local_rel_error,2)));
      delete h_local_yz[iz][iy];
      h_local_yz[iz][iy] = NULL;
    }
    int per = (double)(iz+1)/nbinsz*100;
    if(per%10 == 0)std::cout << per << "% of local fits completed" << std::endl;
  }
  
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Fit time = " << t2-t1 << std::endl;
  
  //write histograms and close
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Writing correction histogram to root file" << std::endl;
  correction->Write();
}

//********************************************************************
void dQdxYZCalibration::FillToyHistograms(){
//********************************************************************

  std::cout << "----------------------------------" << std::endl;
  std::cout << "Computing toy corrections"          << std::endl;
  std::cout << "----------------------------------" << std::endl;
  
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
    hproj=h_global_yz_toy->ProjectionY("hproj",itoy+1,itoy+2);
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
  delete h_global_yz_toy;
  h_global_yz_toy = NULL;

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
	hproj=h_local_yz_toy[iz][iy]->ProjectionY("hproj",itoy+1,itoy+2);
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

      delete h_local_yz_toy[iz][iy];
      h_local_yz_toy[iz][iy] = NULL;
    }
    int per = (double)(iz+1)/nbinsz*100;
    if(per%10 == 0)std::cout << per << "% of local fits completed" << std::endl;
  }

  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Fit time = " << t2-t1 << std::endl;
  
  //write histograms and close
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Writing toy histograms to root file and finishing" << std::endl;
  error_map->Write();
  toy_correction->Write();
}
