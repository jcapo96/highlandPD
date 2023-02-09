#include "dQdxXCalibration.hxx"

#include <sys/time.h>
#include <iomanip>

#include "Math/MinimizerOptions.h"

#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"
#include "dEdxUtils.hxx"

#include "HighlandMiniTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"

#include "dEdxTrackSelection.hxx"

#include "HitPositionSCECorrection.hxx"
#include "HitPitchSCECorrection.hxx"
#include "SCECorrection.hxx"
#include "LifetimeCorrection.hxx"

#include "SCEVariation.hxx"
#include "LifetimeVariation.hxx"

//********************************************************************
dQdxXCalibration::dQdxXCalibration(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  h_global_x     = NULL;
  h_global_x_toy = NULL;
    
  for(int ix = 0; ix < nbinsx; ix++){
    h_local_x[ix]     = NULL;
    h_local_x_toy[ix] = NULL;
  }
  
  yz_correction     = NULL;
  yz_correction_toy = NULL;

  _sce = NULL;
}

//********************************************************************
bool dQdxXCalibration::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("dEdxCalibration.X.MinAccumLevelToSave"));

  _SaveAna = ND::params().GetParameterI("dEdxCalibration.X.SaveAna");
  _SaveToy = ND::params().GetParameterI("dEdxCalibration.X.SaveToy");

  _ApplySCEPositionCorrection = ND::params().GetParameterI("dEdxCalibration.X.ApplySCEPositionCorrection");
  _ApplySCEPitchCorrection = ND::params().GetParameterI("dEdxCalibration.X.ApplySCEPitchCorrection");
  _ApplyLifetimeCorrection = ND::params().GetParameterI("dEdxCalibration.X.ApplyLifetimeCorrection");

  _ApplySCESystematic = ND::params().GetParameterI("dEdxCalibration.X.ApplySCESystematic");
  _ApplyLifetimeSystematic = ND::params().GetParameterI("dEdxCalibration.X.ApplyLifetimeSystematic");

  //initialize histograms
  h_global_x     = new TH1F("x_global"    ,"x_global"    ,              100,0,100);
  h_global_x_toy = new TH2F("x_global_toy","x_global_toy",NTOYS,0,NTOYS,100,0,100);

  for(int ix = 0; ix < nbinsx; ix++){
    std::stringstream ssx;
    ssx << ix;
    h_local_x[ix] = new TH1F(("x_local_"+ssx.str()+"").c_str(),
			     ("x_local_"+ssx.str()+"").c_str(),
			     100,0,100);
    h_local_x_toy[ix] = new TH2F(("x_local_toy_"+ssx.str()+"").c_str(),
				 ("x_local_toy_"+ssx.str()+"").c_str(),
				 NTOYS,0,NTOYS,100,0,100);
  }
  
  //get yz correction histogram
  TFile* yzfile = TFile::Open((std::string(getenv("DEDXCALIBRATIONROOT"))+"/data/yz_lowerangles.root").c_str());
  yz_correction     = (TH2F*)yzfile->Get("correction_yz");
  yz_correction->SetDirectory(0);
  yz_correction_toy = (TH3F*)yzfile->Get("toy_correction_yz");
  if(yz_correction_toy)
    yz_correction_toy->SetDirectory(0);
  yzfile->Close();

  //always initialize space charge for ana configuration
  _sce = new SpaceCharge();
  _sce->Initialize();
  
  return true;
}

//********************************************************************
void dQdxXCalibration::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  _MichelRemovingTree = ND::params().GetParameterI("dEdxCalibration.X.MichelRemovingTree");
  std::stringstream treenumber;
  if(_MichelRemovingTree==0)treenumber.str("");
  else                      treenumber << _MichelRemovingTree;

  input().AddConverter("MichelRemovingTree", new MichelRemovingTreeConverter(("michelremoving"+treenumber.str()+"/Event").c_str()));
  input().AddConverter("minitreefiltered",   new HighlandMiniTreeConverter("MiniTree"));
}

//********************************************************************
void dQdxXCalibration::DefineSelections(){
//********************************************************************

  sel().AddSelection("non stopping crossing tracks", "non stopping crossing tracks", new dEdxTrackSelection(false));  
}

//********************************************************************
void dQdxXCalibration::DefineCorrections(){
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
void dQdxXCalibration::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  if(_ApplySCESystematic)     evar().AddEventVariation(kSCE_syst     , "SCE variation"     , new SCEVariation());
  if(_ApplyLifetimeSystematic)evar().AddEventVariation(kLifetime_syst, "lifetime variation", new LifetimeVariation());
}

//********************************************************************
void dQdxXCalibration::DefineConfigurations(){
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
void dQdxXCalibration::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------
}

//********************************************************************
void dQdxXCalibration::DefineTruthTree(){
//********************************************************************

}

//********************************************************************
void dQdxXCalibration::FillMicroTrees(bool addBase){
//********************************************************************

  
  
  //fill histograms
  int ntracks = box().Tracks.size();
  if(!(ntracks > 0))return; 
  
  int zbin,ybin,xbin;
  for(int itrk = 0; itrk < ntracks; itrk++){
    AnaParticlePD* part = box().Tracks[itrk];
    int nhits = part->Hits[2].size();
    for(int ihit = 0; ihit < nhits; ihit++){
      if(!dEdxUtils::IsInterestingHit(part->Hits[2][ihit]))continue;
      //get corrections
      zbin = part->Hits[2][ihit].Position.Z()/STEP;
      ybin = part->Hits[2][ihit].Position.Y()/STEP;
      double yzcorr = yz_correction->GetBinContent(zbin+1,ybin+1);
      double efield = dEdxUtils::ElectricField(GetSCE(),
					       part->Hits[2][ihit].Position.X(),
					       part->Hits[2][ihit].Position.Y(),
					       part->Hits[2][ihit].Position.Z());
      double recfact = dEdxUtils::Recombination(efield);
      //fill global histogram
      h_global_x->Fill(part->Hits[2][ihit].dQdx_elife*yzcorr*recfact);
      //fill local histogram
      xbin = part->Hits[2][ihit].Position.X()/STEP+nbinsx;
      h_local_x[xbin]->Fill(part->Hits[2][ihit].dQdx_elife*yzcorr*recfact);
    }
  }
}

//********************************************************************
void dQdxXCalibration::FillToyVarsInMicroTrees(bool addBase){
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
  
  int zbin,ybin,xbin;
  for(int itrk = 0; itrk < ntracks; itrk++){
    AnaParticlePD* part = box().Tracks[itrk];
    int nhits = part->Hits[2].size();
    for(int ihit = 0; ihit < nhits; ihit++){
      if(!dEdxUtils::IsInterestingHit(part->Hits[2][ihit]))continue;
      //get corrections
      zbin = part->Hits[2][ihit].Position.Z()/STEP;
      ybin = part->Hits[2][ihit].Position.Y()/STEP;
      double yzcorr = yz_correction->GetBinContent(zbin+1,ybin+1,itoy+1);
      double efield = dEdxUtils::ElectricField(GetSCE(),
					       part->Hits[2][ihit].Position.X(),
					       part->Hits[2][ihit].Position.Y(),
					       part->Hits[2][ihit].Position.Z());
      double recfact = dEdxUtils::Recombination(efield);
      //fill global histogram
      h_global_x_toy->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife*yzcorr*recfact);
      //fill local histogram
      xbin = part->Hits[2][ihit].Position.X()/STEP+nbinsx;
      h_local_x_toy[xbin]->Fill(itoy+1,part->Hits[2][ihit].dQdx_elife*yzcorr*recfact);
    }
  }
}

//********************************************************************
SpaceCharge* dQdxXCalibration::GetSCE(){
//********************************************************************

  if(conf().GetCurrentConfigurationIndex() == 0)return _sce;
  else{
    if(_ApplySCESystematic){
      Int_t itoy = conf().GetToyIndex();
      return static_cast<SCEVariation*>(evar().GetEventVariation("SCE variation"))->GetToySCE(itoy);
    }
    else return _sce;
  }
}

//********************************************************************
void dQdxXCalibration::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

}

//********************************************************************
void dQdxXCalibration::FillCategories(){
//********************************************************************

}

//********************************************************************
void dQdxXCalibration::Finalize(){
//********************************************************************
  
  //speed up fitting procedure
  ROOT::EnableImplicitMT();
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);

  FillHistograms();
  if(_enableAllSystConfig)
    FillToyHistograms();
}

//********************************************************************
void dQdxXCalibration::FillHistograms(){
//********************************************************************

  std::cout << "----------------------------------" << std::endl;
  std::cout << "Computing corrections             " << std::endl;
  std::cout << "----------------------------------" << std::endl;

  //compute means and erros, and fill final histograms
  TH1F* correction = new TH1F("correction_x","correction_x",nbinsx,XMIN,XMAX);
  
  //global values
  std::cout << "Fitting global value" << std::endl;
  
  //initialize clock
  timeval tim;
  gettimeofday(&tim, NULL);
  double t0=tim.tv_sec+(tim.tv_usec/1000000.0);
  
  TF1* f = new TF1("f",dEdxUtils::Langaus,0,100,4);
  double hmax = h_global_x->GetBinCenter(h_global_x->GetMaximumBin());
  dEdxUtils::SetFitParameters(f,hmax);
  h_global_x->Fit("f","QNOR");
  double global_mpv = f->GetParameter(1);
  double global_mpv_error = f->GetParError(1);
  double global_rel_error = global_mpv_error/global_mpv;
  delete h_global_x;

  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Global MPV = " << global_mpv << " +/- " << global_mpv_error << std::endl;
  std::cout << "Fit time = " << t1-t0 << std::endl;

  //local values
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Fitting local values" << std::endl;
  //loop over voxels
  for(int ix = 0; ix < nbinsx; ix++){
    double hmax = h_local_x[ix]->GetBinCenter(h_local_x[ix]->GetMaximumBin());
    dEdxUtils::SetFitParameters(f,hmax);
    h_local_x[ix]->Fit("f","QNOR");
    double local_mpv = f->GetParameter(1);
    double local_mpv_error = f->GetParError(1);
    double local_rel_error = local_mpv_error/local_mpv;
    correction->SetBinContent(ix+1,global_mpv/local_mpv);
    correction->SetBinError(ix+1,global_mpv/local_mpv*sqrt(pow(global_rel_error,2)+pow(local_rel_error,2)));
    delete h_local_x[ix];
  }
  
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Fit time = " << t2-t1 << std::endl;
  
  //write histograms and close
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Writing correction histogram to root" << std::endl;
  correction->Write();
}

//********************************************************************
void dQdxXCalibration::FillToyHistograms(){
//********************************************************************

  std::cout << "----------------------------------" << std::endl;
  std::cout << "Computing toy corrections         " << std::endl;
  std::cout << "----------------------------------" << std::endl;

  //compute means and erros, and fill final histograms
  TH1F* error_map = new TH1F("error_map_x","error_map_x",nbinsx,XMIN,XMAX);
  TH2F* toy_correction = new TH2F("toy_correction_x","toy_correction_x",nbinsx,XMIN,XMAX,NTOYS,0,NTOYS);
  
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
  f->SetParameters(1,60,5000,5);
  f->SetParLimits(0,1,5);
  f->SetParLimits(3,1,5);
  TH1D* hproj = new TH1D("hproj","hproj",100,0,100);
  for(int itoy = 0; itoy < NTOYS; itoy++){
    hproj=h_global_x_toy->ProjectionY("hproj",itoy+1,itoy+2);
    double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
    f->SetParameter(1,hmax);
    f->SetParLimits(1,hmax-5,hmax+10);
    f->SetRange(hmax-7.5,hmax+10);
    hproj->Fit("f","QNOR");
    hsyst->Fill(f->GetParameter(1));
    global_mpv_vector.push_back(f->GetParameter(1));
    hproj->Reset();
  }
  global_mpv = hsyst->GetMean();
  global_syst_error = hsyst->GetRMS();
  hsyst->Reset();
  delete h_global_x_toy;

  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  std::cout << "Global MPV = " << global_mpv << " +/- " << global_syst_error << std::endl;
  std::cout << "Fit time = " << t1-t0 << std::endl;

  //local values
  std::cout << "----------------------------------" << std::endl;
  std::cout << "Fitting local values" << std::endl;
  double local_mpv, local_syst_error;
  //loop over voxels
  for(int ix = 0; ix < nbinsx; ix++){
    //loop over toys
    for(int itoy = 0; itoy < NTOYS; itoy++){
      hproj=h_local_x_toy[ix]->ProjectionY("hproj",itoy+1,itoy+2);
      double hmax = hproj->GetBinCenter(hproj->GetMaximumBin());
      f->SetParameter(1,hmax);
      f->SetParLimits(1,hmax-5,hmax+10);
      f->SetRange(hmax-7.5,hmax+10);
      hproj->Fit("f","QNOR");
      hsyst->Fill(f->GetParameter(1));
      toy_correction->SetBinContent(ix+1,itoy+1,global_mpv_vector[itoy]/f->GetParameter(1));
      hproj->Reset();
    }
    local_mpv = hsyst->GetMean();
    local_syst_error = hsyst->GetRMS();
    hsyst->Reset();
    error_map->SetBinContent(ix+1,sqrt(pow(global_syst_error/global_mpv,2)+pow(local_syst_error/local_mpv,2))*100);

    delete h_local_x_toy[ix];
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
