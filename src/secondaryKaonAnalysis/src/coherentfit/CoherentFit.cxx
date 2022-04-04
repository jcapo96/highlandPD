#include "CoherentFit.hxx"
#include "CoherentFitUtils.hxx"

#include "TPad.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include <iostream>

//********************************************************************
CoherentFit::CoherentFit(){
//********************************************************************

  fFilename = "";
  fFile = NULL;  
  fMiniTree = NULL;
  fSpill = NULL;
  fSignalPlusBackground = NULL;
  fSignal = NULL;
  fBackground = NULL;
  fTrueSignal = NULL;
  fTrueBackground = NULL;
}

//********************************************************************
CoherentFit::CoherentFit(const std::string& filename){
//********************************************************************
  
  //get MiniTree
  fFilename = filename;
  fFile = TFile::Open(fFilename.c_str());
  
  if(!fFile){
    std::cout << "couldn't find " << fFilename << std::endl;
    std::exit(1);
  }
  
  fMiniTree = (TTree*)fFile->Get("MiniTree");
  fSpill = new AnaSpillB();
  fMiniTree->SetBranchAddress("Spill", &fSpill);
  
  fSignalPlusBackground = NULL;
  fSignal = NULL;
  fBackground = NULL;
  fTrueSignal = NULL;
  fTrueBackground = NULL;
}

//********************************************************************
void CoherentFit::WriteToRootFile(const std::string& filename){
//********************************************************************

  fSignalPlusBackground->WriteToRootFile(filename+"_SB");
  fSignal->WriteToRootFile(filename+"_S");
  fBackground->WriteToRootFile(filename+"_B");
  fTrueSignal->WriteToRootFile(filename+"_TS");
  fTrueBackground->WriteToRootFile(filename+"_TB");
  
}
  
//********************************************************************
void CoherentFit::CreateCoherentSamples(const Double_t Chi2Cut){
//********************************************************************
  
  fSignalPlusBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kSignalPlusBackground);
  fSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kSignal);
  fBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kBackground);
  fTrueSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSignal);
  fTrueBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueBackground);
  
  fSignalPlusBackground->SetSignal(fSignal);
  fSignalPlusBackground->SetBackground(fBackground);
  fSignal->SetBackground(fBackground);
  fBackground->SetSignal(fSignal);
  fSignalPlusBackground->SetTrueSignal(fTrueSignal);
  fSignalPlusBackground->SetTrueBackground(fTrueBackground);
  fTrueSignal->SetTrueBackground(fTrueBackground);
  fTrueBackground->SetTrueSignal(fTrueSignal);  
  
  fSignalPlusBackground->SetChi2Cut(Chi2Cut);
  fSignal->SetChi2Cut(Chi2Cut);
  fBackground->SetChi2Cut(Chi2Cut);
  fSignalPlusBackground->SetChi2Cut(Chi2Cut);
  fTrueSignal->SetChi2Cut(Chi2Cut);
  fTrueBackground->SetChi2Cut(Chi2Cut);  
}

//********************************************************************
void CoherentFit::SetParametersFromMCFit(const CoherentFit* c){
//********************************************************************

  fSignal->SetCFitParameters(c->GetSignalSample());
  fBackground->SetCFitParameters(c->GetBackgroundSample());
}
  
//********************************************************************
void CoherentFit::GenerateTrueMCHistograms(const double RMIN, const double RMAX, const double STEP,
					   const double Chi2Cut,
					   const double bin_min, const double bin_max, const double bin_width,
					   const bool normalize, const bool resize){
//********************************************************************

  fMiniTree->GetEntry(0);
  if(!fSpill->GetIsMC()){
    std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
    std::exit(1);
  }
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  
  while(rmax <= RMAX){
    fSignalPlusBackground->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    fSignalPlusBackground->AddToHistVector(CoherentFitUtils::GetHistogramFromResRangeSlice(fMiniTree, fSpill, rmin, rmax, Chi2Cut, bin_min, bin_max, bin_width));
    fTrueSignal->AddToHistVector(CoherentFitUtils::GetSignalHistogramFromResRangeSlice(fMiniTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
    fTrueBackground->AddToHistVector(CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(fMiniTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
    rmin += STEP;
    rmax += STEP;
  }

  fTrueSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fTrueBackground->SetRRVector(fSignalPlusBackground->GetRRVector());
  fSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fBackground->SetRRVector(fSignalPlusBackground->GetRRVector());
  
  if(normalize)NormalizeHistograms();
  //if(resize)ResizeHistograms();
}

//********************************************************************
void CoherentFit::GenerateTrueBackgroundHistogramsFromTrueSignal(double RMAX){
//********************************************************************

  fMiniTree->GetEntry(0);
  if(!fSpill->GetIsMC()){
    std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
    std::exit(1);
  }

  if(!fTrueSignal){
    std::cout << "no true signal sample exists!" << std::endl;
    std::exit(1);
  }
    
  if(fTrueSignal->GetHistVector().size() < 2){
    std::cout << "no true signal histograms exist!" << std::endl;
    std::exit(1);
  } 

  double STEP = fTrueSignal->GetRRVector()[1].first-fTrueSignal->GetRRVector()[0].first;
  double rmin = fTrueSignal->GetRRVector()[0].first - STEP/2;
  double rmax = fTrueSignal->GetRRVector()[0].first + STEP/2;
  if(RMAX == -1)RMAX = fTrueSignal->GetRRVector().back().first + STEP/2;
  double Chi2Cut = fTrueSignal->GetChi2Cut();
  
  int counter = 0;
  while(rmax <= RMAX && counter < (int)fTrueSignal->GetHistVector().size()){
    fTrueBackground->AddToHistVector(CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(fMiniTree, fSpill, fTrueSignal->GetHistVector()[counter], rmin, rmax, Chi2Cut));
    rmin += STEP;
    rmax += STEP;
    counter++;
  }

  fTrueBackground->SetRRVector(fTrueSignal->GetRRVector());
  fTrueBackground->NormalizeHistograms(fTrueSignal->GetIntegralVector());
}

//********************************************************************
void CoherentFit::GenerateFakeBackground(const double RMIN, const double RMAX, const double STEP,
					 const double Chi2Cut,
					 const double shift_mean, const double shift_sigma){
//********************************************************************
  
  fMiniTree->GetEntry(0);
  if(!fSpill->GetIsMC()){
    std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
    std::exit(1);
  }

  //randomizer
  TRandom3* r  = new TRandom3();
  r->SetSeed(1);
  double shift = 0;
  //generate a shift for each particle, which should be used for filling every histogram
  std::vector<double> shift_vector;
  shift_vector.clear();
  for(int ientry = 0; ientry < fMiniTree->GetEntries(); ientry++){
    fMiniTree->GetEntry(ientry);
    AnaBunchB* bunch = static_cast<AnaBunchB*>(fSpill->Bunches[0]);
    for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
      //while(shift <= 0)
      shift = r->Gaus(shift_mean,shift_sigma);
      if(shift < 0)shift = 0;
      shift_vector.push_back(shift);
      shift = 0;
    }
  }
  //reset spill to entry zero
  fMiniTree->GetEntry(0);
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  int counter = 0;
  
  while(rmax <= RMAX){
    fTrueBackground->AddToHistVector(CoherentFitUtils::GenerateBackgroundHistogramFromTrueSignal(fMiniTree, fSpill, fTrueSignal->GetHistVector()[counter], rmin, rmax, Chi2Cut, shift_vector));
    counter++;
    rmin += STEP;
    rmax += STEP;
  }
  
  fTrueBackground->SetRRVector(fTrueSignal->GetRRVector());
  fTrueBackground->NormalizeHistograms(fTrueSignal->GetIntegralVector());
}

//********************************************************************
void CoherentFit::GenerateHistograms(const double RMIN, const double RMAX, const double STEP,
				     const double Chi2Cut,
				     const double bin_min, const double bin_max, const double bin_width,
				     const bool normalize, const bool resize){
//********************************************************************
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  
  while(rmax <= RMAX){
    fSignalPlusBackground->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    fSignalPlusBackground->AddToHistVector(CoherentFitUtils::GetHistogramFromResRangeSlice(fMiniTree, fSpill, rmin, rmax, Chi2Cut, bin_min, bin_max, bin_width));
    rmin += STEP;
    rmax += STEP;
  }

  fSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fBackground->SetRRVector(fSignalPlusBackground->GetRRVector());
  
  if(normalize)NormalizeHistograms();
  //if(resize)ResizeHistograms();
}

//********************************************************************
void CoherentFit::NormalizeHistograms(){
//********************************************************************

  fSignalPlusBackground->NormalizeHistograms();
  if(GetIsMC()){
    fTrueSignal->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
    fTrueBackground->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
  }
}

//********************************************************************
void CoherentFit::SequentialCoherentFit(){
//********************************************************************
  
  fTrueSignal->SequentialCoherentFit();
  fTrueBackground->SequentialCoherentFit();
  fSignalPlusBackground->SequentialCoherentFit();
}

//********************************************************************
void CoherentFit::DataCoherentFit(const CoherentFit* c){
//********************************************************************
  
  SetParametersFromMCFit(c);
  fSignalPlusBackground->CoherentFit();
}
