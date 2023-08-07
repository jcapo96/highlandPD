#include "CoherentFit.hxx"
#include "CoherentFitUtils.hxx"
#include "TPad.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include <iostream>

//********************************************************************
CoherentFit::CoherentFit(){
//********************************************************************

  fFilename = "";
  fFile = NULL;  
  fTree = NULL;
  fTreeSystematics = NULL;
  fIsMiniTree = false;
  fIsSystTree = false;
  fIsMC = false;
  fSpill = NULL;
  fSignalPlusBackground = NULL;
  fSignal = NULL;
  fSemiSignal = NULL;
  fBackground = NULL;
  fSemiBackground = NULL;
  fTrueSignal = NULL;
  fTrueSemiSignal = NULL;
  fTrueBackground = NULL;
  fTrueSemiBackground = NULL;
  fToySample = NULL;
  h_toy_A = NULL;
  h_toy_B = NULL;
  h_toy_C = NULL;
  h_toy_Eff = NULL;
  h_toy_dEdx_RR = NULL;
}

//********************************************************************
CoherentFit::CoherentFit(const std::string& filename, bool IsMC){
//********************************************************************
  
  //get MiniTree
  fFilename = filename;
  fFile = TFile::Open(fFilename.c_str());
  
  if(!fFile){
    std::cout << "couldn't find " << fFilename << std::endl;
    std::exit(1);
  }

  fIsMC = IsMC;
  
  fTree = GetTreeFromRootFile();//(TTree*)fFile->Get("MiniTree");
  if(fIsMiniTree){
    fSpill = new AnaSpillB();
    fTree->SetBranchAddress("Spill", &fSpill);
    fTree->GetEntry(0);
    fIsMC = fSpill->GetIsMC();
  }
  else fSpill = NULL;

  fSignalPlusBackground = NULL;
  fSignal = NULL;
  fSemiSignal = NULL;
  fBackground = NULL;
  fSemiBackground = NULL;
  fTrueSignal = NULL;
  fTrueSemiSignal = NULL;
  fTrueBackground = NULL;
  fTrueSemiBackground = NULL;
  fToySample = NULL;

  h_toy_A = NULL;
  h_toy_B = NULL;
  h_toy_C = NULL;
  h_toy_Eff = NULL;
  h_toy_dEdx_RR = NULL;
}

//********************************************************************
TTree* CoherentFit::GetTreeFromRootFile(){
//********************************************************************

  TTree* t = (TTree*)fFile->Get("MiniTree");
  if(t){
    fIsMiniTree = true;
    fIsSystTree = false;
    std::cout << "Input tree is MiniTree" << std::endl;
    return t;
  }

  t = (TTree*)fFile->Get("ana");
  //t = (TTree*)fFile->Get("dQdx_YZcal");
  if(t){
    fIsMiniTree = false;
    fIsSystTree = false;
    std::cout << "Input tree is flat tree" << std::endl;
    return t;
  }

  //if here, no known type of tree provided
  std::cout << "ERROR! TTree type not known. Please check your input file" << std::endl;
  std::exit(1);
}

//********************************************************************
TTree* CoherentFit::GetSystTree(){
//********************************************************************

  TFile* f = TFile::Open(fFilename.c_str());
  
  TTree* t = (TTree*)fFile->Get("all_syst");
  if(t){
    fIsMiniTree = false;
    fIsSystTree = true;
    std::cout << "Tree with systematics found" << std::endl;
    f->Close();
    return t;
  }

  //if here, no known type of tree provided
  std::cout << "ERROR! No TTree with systematics found!" << std::endl;
  std::exit(1);
}

//********************************************************************
void CoherentFit::WriteToRootFile(const std::string& filename){
//********************************************************************

  fSignalPlusBackground->WriteToRootFile(filename+"_SB.root");
  fSignal->WriteToRootFile(filename+"_S.root");
  fBackground->WriteToRootFile(filename+"_B.root");
  fTrueSignal->WriteToRootFile(filename+"_TS.root");
  fTrueBackground->WriteToRootFile(filename+"_TB.root");
  
}
  
//********************************************************************
void CoherentFit::CreateCoherentSamples(const Double_t Chi2Cut){
//********************************************************************
  
  fSignalPlusBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kSignalPlusBackground);
  fSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kSignal);
  fSemiSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kSemiSignal);
  fBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kBackground);
  fSemiBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kSemiBackground);
  fTrueSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSignal);
  fTrueSemiSignal = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSemiSignal);
  fTrueBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueBackground);
  fTrueSemiBackground = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSemiBackground);
  
  CreateSampleLinks();
  
  fSignalPlusBackground->SetChi2Cut(Chi2Cut);
  fSignal->SetChi2Cut(Chi2Cut);
  fBackground->SetChi2Cut(Chi2Cut);
  fSignalPlusBackground->SetChi2Cut(Chi2Cut);
  fTrueSignal->SetChi2Cut(Chi2Cut);
  fTrueBackground->SetChi2Cut(Chi2Cut);  
}

//********************************************************************
void CoherentFit::CreateSampleLinks(){
//********************************************************************
  
  fSignalPlusBackground->SetSignal(fSignal);
  fSignalPlusBackground->SetBackground(fBackground);
  fSignal->SetSemiSignal(fSemiSignal);
  fSignal->SetBackground(fBackground);
  fSemiSignal->SetSignal(fSignal);
  fBackground->SetSignal(fSignal);
  fBackground->SetSemiBackground(fSemiBackground);
  fSemiBackground->SetBackground(fBackground);
  fSignalPlusBackground->SetTrueSignal(fTrueSignal);
  fSignalPlusBackground->SetTrueBackground(fTrueBackground);
  fTrueSignal->SetTrueSemiSignal(fTrueSemiSignal);
  fTrueSignal->SetTrueBackground(fTrueBackground);
  fTrueSemiSignal->SetTrueSignal(fTrueSignal);
  fTrueBackground->SetTrueSignal(fTrueSignal);
  fTrueBackground->SetTrueSemiBackground(fTrueSemiBackground);
  fTrueSemiBackground->SetTrueBackground(fTrueBackground);
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

  if(fIsMiniTree){
    fTree->GetEntry(0);
    if(!fSpill->GetIsMC()){
      std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
      std::exit(1);
    }
  }

  // double min[15] = {5 ,4.5,4.5,4,3.5,3.5,3  ,3,3,2.5,2.5,2.5,2  ,2  ,2};
  // double max[15] = {14,10 ,8  ,7,6.5,6  ,5.5,5,5,5  ,5  ,4.5,4.5,4.5,4.5};
  // int counter = 0;
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;

  while(rmax <= RMAX){
    fSignalPlusBackground->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    fSignalPlusBackground->AddToHistVector(CoherentFitUtils::GetHistogramFromResRangeSlice(fTree, fSpill, rmin, rmax, Chi2Cut, bin_min, bin_max, bin_width));
    fTrueSignal->AddToHistVector(CoherentFitUtils::GetSignalHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
    fTrueBackground->AddToHistVector(CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
    rmin += STEP;
    rmax += STEP;
    // counter++;
  }
  
  fTrueSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fTrueSemiSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fTrueBackground->SetRRVector(fSignalPlusBackground->GetRRVector());
  fSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fBackground->SetRRVector(fSignalPlusBackground->GetRRVector());

  for(int ihist = 0; ihist < (int)fSignalPlusBackground->GetHistVector().size(); ihist++){
    fSignalPlusBackground->GetHistVector()[ihist]->SetBinErrorOption(TH1::kPoisson);
    fSignalPlusBackground->GetHistVector()[ihist]->Sumw2(kTRUE);
    fTrueSignal->GetHistVector()[ihist]->SetBinErrorOption(TH1::kPoisson);
    fTrueSignal->GetHistVector()[ihist]->Sumw2(kTRUE);
    fTrueBackground->GetHistVector()[ihist]->SetBinErrorOption(TH1::kPoisson);
    fTrueBackground->GetHistVector()[ihist]->Sumw2(kTRUE);
  }
  
  if(normalize)NormalizeHistograms();
}

//********************************************************************
void CoherentFit::GenerateTrueBackgroundHistogramsFromTrueSignal(double RMAX){
//********************************************************************

  fTree->GetEntry(0);
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
    fTrueBackground->AddToHistVector(CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(fTree, fSpill, fTrueSignal->GetHistVector()[counter], rmin, rmax, Chi2Cut));
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
  
  fTree->GetEntry(0);
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
  for(int ientry = 0; ientry < fTree->GetEntries(); ientry++){
    fTree->GetEntry(ientry);
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
  fTree->GetEntry(0);
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  int counter = 0;
  
  while(rmax <= RMAX){
    fTrueBackground->AddToHistVector(CoherentFitUtils::GenerateBackgroundHistogramFromTrueSignal(fTree, fSpill, fTrueSignal->GetHistVector()[counter], rmin, rmax, Chi2Cut, shift_vector));
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
    fSignalPlusBackground->AddToHistVector(CoherentFitUtils::GetHistogramFromResRangeSlice(fTree, fSpill, rmin, rmax, Chi2Cut, bin_min, bin_max, bin_width));
    rmin += STEP;
    rmax += STEP;
  }

  fSignal->SetRRVector(fSignalPlusBackground->GetRRVector());
  fBackground->SetRRVector(fSignalPlusBackground->GetRRVector());
  
  if(normalize)NormalizeHistograms();
  //if(resize)ResizeHistograms();
}

//********************************************************************
CoherentSample* CoherentFit::CreateTrueStoppingKaonsSample(const double RMIN, const double RMAX, const double STEP,
							   const double Chi2Cut,
							   const double bin_min, const double bin_max, const double bin_width,
							   const bool normalize, const bool resize) const {
//********************************************************************

  if(fIsMiniTree){
    fTree->GetEntry(0);
    if(!fSpill->GetIsMC()){
      std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
      std::exit(1);
    }
  }
  
  CoherentSample* cs = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSignal);
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  
  int irr = 0;
  while(rmax <= RMAX){
    cs->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    cs->AddToHistVector(CoherentFitUtils::GetStoppingKaonsHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector()[irr], rmin, rmax, Chi2Cut));
    rmin += STEP;
    rmax += STEP;
    irr++;
  }
  
  if(normalize)cs->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
    
  return cs;
}

//********************************************************************
CoherentSample* CoherentFit::CreateTrueAllKaonsSample(const double RMIN, const double RMAX, const double STEP,
						      const double Chi2Cut,
						      const double bin_min, const double bin_max, const double bin_width,
						      const bool normalize, const bool resize) const {
//********************************************************************

  if(fIsMiniTree){
    fTree->GetEntry(0);
    if(!fSpill->GetIsMC()){
      std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
      std::exit(1);
    }
  }

  CoherentSample* cs = new CoherentSample(CoherentSample::SampleTypeEnum::kTrueSignal);
  
  double rmin = RMIN;
  double rmax = RMIN+STEP;

  int irr = 0;
  while(rmax <= RMAX){
    cs->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    cs->AddToHistVector(CoherentFitUtils::GetAllKaonsHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector()[irr], rmin, rmax, Chi2Cut));
    rmin += STEP;
    rmax += STEP;
    irr++;
  }
  
  if(normalize)cs->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
    
  return cs;
}

//********************************************************************
void CoherentFit::ComputeSelfSystematicError(){
//********************************************************************

  //InitializeHistogramsForSystematicErrors();
  std::vector<double> a,b,c,eff;
  
  TRandom3* r = new TRandom3();
  r->SetSeed(1);

  fToySample = fSignalPlusBackground->Clone();
  CoherentSample* ClonedSignal = fSignal->Clone();
  CoherentSample* ClonedBackground = fBackground->Clone();
  fToySample->SetSignal(ClonedSignal);
  fToySample->SetBackground(ClonedBackground);

  //generate toy samples
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::cout << "toy sample " << itoy << std::endl;
    //set initial parameters values with variations
    fToySample->GetSignal()->SetCFitParametersWithVariations(fSignal,r);
    fToySample->GetBackground()->SetCFitParametersWithVariations(fBackground,r);

    //coherent fit
    bool result = false;
    fToySample->CoherentFitSignalPlusBackgroundToy(result);
    //store results
    a.push_back(fToySample->GetSignal()->GetCmpvA().first);
    b.push_back(fToySample->GetSignal()->GetCmpvB().first);
    c.push_back(fToySample->GetSignal()->GetCmpvC().first);
    eff.push_back(result);
  }

  InitializeHistogramsForSystematicErrors(a,b,c,eff);
  FillSystematicHistograms(a,b,c,eff);
  
  TFile* rfile = new TFile("selfsyst.root","NEW");
  h_toy_A->Write();
  h_toy_B->Write();
  h_toy_C->Write();
  h_toy_Eff->Write();
  rfile->Close();
}

//********************************************************************
void CoherentFit::GenerateToySample(const bool apply_toy_weights, const bool apply_toy_variations, const int itoy){
//********************************************************************

  if(!fIsSystTree){
    std::cout << "ERROR! TTree is not a flat tree with systematic variations. Please check your input file" << std::endl;
    std::exit(1);
  }

  //delete previous one
  if(fToySample){
    delete fToySample;
    fToySample = NULL;
  }

  fToySample = new CoherentSample(CoherentSample::SampleTypeEnum::kSignalPlusBackground);
  GenerateToyHistograms(fToySample,apply_toy_weights,apply_toy_variations,itoy);
}

//********************************************************************
void CoherentFit::GenerateToyHistograms(CoherentSample* ToySample, const bool apply_toy_weights, const bool apply_toy_variations, const int itoy){
//********************************************************************

  if(!fIsSystTree){
    std::cout << "ERROR! TTree is not a flat tree with systematic errors. Please check your input file" << std::endl;
    std::exit(1);
  }

  if(!fSignalPlusBackground){
    std::cout << "ERROR! No original sample found. Please check your algorithm" << std::endl;
    std::exit(1);
  }
  
  double step = 2*fSignalPlusBackground->GetRRVector()[0].second;
  double rmin = fSignalPlusBackground->GetRRVector()[0].first - step/2;
  double rmax = fSignalPlusBackground->GetRRVector()[0].first + step/2;

  for(int irr = 0; irr < (int)fSignalPlusBackground->GetRRVector().size(); irr++){
    double bin_min   = fSignalPlusBackground->GetHistVector()[irr]->GetXaxis()->GetXmin();
    double bin_max   = fSignalPlusBackground->GetHistVector()[irr]->GetXaxis()->GetXmax();
    double bin_width = fSignalPlusBackground->GetHistVector()[irr]->GetXaxis()->GetBinWidth(0);
    ToySample->AddToHistVector(CoherentFitUtils::GetToyHistogramFromResRangeSlice(fTreeSystematics,
										  rmin, rmax,
										  bin_min, bin_max, bin_width,
										  apply_toy_weights, apply_toy_variations, itoy));
    rmin += step;
    rmax += step;
  }
  ToySample->SetRRVector(fSignalPlusBackground->GetRRVector());
  NormalizeHistograms(ToySample);
}

//********************************************************************
void CoherentFit::InitializeHistogramsForSystematicErrors(){
//********************************************************************

  double mpvA_mean = fSignalPlusBackground->GetSignal()->GetCmpvA().first;
  double mpvB_mean = fSignalPlusBackground->GetSignal()->GetCmpvB().first;
  double mpvC_mean = fSignalPlusBackground->GetSignal()->GetCmpvC().first;
  h_toy_A = new TH1F("h_toy_A","h_toy_A",500,0.1*mpvA_mean,1.9*mpvA_mean);
  h_toy_B = new TH1F("h_toy_B","h_toy_B",5000,0.1*mpvB_mean,1.9*mpvB_mean);
  h_toy_C = new TH1F("h_toy_C","h_toy_C",500,0.1*mpvC_mean,1.9*mpvC_mean);
  h_toy_Eff = new TH1F("h_toy_Eff","h_toy_Eff",2,0,2);
  h_toy_dEdx_RR = new TH2F("h_toy_dEdx_RR","h_toy_dEdx_RR",5900,1,60,2000,0,20);
}

//********************************************************************
void CoherentFit::InitializeHistogramsForSystematicErrors(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> eff){
//********************************************************************

  double amin = *std::min_element(a.begin(), a.end())*1.2; //a is negative
  double amax = *std::max_element(a.begin(), a.end())*0.8; //a is negative
  h_toy_A = new TH1F("h_toy_A","h_toy_A",5000,amin,amax); 

  double bmin = *std::min_element(b.begin(), b.end())*0.8;
  double bmax = *std::max_element(b.begin(), b.end())*1.2;
  h_toy_B = new TH1F("h_toy_B","h_toy_B",5000,bmin,bmax);

  double cmin = *std::min_element(c.begin(), c.end())*1.2; //c is negative
  double cmax = *std::max_element(c.begin(), c.end())*0.8; //c is negative
  h_toy_C = new TH1F("h_toy_C","h_toy_C",5000,cmin,cmax);

  h_toy_Eff = new TH1F("h_toy_Eff","h_toy_Eff",2,0,2);
}

//********************************************************************
void CoherentFit::FillSystematicHistograms(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> eff){
//********************************************************************

  for(int i = 0; i < (int)a.size(); i++){
    h_toy_Eff->Fill(eff[i]);
    if(eff[i] == 1){
      h_toy_A->Fill(a[i]);
      h_toy_B->Fill(b[i]);
      h_toy_C->Fill(c[i]);
    }
  }
}

//********************************************************************
void CoherentFit::NormalizeHistograms(){
//********************************************************************

  fSignalPlusBackground->NormalizeHistograms();
  if(GetIsMC()){
    fTrueSignal->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
    // fTrueSemiSignal->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
    fTrueBackground->NormalizeHistograms(fSignalPlusBackground->GetIntegralVector());
  }
}

//********************************************************************
void CoherentFit::NormalizeHistograms(CoherentSample* cs){
//********************************************************************

  cs->NormalizeHistograms();
}

//********************************************************************
void CoherentFit::FitToySample(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &eff){
//********************************************************************

  //create toy signal and background
  //this is done so the seed is never changed
  CoherentSample* toySignal = fSignal->Clone(); 
  CoherentSample* toyBackground = fBackground->Clone();
  CoherentSample* toySemiBackground = fBackground->GetSemiBackground()->Clone();
  toyBackground->SetSemiBackground(toySemiBackground);
  
  //make toy sample point to toysignal and toybackground
  fToySample->SetSignal(toySignal);
  fToySample->SetBackground(toyBackground);
  //initialize seed to original
  fToySample->GetSignal()->SetCFitParameters(fSignal);
  fToySample->GetBackground()->SetCFitParameters(fBackground);

  //coherent fit
  bool result = false;
  fToySample->CoherentFitSignalPlusBackgroundToy(result);

  //store results
  a.push_back(fToySample->GetSignal()->GetCmpvA().first);
  b.push_back(fToySample->GetSignal()->GetCmpvB().first);
  c.push_back(fToySample->GetSignal()->GetCmpvC().first);
  eff.push_back(result);
  
  std::cout << fSignalPlusBackground->GetSignal()->GetCmpvA().first << " " << fToySample->GetSignal()->GetCmpvA().first << std::endl;
  std::cout << fSignalPlusBackground->GetSignal()->GetCmpvB().first << " " << fToySample->GetSignal()->GetCmpvB().first << std::endl;
  std::cout << fSignalPlusBackground->GetSignal()->GetCmpvC().first << " " << fToySample->GetSignal()->GetCmpvC().first << std::endl;
}

//********************************************************************
void CoherentFit::SequentialCoherentFit(bool minos){
//********************************************************************
  
  fTrueSignal->SequentialCoherentFit();
  fTrueBackground->SequentialCoherentFit();
  fSignalPlusBackground->SequentialCoherentFit(minos);
}

//********************************************************************
void CoherentFit::DataCoherentFit(const CoherentFit* c){
//********************************************************************
  
  SetParametersFromMCFit(c);
  ScaleParameters();
  fSignalPlusBackground->CoherentFitSignalPlusBackground(true,1);
}

//********************************************************************
void CoherentFit::ScaleParameters(){
//********************************************************************

  double scale = 1.4;
  fSignal->SetClwA(std::make_pair(fSignal->GetClwA().first*scale,fSignal->GetClwA().second*scale));
  fSignal->SetCgwA(std::make_pair(fSignal->GetCgwA().first*scale,fSignal->GetCgwA().second*scale));
  fBackground->SetClwA(std::make_pair(fBackground->GetClwA().first*scale,fBackground->GetClwA().second*scale));
}

//********************************************************************
void CoherentFit::ToyLoop(const bool apply_toy_weights, const bool apply_toy_variations){
//********************************************************************

  //create vectors for storing results
  std::vector<double> a, b, c, eff;
  a.clear(); b.clear(); c.clear(); eff.clear();
  
  //loop over toys
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::cout << "---------------------" << std::endl;
    std::cout << "TOY " << itoy << std::endl;
    std::cout << "---------------------" << std::endl;
    GenerateToySample(apply_toy_weights, apply_toy_variations, itoy);
    FitToySample(a,b,c,eff);
  }

  InitializeHistogramsForSystematicErrors(a,b,c,eff);
  FillSystematicHistograms(a,b,c,eff);
}

//********************************************************************
void CoherentFit::WriteSystematicHistograms(const std::string& filename){
//********************************************************************

  TFile* rfile = new TFile(filename.c_str(),"NEW");
  h_toy_A->Write();
  h_toy_B->Write();
  h_toy_C->Write();
  h_toy_Eff->Write();
  // h_toy_dEdx_RR->Write();
  rfile->Close();
}

//********************************************************************
void CoherentFit::PropagateSystematicErrors(const std::string& filename, const bool apply_toy_weights, const bool apply_toy_variations){
//********************************************************************

  if(apply_toy_weights && !apply_toy_variations)std::cout << "propagating weight systematics" << std::endl;
  if(!apply_toy_weights && apply_toy_variations)std::cout << "propagating variation systematics" << std::endl;
  if(apply_toy_weights && apply_toy_variations) std::cout << "propagating all systematics" << std::endl;
  if(!fTreeSystematics){
    std::cout << "ERROR! no tree with systematics found" << std::endl;
    std::exit(1);
  }
  fTreeSystematics = GetSystTree();
  // InitializeHistogramsForSystematicErrors();
  ToyLoop(apply_toy_weights,apply_toy_variations);
  WriteSystematicHistograms(filename);
}


