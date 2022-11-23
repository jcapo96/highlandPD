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
  fIsMiniTree = false;
  fIsSystTree = false;
  fIsMC = false;
  fSpill = NULL;
  fSignalPlusBackground = NULL;
  fSignal = NULL;
  fBackground = NULL;
  fTrueSignal = NULL;
  fTrueBackground = NULL;
  fToySamples.clear();
  h_toy_A = NULL;
  h_toy_B = NULL;
  h_toy_C = NULL;
  h_toy_D = NULL;
  h_toy_R = NULL;
  h_toy_dEdx_RR = NULL;
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
  fBackground = NULL;
  fTrueSignal = NULL;
  fTrueBackground = NULL;
  fToySamples.clear();

  h_toy_A = NULL;
  h_toy_B = NULL;
  h_toy_C = NULL;
  h_toy_D = NULL;
  h_toy_R = NULL;
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

  t = (TTree*)fFile->Get("all_syst");
  //t = (TTree*)fFile->Get("dQdx_YZcal");
  if(t){
    fIsMiniTree = false;
    fIsSystTree = true;
    fIsMC = true;
    std::cout << "Input tree is flat tree with systematics" << std::endl;
    return t;
  }

  //if here, no known type of tree provided
  std::cout << "ERROR! TTree type not known. Please check your input file" << std::endl;
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

  if(fIsMiniTree){
    fTree->GetEntry(0);
    if(!fSpill->GetIsMC()){
      std::cout << "this tree is not MC!!! cannot create signal and background histograms" << std::endl;
      std::exit(1);
    }
  }    
  double rmin = RMIN;
  double rmax = RMIN+STEP;
  
  while(rmax <= RMAX){
    fSignalPlusBackground->AddToRRVector(std::make_pair((rmin+rmax)/2,(rmax-rmin)/2));
    fSignalPlusBackground->AddToHistVector(CoherentFitUtils::GetHistogramFromResRangeSlice(fTree, fSpill, rmin, rmax, Chi2Cut, bin_min, bin_max, bin_width));
    fTrueSignal->AddToHistVector(CoherentFitUtils::GetSignalHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
    fTrueBackground->AddToHistVector(CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(fTree, fSpill, fSignalPlusBackground->GetHistVector().back(), rmin, rmax, Chi2Cut));
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
void CoherentFit::ComputeSelfSystematicError(const bool apply_all_var,
					     const bool apply_only_s_lw_var,
					     const bool apply_only_mpv_var,
					     const bool apply_only_s_gw_var,
					     const bool apply_only_b_lw_var,
					     const bool apply_only_b_gw_var,
					     const bool apply_only_shift_var,
					     const bool apply_only_norm_var){
//********************************************************************

  TRandom3* r = new TRandom3();
  r->SetSeed(1);

  CoherentSample* ClonedSample = fSignalPlusBackground->Clone();
  ClonedSample->SetSignal(fSignal->Clone());
  ClonedSample->SetBackground(fBackground->Clone());

  //generate toy samples
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::cout << "toy sample " << itoy << std::endl;
    fToySamples.push_back(ClonedSample);
    //generate seed values for fitting
    fToySamples[itoy]->GetSignal()->SetCFitParametersWithVariations(fSignal,r,0.05,
								    apply_all_var,
								    apply_only_s_lw_var,apply_only_mpv_var,apply_only_s_gw_var,
								    apply_only_shift_var,apply_only_norm_var);
    fToySamples[itoy]->GetBackground()->SetCFitParametersWithVariations(fBackground,r,0.05,
									apply_all_var,
									apply_only_b_lw_var,apply_only_mpv_var,apply_only_b_gw_var,
									apply_only_shift_var,apply_only_norm_var);
    //fit
    fToySamples[itoy]->CoherentFit();
    //store results
    for(int ibin = 0; ibin < 5900; ibin++)
      h_toy_dEdx_RR->Fill(h_toy_dEdx_RR->GetXaxis()->GetBinCenter(ibin+1),
			  fToySamples[itoy]->GetSignal()->GetCmpvFit()->Eval(h_toy_dEdx_RR->GetXaxis()->GetBinCenter(ibin+1)));
    
    /*for(int ihist = 0; ihist < fToySamples[itoy]->GetHistVector().size(); ihist++){
      fToySamples[itoy]->GetHistVector()[ihist]->Draw("histo");
      fToySamples[itoy]->GetCFitVector()[ihist]->Draw("same");
      gPad->Update();gSystem->Sleep(10000);
      }*/
  }
  
  TFile* rfile = new TFile("selfsyst_norm.root","NEW");
  h_toy_dEdx_RR->Write();
  rfile->Close();
}

//********************************************************************
void CoherentFit::GenerateToySamples(const bool apply_toy_weights, const bool apply_toy_variations){
//********************************************************************

  if(!fIsSystTree){
    std::cout << "ERROR! TTree is not a flat tree with systematic variations. Please check your input file" << std::endl;
    std::exit(1);
  }
  
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::cout << "toy sample " << itoy << std::endl;
    fToySamples.push_back(new CoherentSample(CoherentSample::SampleTypeEnum::kSignalPlusBackground));
    GenerateToyHistograms(fToySamples.back(),apply_toy_weights,apply_toy_variations,itoy);
  }

  /*for(int i = 0; i < fSignalPlusBackground->GetRRVector().size(); i++){
    for(int j = 0; j < fToySamples.size(); j++){
      fToySamples[j]->SetHistogramColor(j+1);
      fToySamples[j]->SetHistogramLineWidth(2);
      std::stringstream ssj;
      ssj << j+1;
      fToySamples[j]->GetHistVector()[i]->SetTitle(("Toy "+ssj.str()+"").c_str());
      if(j==0){
	fToySamples[j]->GetHistVector()[i]->GetXaxis()->SetTitle("dEdx [MeV/cm]");
	fToySamples[j]->GetHistVector()[i]->Draw("histo");
      }
      else    fToySamples[j]->GetHistVector()[i]->Draw("histosame");
    }
    gPad->BuildLegend();
    gPad->Update();gSystem->Sleep(20000);//gPad->WaitPrimitive();
    }*/
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
    ToySample->AddToHistVector(CoherentFitUtils::GetToyHistogramFromResRangeSlice(fTree,
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
  double mpvD_mean = fSignalPlusBackground->GetSignal()->GetCmpvD().first;
  double mpvR_mean = fSignalPlusBackground->GetSignal()->GetCmpvR().first;
  h_toy_A = new TH1F("h_toy_A","h_toy_A",500,0.1*mpvA_mean,1.9*mpvA_mean);
  h_toy_B = new TH1F("h_toy_B","h_toy_B",500,0.1*mpvB_mean,1.9*mpvB_mean);
  h_toy_C = new TH1F("h_toy_C","h_toy_C",500,0.1*mpvC_mean,1.9*mpvC_mean);
  h_toy_D = new TH1F("h_toy_D","h_toy_D",500,0.1*mpvD_mean,1.9*mpvD_mean);
  h_toy_R = new TH1F("h_toy_R","h_toy_R",500,0.1*mpvR_mean,1.9*mpvR_mean);
  h_toy_dEdx_RR = new TH2F("h_toy_dEdx_RR","h_toy_dEdx_RR",5900,1,60,2000,0,20);
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
void CoherentFit::NormalizeHistograms(CoherentSample* cs){
//********************************************************************

  cs->NormalizeHistograms();
}

//********************************************************************
void CoherentFit::FitToySamples(){
//********************************************************************

  //create toy signal and background
  CoherentSample* toySignal = fSignal->Clone();
  CoherentSample* toyBackground = fBackground->Clone();
  
  for(int itoy = 0; itoy < NTOYS; itoy++){
    std::cout << "fitting toy sample " << itoy << std::endl;
    //make toy sample point to toysignal and toybackground
    fToySamples[itoy]->SetSignal(toySignal);
    fToySamples[itoy]->SetBackground(toyBackground);
    fToySamples[itoy]->SetBackgroundModel(fSignalPlusBackground->GetBackgroundModel());
    //initialize seed to original
    fToySamples[itoy]->GetSignal()->SetCFitParameters(fSignal);
    fToySamples[itoy]->GetBackground()->SetCFitParameters(fBackground);
    //coherent fit
    fToySamples[itoy]->CoherentFit();
    //store results
    h_toy_A->Fill(fToySamples[itoy]->GetSignal()->GetCmpvA().first);
    h_toy_B->Fill(fToySamples[itoy]->GetSignal()->GetCmpvB().first);
    h_toy_C->Fill(fToySamples[itoy]->GetSignal()->GetCmpvC().first);
    h_toy_D->Fill(fToySamples[itoy]->GetSignal()->GetCmpvD().first);
    h_toy_R->Fill(fToySamples[itoy]->GetSignal()->GetCmpvR().first);
    for(int ibin = 0; ibin < 5900; ibin++)
      h_toy_dEdx_RR->Fill(h_toy_dEdx_RR->GetXaxis()->GetBinCenter(ibin+1),
			  fToySamples[itoy]->GetSignal()->GetCmpvFit()->Eval(h_toy_dEdx_RR->GetXaxis()->GetBinCenter(ibin+1)));
  }

  /*for(int i = 0; i < fSignalPlusBackground->GetRRVector().size(); i++){
    for(int j = 0; j < fToySamples.size(); j++){
      fToySamples[j]->SetHistogramColor(j+1);
      fToySamples[j]->SetCFitColor(j+1);
      fToySamples[j]->SetHistogramLineWidth(2);
      std::stringstream ssj;
      ssj << j+1;
      fToySamples[j]->GetHistVector()[i]->SetTitle(("Toy "+ssj.str()+"").c_str());
      fToySamples[j]->GetCFitVector()[i]->SetTitle(("Toy "+ssj.str()+" Fit").c_str());
      if(j==0){
	fToySamples[j]->GetHistVector()[i]->GetXaxis()->SetTitle("dEdx [MeV/cm]");
	fToySamples[j]->GetHistVector()[i]->Draw("histo");
	fToySamples[j]->GetCFitVector()[i]->Draw("same");
      }
      else{
	fToySamples[j]->GetHistVector()[i]->Draw("histosame");
	fToySamples[j]->GetCFitVector()[i]->Draw("same");
      }
    }
    gPad->BuildLegend();
    gPad->Update();gPad->WaitPrimitive();
  }*/
  
  TFile* rfile = new TFile("syst_histos_Recombination.root","NEW");
  h_toy_A->Write();
  h_toy_B->Write();
  h_toy_C->Write();
  h_toy_D->Write();
  h_toy_R->Write();
  h_toy_dEdx_RR->Write();
  rfile->Close();
  
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

//********************************************************************
void CoherentFit::PropagateSystematicErrors(const bool apply_toy_weights, const bool apply_toy_variations){
//********************************************************************

  if(apply_toy_weights && !apply_toy_variations)std::cout << "propagating weight systematics" << std::endl;
  if(!apply_toy_weights && apply_toy_variations)std::cout << "propagating variation systematics" << std::endl;
  if(apply_toy_weights && apply_toy_variations) std::cout << "propagating all systematics" << std::endl;
  GenerateToySamples(apply_toy_weights,apply_toy_variations);
  InitializeHistogramsForSystematicErrors();
  FitToySamples();
}


