#include "CoherentFitUtils.hxx"

#include <iostream>

#include "TPad.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"

//********************************************
TH1F* CoherentFitUtils::CreateHistogram(double bin_min, double bin_max, double bin_width){
//********************************************

  if(bin_min>=bin_max || bin_width<=0){
    std::cout << "Invalid histogram binning" << std::endl;
    std::exit(1);
  }

  int nbins = ceil((bin_max-bin_min)/bin_width);
  bin_max = bin_min + nbins*bin_width;
  return new TH1F("h","h",nbins,bin_min,bin_max);

}

//********************************************
TH1F* CoherentFitUtils::GetHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill,
						      const double RMIN, const double RMAX,
						      const double Chi2Cut,
						      double bin_min, double bin_max, double bin_width){
//********************************************

  if(spill){
    TH1F* h_dummy = CreateHistogram(bin_min,bin_max,bin_width);
    std::vector<double> dEdx;
    dEdx.clear();
    
    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	    if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r){
	      h_dummy->Fill(part->Hits[2][ihit].dEdx);
	      dEdx.push_back(part->Hits[2][ihit].dEdx);
	    }
	  }
	}
      }
    }

    double new_bin_min = 0;
    for(int ibin = 1; ibin < h_dummy->GetNbinsX()-1; ibin++){
      if(h_dummy->GetBinContent(ibin+1)>5 && h_dummy->GetBinContent(ibin+2)>5){
	new_bin_min = h_dummy->GetBinLowEdge(ibin);
	break;
      }
    }
    new_bin_min = std::max(bin_min,new_bin_min);
    
    double new_bin_max = 0;
    for(int ibin = h_dummy->GetMaximumBin(); ibin < h_dummy->GetNbinsX()-1; ibin++){
      if(h_dummy->GetBinContent(ibin+1)<2 && h_dummy->GetBinContent(ibin+2)<2){
	new_bin_max = h_dummy->GetBinLowEdge(ibin);
	break;
      }
    }
    new_bin_max = std::min(new_bin_max,bin_max);
    
    TH1F* h = CreateHistogram(new_bin_min,new_bin_max,bin_width);
    for(int i = 0; i < (int)dEdx.size(); i++)h->Fill(dEdx[i]);
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("h_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;//ChangeHistogramToVariableBinning(h,1,10);
  }
  else{
    //fill histogram
    TH1F* h_dummy = GetHistogramFromResRangeSliceFromFlatTree(t,RMIN,RMAX,bin_min,bin_max,bin_width);

    double new_bin_min = 0;
    for(int ibin = 1; ibin < h_dummy->GetNbinsX()-1; ibin++){
      if(h_dummy->GetBinContent(ibin+1)>5 && h_dummy->GetBinContent(ibin+2)>5){
	new_bin_min = h_dummy->GetBinLowEdge(ibin);
	break;
      }
    }
    new_bin_min = std::max(bin_min,new_bin_min);

    double new_bin_max = 0;
    for(int ibin = h_dummy->GetMaximumBin(); ibin < h_dummy->GetNbinsX()-1; ibin++){
      new_bin_max = h_dummy->GetBinLowEdge(ibin);
      if(h_dummy->GetBinContent(ibin+1)<2 && h_dummy->GetBinContent(ibin+2)<2){
	new_bin_max = h_dummy->GetBinLowEdge(ibin);
	break;
      }
    }
    new_bin_max = std::min(new_bin_max,bin_max);

    int new_nbins = (new_bin_max-new_bin_min)/bin_width;
    std::vector<double> bins;
    bins.clear();
    for(int i = 0; i < new_nbins+1; i++)bins.push_back(new_bin_min+i*bin_width);
    TH1F* h = (TH1F*)h_dummy->Rebin(new_nbins,"rebin",&bins[0]);
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("h_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GetSignalHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
							    const double RMIN, const double RMAX,
							    const double Chi2Cut){
//********************************************

  if(spill){
    if(!spill->GetIsMC()){
      std::cout << "this MiniTree is not MC, can't get signal histogram!" << std::endl;
      std::exit(1);
    }
    
    TH1F* h = (TH1F*)ha->Clone();
    h->Reset();

    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2 || abs(part->PositionEnd[2]-230)<20)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	  if(!truePart)continue;
	  //if(truePart->PDG==321 && truePart->ProcessEnd == 2 && part->Chi2Proton/part->Chi2ndf<100){
	  //if(truePart->PDG==321 && truePart->ProcessEnd == 2 && part->Chi2Proton/part->Chi2ndf<90){
	  // if(truePart->PDG==321 && truePart->ProcessEnd == 2 && part->Chi2Proton/part->Chi2ndf<80){
	  //if(truePart->PDG==321 && truePart->ProcessEnd == 2){
	  if(truePart->PDG==321){
	    for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	      if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r)
		h->Fill(part->Hits[2][ihit].dEdx);
	    }
	  }
	}
      }
    }

    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
  else{
    //fill histogram
    TH1F* h = GetHistogramFromResRangeSliceFromFlatTree(t,
							RMIN,RMAX,
							ha->GetXaxis()->GetXmin(), ha->GetXaxis()->GetXmax(), ha->GetXaxis()->GetBinWidth(0),
							"&& (bestcandidate_truepdg==321)"); 
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GetBackgroundHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
								const double RMIN, const double RMAX,
								const double Chi2Cut){
//********************************************

  if(spill){
    if(!spill->GetIsMC()){
      std::cout << "this MiniTree is not MC, can't get background histogram!" << std::endl;
      std::exit(1);
    }
    
    TH1F* h = (TH1F*)ha->Clone();
    h->Reset();
    
    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	  if(!truePart)continue;
	  if(truePart->PDG==321)continue;
	  //if(truePart->PDG==321 && truePart->ProcessEnd == 2)continue;
	  //if(truePart->PDG==321 && truePart->ProcessEnd == 2 && part->Chi2Proton/part->Chi2ndf<90)continue;
	  //if(truePart->PDG==321 && truePart->ProcessEnd==2 && part->Chi2Proton/part->Chi2ndf<90)continue;
	  //if(truePart->PDG==321 && truePart->ProcessEnd != 2){
	  //if(truePart->PDG==321 && truePart->ProcessEnd != 2){
	  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	    if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r)
	      h->Fill(part->Hits[2][ihit].dEdx);
	  }
	}
      }
    }
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hb_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
  else{
    //fill histogram
    TH1F* h = GetHistogramFromResRangeSliceFromFlatTree(t,
							RMIN,RMAX,
							ha->GetXaxis()->GetXmin(), ha->GetXaxis()->GetXmax(), ha->GetXaxis()->GetBinWidth(0),
							"&& !(bestcandidate_truepdg==321)");
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hb_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GetSemiSignalHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
								const double RMIN, const double RMAX,
								const double Chi2Cut){
//********************************************

  if(spill){
    if(!spill->GetIsMC()){
      std::cout << "this MiniTree is not MC, can't get background histogram!" << std::endl;
      std::exit(1);
    }
    
    TH1F* h = (TH1F*)ha->Clone();
    h->Reset();
    
    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	  if(!truePart)continue;
	  if((truePart->PDG==321 && truePart->ProcessEnd != 2 && part->Chi2Proton/part->Chi2ndf<100) ||
	     (truePart->PDG==321 && truePart->ProcessEnd == 2 && part->Chi2Proton/part->Chi2ndf>100)){
	    for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	      if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r)
		h->Fill(part->Hits[2][ihit].dEdx);
	    }
	  }
	}
      }
    }
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hb_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
  else{
    //fill histogram
    TH1F* h = GetHistogramFromResRangeSliceFromFlatTree(t,
							RMIN,RMAX,
							ha->GetXaxis()->GetXmin(), ha->GetXaxis()->GetXmax(), ha->GetXaxis()->GetBinWidth(0),
							"&& !(bestcandidate_truepdg==321)"); 
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hss_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GetStoppingKaonsHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
								   const double RMIN, const double RMAX,
								   const double Chi2Cut){
//********************************************

  if(spill){
    if(!spill->GetIsMC()){
      std::cout << "this MiniTree is not MC, can't get signal histogram!" << std::endl;
      std::exit(1);
    }
    
    TH1F* h = (TH1F*)ha->Clone();
    h->Reset();
    
    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	  if(!truePart)continue;
	  if(truePart->PDG==321 && truePart->ProcessEnd == 2){
	    for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	      if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r)
		h->Fill(part->Hits[2][ihit].dEdx);
	    }
	  }
	}
      }
    }
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
  else{
    //fill histogram
    std::stringstream sschi;
    sschi << Chi2Cut;
    TH1F* h = GetHistogramFromResRangeSliceFromFlatTree(t,
							RMIN,RMAX,
							ha->GetXaxis()->GetXmin(), ha->GetXaxis()->GetXmax(), ha->GetXaxis()->GetBinWidth(0),
							("&& (bestcandidate_truepdg==321) && (bestcandidate_trueendproc==2) && (bestcandidate_chi2_prot/bestcandidate_chi2_ndf<"+sschi.str()+")").c_str()); 
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GetAllKaonsHistogramFromResRangeSlice(TTree* t, AnaSpillB* spill, TH1F* ha,
							      const double RMIN, const double RMAX,
							      const double Chi2Cut){
//********************************************

  if(spill){
    if(!spill->GetIsMC()){
      std::cout << "this MiniTree is not MC, can't get signal histogram!" << std::endl;
      std::exit(1);
    }
    
    TH1F* h = (TH1F*)ha->Clone();
    h->Reset();
    
    //fill histogram
    double rr_0 = (RMIN+RMAX)/2;
    double rr_r = (RMAX-RMIN)/2;
    
    //loop over candidates
    for(int ientry = 0; ientry < t->GetEntries(); ientry++){
      t->GetEntry(ientry);
      AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
      for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
	AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
	if(!part || (int)part->Hits[2].size()<2)continue;
	if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	  if(!truePart)continue;
	  if(truePart->PDG==321){
	    for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	      if(abs(part->Hits[2][ihit].ResidualRange-rr_0)<rr_r)
		h->Fill(part->Hits[2][ihit].dEdx);
	    }
	  }
	}
      }
    }
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
  else{
    //fill histogram
    std::stringstream sschi;
    sschi << Chi2Cut;
    TH1F* h = GetHistogramFromResRangeSliceFromFlatTree(t,
							RMIN,RMAX,
							ha->GetXaxis()->GetXmin(), ha->GetXaxis()->GetXmax(), ha->GetXaxis()->GetBinWidth(0),
							("&& (bestcandidate_truepdg==321) && (bestcandidate_chi2_prot/bestcandidate_chi2_ndf<"+sschi.str()+")").c_str()); 
    
    //set histogram title
    std::stringstream srl, srh;
    srl << RMIN;
    srh << RMAX;
    std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
    h->SetName(("hs_"+srl.str()+"_"+srh.str()+"").c_str());
    h->SetTitle((st).c_str());
    
    return h;
  }
}

//********************************************
TH1F* CoherentFitUtils::GenerateBackgroundHistogramFromTrueSignal(TTree* t, AnaSpillB* spill, TH1F* ha,
								  const double RMIN, const double RMAX,
								  const double Chi2Cut,
								  const std::vector<double> shift){
//********************************************
  
  if(!spill->GetIsMC()){
    std::cout << "this MiniTree is not MC, can't get background histogram!" << std::endl;
    std::exit(1);
  }
  
  TH1F* h = (TH1F*)ha->Clone();
  h->Reset();
  
  //fill histogram
  double rr_0 = (RMIN+RMAX)/2;
  double rr_r = (RMAX-RMIN)/2;
  
  //loop over candidates
  for(int ientry = 0; ientry < t->GetEntries(); ientry++){
    t->GetEntry(ientry);
    AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
    for(int ipart = 0; ipart < (int)bunch->Particles.size(); ipart++){
      if(shift[ientry+ipart]<=0)continue;
      AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
      if(!part || (int)part->Hits[2].size()<2)continue;
      if(part->Chi2Proton/part->Chi2ndf<Chi2Cut && part->Chi2Proton>0){
	AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
	if(!truePart)continue;
	if(truePart->PDG==321 && truePart->ProcessEnd == 2){
	  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){
	    if(abs((part->Hits[2][ihit].ResidualRange-shift[ientry+ipart])-rr_0)<rr_r)
	      h->Fill(part->Hits[2][ihit].dEdx);
	  }
	}
      }
    }
  }
  
  //set histogram title
  std::stringstream srl, srh;
  srl << RMIN;
  srh << RMAX;
  std::string st = srl.str()+" < Residual Range [cm] < "+srh.str();
  h->SetName(("h_fakeb_"+srl.str()+"_"+srh.str()+"").c_str());
  h->SetTitle((st).c_str());

  return h;
}

//********************************************
TH1F* CoherentFitUtils::GetHistogramFromResRangeSliceFromFlatTree(TTree* t,
								  const double RMIN, const double RMAX,
								  double bin_min, double bin_max, double bin_width,
								  std::string cut,
								  const bool apply_toy_weights, const bool apply_toy_variations, int itoy){
//********************************************

  std::stringstream ssr0, ssrr, ssitoy;
  ssr0 << (RMIN+RMAX)/2;
  ssrr << (RMAX-RMIN)/2;
  ssitoy << itoy;
  
  // TH1F* h = CreateHistogram(bin_min,bin_max,bin_width);
  // TH1F* h_dummy = CreateHistogram(bin_min,bin_max,bin_width);
  // h_dummy->SetTitle("h_dummy");
  // h_dummy->SetName("h_dummy");

  // for(int i = 0; i < NHITS; i++){
  //   std::stringstream ssi;
  //   ssi << i;

  //   if(!apply_toy_weights && !apply_toy_variations){
  //     t->Project("h_dummy",
  // 		 ("bestcandidate_hit_dedx["+ssi.str()+"]").c_str(),
  // 		 ("accum_level[0][]>6 && abs(bestcandidate_hit_resrange["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  //   }
  //   else if(!apply_toy_weights && apply_toy_variations){
  //     t->Project("h_dummy",
  // 		 ("bestcandidate_hit_dedx_toy["+ssitoy.str()+"]["+ssi.str()+"]").c_str(),
  // 		 //("accum_level["+ssitoy.str()+"][0]>6 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"]["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  // 		 ("accum_level["+ssitoy.str()+"][0]>6 && abs(bestcandidate_hit_resrange["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  //   } 
  //   else if(apply_toy_weights && !apply_toy_variations){
  //     t->Project("h_dummy",
  // 		 ("bestcandidate_hit_dedx["+ssi.str()+"]").c_str(),
  // 		 ("(accum_level["+ssitoy.str()+"][0]>6 && abs(bestcandidate_hit_resrange["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+")*weight_syst_total["+ssitoy.str()+"]").c_str(),"");
  //   }
  //   else{
  //     t->Project("h_dummy",
  // 		 ("bestcandidate_hit_dedx["+ssi.str()+"]").c_str(),
  // 		 ("(accum_level["+ssitoy.str()+"][0]>6 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"]["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+")*weight_syst_total["+ssitoy.str()+"]").c_str(),"");
  //   }
  //   h->Add(h_dummy);
  //   h_dummy->Reset();
  // }
  // delete h_dummy;
  // return h;
  
  TH1F* h = CreateHistogram(bin_min,bin_max,bin_width);
  TH1F* h_dummy = CreateHistogram(bin_min,bin_max,bin_width);
  h_dummy->SetTitle("h_dummy");
  h_dummy->SetName("h_dummy");

  if(!apply_toy_weights && !apply_toy_variations){
    t->Project("h_dummy",
  	       "bestcandidate_hit_dedx[]",
  	       ("accum_level[0]>2 && abs(bestcandidate_hit_resrange[]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  }
  else if(!apply_toy_weights && apply_toy_variations){
    t->Project("h_dummy",
  	       ("bestcandidate_hit_dedx_toy["+ssitoy.str()+"][]").c_str(),
  	       //("accum_level["+ssitoy.str()+"][0]>8 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"]["+ssi.str()+"]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  	       ("accum_level["+ssitoy.str()+"]>2 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"][]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+"").c_str(),"");
  } 
  else if(apply_toy_weights && !apply_toy_variations){
    t->Project("h_dummy",
  	       ("bestcandidate_hit_dedx_toy["+ssitoy.str()+"][]").c_str(),
  	       ("(accum_level["+ssitoy.str()+"]>2 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"][]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+")*weight_syst_total["+ssitoy.str()+"]").c_str(),"");
  }
  else{
    t->Project("h_dummy",
  	       ("bestcandidate_hit_dedx_toy["+ssitoy.str()+"][]").c_str(),
  	       ("(accum_level["+ssitoy.str()+"]>2 && abs(bestcandidate_hit_resrange_toy["+ssitoy.str()+"][]-"+ssr0.str()+")<"+ssrr.str()+" "+cut+")*weight_syst_total["+ssitoy.str()+"]").c_str(),"");
  }

  h->Add(h_dummy);
  delete h_dummy;
  return h;
}

//********************************************
TH1F* CoherentFitUtils::GetToyHistogramFromResRangeSlice(TTree* t,
							 const double RMIN, const double RMAX,
							 double bin_min, double bin_max, double bin_width,
							 const bool apply_toy_weights, const bool apply_toy_variations,
							 const int itoy){
//********************************************

  return GetHistogramFromResRangeSliceFromFlatTree(t,
						   RMIN,RMAX,
						   bin_min,bin_max,bin_width,
						   "",
						   apply_toy_weights,apply_toy_variations,itoy);
}

//********************************************
TH1F* CoherentFitUtils::ChangeHistogramToVariableBinning(TH1F* h_original, const double I, const int min){
//********************************************

  if(!h_original){
    std::cout << "invalid pointer to histogram" << std::endl;
    std::exit(1);
  }

  if(min <= 1){
    std::cout << "minimum value per bin should be larger than one" << std::endl;
    std::exit(1);
  }

  int counts = 0;
  std::vector<double> edges, contents;
  edges.clear(); contents.clear();
  edges.push_back(h_original->GetBinLowEdge(1));
  for(int ibin = 0; ibin < h_original->GetNbinsX(); ibin++){
    counts += h_original->GetBinContent(ibin+1)*I;
    if(counts>=min){
      edges.push_back(h_original->GetBinLowEdge(ibin+1)+h_original->GetBinWidth(ibin+1));
      contents.push_back(counts);
      counts = 0;
    }
  }
  if(counts != 0){
    edges.pop_back();
    edges.push_back(h_original->GetBinLowEdge(h_original->GetNbinsX()+1)+h_original->GetBinWidth(h_original->GetNbinsX()+1));
    contents.back() += counts;
  }

  TH1F* h_new = new TH1F(h_original->GetName(), h_original->GetTitle(), edges.size()-1, &edges[0]);
  //for(int ibin = 0; ibin < h_original->GetNbinsX(); ibin++)
  //  h_new->Fill(h_original->GetBinCenter(ibin+1),h_original->GetBinContent(ibin+1));
  for(int ibin = 0; ibin < (int)contents.size(); ibin++){
    h_new->SetBinContent(ibin+1,contents[ibin]/I);
    h_new->SetBinError(ibin+1,sqrt(contents[ibin])/I);
  }

  //h_original->Draw();gPad->Update();gPad->WaitPrimitive();
  //h_new->Draw("e");gPad->Update();gPad->WaitPrimitive();
  
  delete h_original;
  return h_new;
}

//********************************************
std::vector<TH1F*> CoherentFitUtils::GetReducedHistograms(std::vector<TH1F*> h_originals, const double xmin, const double xmax, const int min_counts){
//********************************************

  std::vector<TH1F*> h;
  
  double bin_min = 0;
  double bin_max = 0;
  double dis_to_min = 100;
  double dis_to_max = 100;
  for(int ihist = 0; ihist < h_originals.size(); ihist++){
    for(int ibin = 0; ibin < h_originals[ihist]->GetNbinsX(); ibin++){
      if(abs(h_originals[ihist]->GetBinLowEdge(ibin+1)-xmin)<dis_to_min){
	dis_to_min = abs(h_originals[ihist]->GetBinLowEdge(ibin+1)-xmin);
	bin_min = h_originals[ihist]->GetBinLowEdge(ibin+1);
      }
      if(abs(h_originals[ihist]->GetBinLowEdge(ibin+1)+h_originals[ihist]->GetBinWidth(ibin+1)-xmax)<dis_to_max){
	dis_to_max = abs(h_originals[ihist]->GetBinLowEdge(ibin+1)+h_originals[ihist]->GetBinWidth(ibin+1)-xmax);
	bin_max = h_originals[ihist]->GetBinLowEdge(ibin+1)+h_originals[ihist]->GetBinWidth(ibin+1);
      }
    }
    std::cout << bin_min << " " << bin_max << std::endl;
    h.push_back(CreateHistogram(bin_min,bin_max,0.1));
    for(int ibin = 0; ibin < h[ihist]->GetNbinsX(); ibin++){
      for(int jbin = 0; jbin < h_originals[ihist]->GetNbinsX(); jbin++){
	if(h[ihist]->GetBinCenter(ibin+1)==h_originals[ihist]->GetBinCenter(jbin+1))
	  h[ihist]->SetBinContent(ibin+1,h_originals[ihist]->GetBinContent(jbin+1));
      }
    }
  }

  return h;
}

//********************************************
void CoherentFitUtils::CopyHistogramBinning(TH1F* h_original, TH1F* h_copy){
//********************************************

  TH1F* h_new = (TH1F*)h_copy->Clone();
  h_new->Reset();

  for(int in = 0; in < h_new->GetNbinsX(); in++)
    for(int io = 0; io < h_original->GetNbinsX(); io++)
      if(h_new->GetBinCenter(in+1) == h_original->GetBinCenter(io+1))
	h_new->SetBinContent(in+1,h_original->GetBinContent(io+1));

  delete h_original;
  h_original = h_new;
}

//********************************************
Double_t CoherentFitUtils::GetFunctionNormalizationInsideHistogramBoundaries(const TH1F* h, const TF1* f){
//********************************************

  if(!h || !f){
    std::cout << "no histo or function provided, returning -1" << std::endl;
    return -1.;
  }

  double norm = 0;
  for(int i = 0; i < h->GetNbinsX(); i++)
    norm += f->Eval(h->GetBinCenter(i+1))*h->GetBinWidth(i+1);

  return norm;
}

//********************************************
Double_t CoherentFitUtils::Langaus(Double_t *x, Double_t *par) {
//********************************************
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  //Double_t np = 100.0;      // number of convolution steps
  Double_t np = 1000.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

//********************************************
Double_t CoherentFitUtils::LangausPlusConstant(Double_t *x, Double_t *par) {
//********************************************
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  //Double_t np = 100.0;      // number of convolution steps
  Double_t np = 1000.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3])+par[4];
}

//********************************************
Double_t CoherentFitUtils::DoubleLangaus(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double blw    = par[4];
  double bmpv   = par[5];
  double bnorm  = par[6];
  double bgw    = par[7];

  double spar[] = {slw,smpv,snorm,sgw};
  double bpar[] = {blw,bmpv,bnorm,bgw};

  float S = CoherentFitUtils::Langaus(x,spar);
  float B = CoherentFitUtils::Langaus(x,bpar);

  return S+B;
  
}

//********************************************
Double_t CoherentFitUtils::TripleLangaus(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double sslw   = par[4];
  double ssmpv  = par[5];
  double ssnorm = par[6];
  double ssgw   = par[7];
  double blw    = par[8];
  double bmpv   = par[9];
  double bnorm  = par[10];
  double bgw    = par[11];

  double spar[]  = {slw ,smpv ,snorm ,sgw};
  double sspar[] = {sslw,ssmpv,ssnorm,ssgw};
  double bpar[]  = {blw ,bmpv ,bnorm ,bgw};

  float S  = CoherentFitUtils::Langaus(x,spar);
  float SS = CoherentFitUtils::Langaus(x,sspar);
  float B  = CoherentFitUtils::Langaus(x,bpar);

  return S+SS+B; 
}

//********************************************
Double_t CoherentFitUtils::DoubleLangausPlusConstant(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double blw    = par[4];
  double bmpv   = par[5];
  double bnorm  = par[6];
  double bgw    = par[7];
  double C      = par[8];

  double spar[] = {slw,smpv,snorm,sgw};
  double bpar[] = {blw,bmpv,bnorm,bgw};

  float S = CoherentFitUtils::Langaus(x,spar);
  float B = CoherentFitUtils::Langaus(x,bpar);

  return S+B+C;
  
}

//********************************************
Double_t CoherentFitUtils::TripleLangausPlusConstant(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double sslw   = par[4];
  double ssmpv  = par[5];
  double ssnorm = par[6];
  double ssgw   = par[7];
  double blw    = par[8];
  double bmpv   = par[9];
  double bnorm  = par[10];
  double bgw    = par[11];
  double C      = par[12];

  double spar[]  = {slw ,smpv ,snorm ,sgw};
  double sspar[] = {sslw,ssmpv,ssnorm,ssgw};
  double bpar[]  = {blw ,bmpv ,bnorm ,bgw};

  float S  = CoherentFitUtils::Langaus(x,spar);
  float SS = CoherentFitUtils::Langaus(x,sspar);
  float B  = CoherentFitUtils::Langaus(x,bpar);

  return S+SS+B+C; 
}

//********************************************
Double_t CoherentFitUtils::LangausPlusLandau(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double bnorm  = par[4];
  double bmpv   = par[5];
  double blw    = par[6];

  double lpar[]  = {slw ,smpv ,snorm ,sgw};

  float L  = CoherentFitUtils::Langaus(x,lpar);
  float G = bnorm*TMath::Gaus(*x,bmpv,blw,true);

  return L+G;
  
}

//********************************************
Double_t CoherentFitUtils::LangausPlusGaus(Double_t *x, Double_t *par) {
//********************************************

  double llw    = par[0];
  double lmpv   = par[1];
  double lnorm  = par[2];
  double lgw    = par[3];
  double gnorm  = par[4];
  double gmu    = par[5];
  double gsigma = par[6];

  double lpar[]  = {llw ,lmpv ,lnorm ,lgw};

  float L  = CoherentFitUtils::Langaus(x,lpar);
  float G = gnorm*TMath::Gaus(*x,gmu,gsigma,true);

  return L+G;
  
}

//********************************************
Double_t CoherentFitUtils::LangausPlusGausPlusConstant(Double_t *x, Double_t *par) {
//********************************************

  double llw    = par[0];
  double lmpv   = par[1];
  double lnorm  = par[2];
  double lgw    = par[3];
  double gnorm  = par[4];
  double gmu    = par[5];
  double gsigma = par[6];
  double C      = par[7];

  double lpar[]  = {llw ,lmpv ,lnorm ,lgw};

  float L  = CoherentFitUtils::Langaus(x,lpar);
  float G = gnorm*TMath::Gaus(*x,gmu,gsigma,true);

  return L+G+C;
  
}

//********************************************
Double_t CoherentFitUtils::DoubleLangausPlusGausPlusConstant(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double blw    = par[4];
  double bmpv   = par[5];
  double bnorm  = par[6];
  double bgw    = par[7];
  double gnorm  = par[8];
  double gmu    = par[9];
  double gsigma = par[10];
  double C      = par[11];

  double spar[]  = {slw ,smpv ,snorm ,sgw};
  double bpar[]  = {blw ,bmpv ,bnorm ,bgw};

  float S  = CoherentFitUtils::Langaus(x,spar);
  float SS = gnorm*TMath::Gaus(*x,gmu,gsigma,true);
  float B  = CoherentFitUtils::Langaus(x,bpar);

  return S+SS+B+C;
  
}

//********************************************
Double_t CoherentFitUtils::LaguerreSum(Double_t *x, Double_t *par) {
//********************************************

  double xt = (x[0]-1)/(x[0]+1);
  double result = par[0]*ROOT::Math::laguerre(0,xt)
    + par[1]*ROOT::Math::laguerre(1,xt)
    + par[2]*ROOT::Math::laguerre(2,xt)
    + par[3]*ROOT::Math::laguerre(3,xt);
  
  return result;
}

//********************************************
TF1* CoherentFitUtils::LangausFit(TH1F* h, bool use_poisson){
//********************************************

  if(h->GetEntries() == 0){
    std::cout << "Empty histogram" << std::endl;
    return NULL;
  }

  TF1* f = new TF1("f",Langaus,h->GetBinCenter(1),h->GetBinCenter(h->GetNbinsX()),4);
  f->SetParameters(0.075,                              //landau width
		   h->GetBinCenter(h->GetMaximumBin()),//landau mpv
		   0.6,                                //normalization
		   0.08);                              //gaussian width
  f->SetParLimits(0,0.03,2);

  h->Fit("f","LNQ","R");

  // h->Draw();
  // gPad->Update();gPad->WaitPrimitive();
  
  return f;
}

//********************************************
TF1* CoherentFitUtils::LangausPlusConstantFit(TH1F* h, bool use_poisson){
//********************************************

  if(h->GetEntries() == 0){
    std::cout << "Empty histogram" << std::endl;
    return NULL;
  }

  TF1* f = new TF1("f",LangausPlusConstant,h->GetBinCenter(1),h->GetBinCenter(h->GetNbinsX()),5);
  f->SetParameters(0.075,                              //landau width
		   h->GetBinCenter(h->GetMaximumBin()),//landau mpv
		   0.8,                                //normalization
		   0.08,                               //gaussian width
		   0.01);                              //constant
  f->SetParLimits(0,0.03,2);
  f->SetParLimits(3,0   ,5);
  f->SetParNames("#sigma_{L,B}","MPV_{B}","N_{B}","#sigma_{G,B}","C");
  
  h->Fit("f","LNQ");
  h->Fit("f","WLQ");

  // h->Draw();
  // gPad->Update();gPad->WaitPrimitive();
  
  return f;
}

//********************************************
TF1* CoherentFitUtils::DoubleLangausFit(TH1F* h, bool use_poisson){
//********************************************

  if(h->GetEntries() == 0){
    std::cout << "Empty histogram" << std::endl;
    return NULL;
  }

  std::string fname(h->GetTitle());
  fname = fname+"_IncoherentFit";
  
  TF1* f = new TF1(fname.c_str(),DoubleLangaus,h->GetBinCenter(1),h->GetBinCenter(h->GetNbinsX()),8);
  f->SetParameters(0.1,                                //landau width
		   h->GetBinCenter(h->GetMaximumBin()),//landau mpv
		   0.4,                                //normalization
		   0.1,                                //gaussian width
		   0.1,                                //landau width
		   2.0,                                //landau mpv
		   0.4,                                //normalization
		   0.1);                               //gaussian width

  f->SetParLimits(0,0.003,1); //landau with below this value gives shaky behaviour
  f->SetParLimits(3,0   ,5); //landau with must be positive
  f->SetParLimits(4,0.03,1); //landau with below this value gives shaky behaviour
  f->SetParLimits(7,0   ,5); //landau with must be positive

  f->SetParNames("#sigma_{L,S}","MPV_{S}","N_{S}","#sigma_{G,S}","#sigma_{L,RE}","MPV_{RE}","N_{RE}","#sigma_{G,RE}");
  
  //h->Fit(fname.c_str(),"LNQ");
  // h->Draw();
  // gPad->Update();gPad->WaitPrimitive();

   h->Fit(fname.c_str(),"WLQ");
  // h->Draw();
  // gPad->Update();//gPad->WaitPrimitive();
  
  return f;
}

//********************************************
TF1* CoherentFitUtils::GausnFit(TH1F* h, bool use_poisson){
//********************************************

  if(h->GetEntries() == 0){
    std::cout << "Empty histogram" << std::endl;
    return NULL;
  }

  TF1* f = new TF1("f","gausn(0)",h->GetBinCenter(1),h->GetBinCenter(h->GetNbinsX()));
  f->SetParameters(0.2,                                //norm
		   h->GetBinCenter(h->GetMaximumBin()),//mean
		   h->GetRMS());                       //sigma
  f->SetParLimits(1,h->GetBinCenter(1),h->GetBinCenter(h->GetNbinsX()));
  f->SetParLimits(2,0,10);
  if(use_poisson)h->Fit("f","QWL");
  else h->Fit("f","Q");

  return f;
}

//********************************************
void CoherentFitUtils::GetABCParametrization(double &A, double &B, double &C,
					     std::vector<std::pair<double,double>> X,
					     std::vector<std::pair<double,double>> Y,
					     bool equal_weights, bool draw_plot){
//********************************************

  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);
  //DeleteOutliers(x,x_error,y,y_error);
  
  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);

  if(draw_plot)tg->Draw("ap");
  
  TF1* f = new TF1("f","[A]*pow(x,[B])+[C]",3,60);
  f->SetParameters(10,-0.4,0);

  if(equal_weights)tg->Fit("f","WQ");
  else             tg->Fit("f","Q");

  A = f->GetParameter(0);
  B = f->GetParameter(1);
  C = f->GetParameter(2);

  tg->Draw("ap");
  gPad->Update();gPad->WaitPrimitive();
  if(draw_plot){gPad->Update();gPad->WaitPrimitive();}
}

//********************************************
void CoherentFitUtils::GetABCDRParametrization(double &A, double &B, double &C, double &D, double &R,
					       std::vector<std::pair<double,double>> X,
					       std::vector<std::pair<double,double>> Y,
					       bool equal_weights, bool draw_plot){
//********************************************

  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);
  //DeleteOutliers(x,x_error,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","[A]*pow(x+[R],[B])*(log(x)+[C])+[D]",3,60);
  f->SetParameters(-0.294,-0.597,-56.1,1.07,1.86);
  if(draw_plot)tg->Draw("ap");
  tg->Fit("f","Q");

  A = f->GetParameter("A");
  B = f->GetParameter("B");
  C = f->GetParameter("C");
  D = f->GetParameter("D");
  R = f->GetParameter("R");

  tg->Draw("ap");
  gPad->Update();gPad->WaitPrimitive();
  if(draw_plot){gPad->Update();gPad->WaitPrimitive();}
}

//********************************************
void CoherentFitUtils::GetLaguerreParametrization(double &l1, double &l2, double &l3, double &l4,
						  double &error_l1, double &error_l2, double &error_l3, double &error_l4,
						  std::vector<std::pair<double,double>> X,
						  std::vector<std::pair<double,double>> Y){
//********************************************

  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f",CoherentFitUtils::LaguerreSum,1,61,4);
  if(l4 == 0)f->FixParameter(3,0);

  tg->Fit("f","WWQ");
  tg->Fit("f","Q");
  
  l1 = f->GetParameter(0);
  l2 = f->GetParameter(1);
  l3 = f->GetParameter(2);
  l4 = f->GetParameter(3);
  error_l1 = f->GetParError(0);
  error_l2 = f->GetParError(1);
  error_l3 = f->GetParError(2);
  error_l4 = f->GetParError(3);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetLinearParametrization(double &l1, double &l2,
						double &error_l1, double &error_l2,
						std::vector<std::pair<double,double>> X,
						std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","pol1(0)",1,61);

  tg->Fit("f","WWQ");
  tg->Fit("f","Q");
  
  l1 = f->GetParameter(0);
  l2 = f->GetParameter(1);
  error_l1 = f->GetParError(0);
  error_l2 = f->GetParError(1);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetParabolicParametrization(double &l1, double &l2,
						   double &error_l1, double &error_l2,
						   std::vector<std::pair<double,double>> X,
						   std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);
  
  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","[0]*pow(x+[1],2)",1,61);
  
  tg->Fit("f","WWQ");
  tg->Fit("f","Q");
  
  l1 = f->GetParameter(0);
  l2 = f->GetParameter(1);
  error_l1 = f->GetParError(0);
  error_l2 = f->GetParError(1);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetDiLogParametrization(double &l1, double &l2, double &l3,
					       double &error_l1, double &error_l2, double &error_l3,
					       std::vector<std::pair<double,double>> X,
					       std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","[0]+[1]*TMath::DiLog(x+[2])",1,61);
  f->SetParameters(2.16533e+00,7.37390e-02,8.85421e-01);

  tg->Fit("f","WWQ");
  tg->Fit("f","Q");
  
  l1 = f->GetParameter(0);
  l2 = f->GetParameter(1);
  l3 = f->GetParameter(2);
  error_l1 = f->GetParError(0);
  error_l2 = f->GetParError(1);
  error_l3 = f->GetParError(2);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetSignalLWEstimation(double &p1, double &p2,
					     double &error_p1, double &error_p2,
					     std::vector<std::pair<double,double>> X,
					     std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","([0]/x-1)/x+[1]",1,61);
  f->SetParameters(4.3,0.09);

  tg->Fit("f","WWQ");
  
  p1 = f->GetParameter(0);
  p2 = f->GetParameter(1);
  error_p1 = f->GetParError(0);
  error_p2 = f->GetParError(1);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetSignalMPVEstimation(double &p1, double &p2, double &p3,
					      double &error_p1, double &error_p2, double &error_p3,
					      std::vector<std::pair<double,double>> X,
					      std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  //TF1* f = new TF1("f","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  TF1* f = new TF1("f","[0]*(([1]*x-1)/([1]*x+1)+[2])",1,61);
  f->SetParameters(-4.09,0.16,-1.4);

  tg->Fit("f","WWQ");
  
  p1 = f->GetParameter(0);
  p2 = f->GetParameter(1);
  p3 = f->GetParameter(2);
  error_p1 = f->GetParError(0);
  error_p2 = f->GetParError(1);
  error_p3 = f->GetParError(2);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
void CoherentFitUtils::GetSignalGWEstimation(double &p1, double &p2,
					     double &error_p1, double &error_p2,
					     std::vector<std::pair<double,double>> X,
					     std::vector<std::pair<double,double>> Y){
//********************************************
  
  std::vector<double> x, x_error, y, y_error;
  SeparatePair(X,x,x_error);
  SeparatePair(Y,y,y_error);

  int n = std::min((int)x.size(),(int)y.size());
  TGraphErrors* tg = new TGraphErrors(n,&x[0],&y[0],&x_error[0],&y_error[0]);
  tg->SetMarkerStyle(20);
  
  TF1* f = new TF1("f","[0]/(x+1)+[1]",1,61);
  f->SetParameters(3.86,0.07);

  tg->Fit("f","WWQ");
  
  p1 = f->GetParameter(0);
  p2 = f->GetParameter(1);
  error_p1 = f->GetParError(0);
  error_p2 = f->GetParError(1);

  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************
double CoherentFitUtils::ComputeLikelihood(TH1F* h, TF1* f, double integral){
//********************************************
  
  double Likelihood = 0;
  for(int i = 0; i < h->GetNbinsX(); i++){
    double hvalue = h->GetBinContent(i+1);
    if(hvalue*integral<1)continue;
    double fvalue = f->Eval(h->GetBinCenter(i+1));
    Likelihood = Likelihood + log(ROOT::Math::poisson_pdf(hvalue*integral,fvalue*integral));
  }

  return Likelihood;
}

//********************************************
double CoherentFitUtils::NormRegularization(const std::vector<double> norm, const std::vector<double> integral){
//********************************************
  
  if(norm.empty())return 0;

  double regularization = 0;

  /*double mean = 0;
  for(int i = 0; i < (int)norm.size(); i++)mean += norm[i];
  mean = mean/norm.size();

  double stddev = 0;
  for(int i = 0; i < (int)norm.size(); i++)stddev += pow(norm[i]-mean,2);
  stddev = stddev/(norm.size()-1);
  stddev = sqrt(stddev);

  for(int i = 0; i < (int)norm.size(); i++)regularization += pow(integral[i]*(norm[i]-mean),2);///stddev,2);*/

  for(int i = 1; i < (int)norm.size()-1; i++){
    double enorm = (norm[i-1]+norm[i+1])/2;
    double diff  = pow(integral[i]*(norm[i]-enorm),2);
    double error = sqrt(integral[i]*(norm[i]+enorm)/2)/integral[i];
    regularization += diff / (2*error);
  }

  return regularization;
}

//********************************************
double CoherentFitUtils::AlphaRegularization(const std::vector<double> alpha, const std::vector<double> fw){
//********************************************
  
  if(alpha.empty())return 0;

  double regularization = 0;

  for(int i = 1; i < (int)alpha.size()-1; i++){
    double ealpha = (alpha[i-1]+alpha[i+1])/2;
    double diff  = pow((alpha[i]-ealpha),2);
    double error = pow(fw[i],2);
    regularization += diff / error;
  }

  return regularization;
}

//********************************************
Double_t CoherentFitUtils::ABCParametrization(Double_t *x, Double_t *par){
//********************************************

  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t Q = par[3];

  double result = 0;
  if(Q == 0)result = A*pow(x[0],B)+C;
  else      result = sqrt(pow(A*pow(x[0],B)+C,2)+pow(Q,2));
  
  return result;
}

//********************************************
Double_t CoherentFitUtils::QuadraticABCParametrization(Double_t *x, Double_t *par){
//********************************************

  Double_t A1 = par[0];
  Double_t B1 = par[1];
  Double_t C1 = par[2];
  Double_t A2 = par[3];
  Double_t B2 = par[4];
  Double_t C2 = par[5];

  double r1 = A1*pow(x[0],B1)+C1;
  double r2 = A2*pow(x[0],B2)+C2;
  
  return sqrt(pow(r1,2)+pow(r2,2));
}

//********************************************
Double_t CoherentFitUtils::ABCDRParametrization(Double_t *x, Double_t *par){
//********************************************

  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t D = par[3];
  Double_t R = par[4];
  Double_t S = par[5];

  Double_t xx = x[0]+S;
  
  return A*pow(xx+R,B)*(log(xx)+C)+D;
}

//********************************************
Double_t CoherentFitUtils::ABCDRDerivativeA(Double_t *x, Double_t *par){
//********************************************

  Double_t B = par[1];
  Double_t C = par[2];
  Double_t R = par[4];
  Double_t S = par[5];

  Double_t xx = x[0]+S;
  
  return pow(xx+R,B)*(log(xx)+C);
}

//********************************************
Double_t CoherentFitUtils::ABCDRDerivativeB(Double_t *x, Double_t *par){
//********************************************

  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t R = par[4];
  Double_t S = par[5];

  Double_t xx = x[0]+S;
  
  return A*(log(xx)+C)*log(xx+R)*pow(xx+R,B);
}

//********************************************
Double_t CoherentFitUtils::ABCDRDerivativeC(Double_t *x, Double_t *par){
//********************************************

  Double_t A = par[0];
  Double_t B = par[1];
  Double_t R = par[4];
  Double_t S = par[5];

  Double_t xx = x[0]+S;
  
  return A*pow(xx+R,B);
}

//********************************************
Double_t CoherentFitUtils::ABCDRDerivativeD(Double_t *x, Double_t *par){
//********************************************

  return 1.;
}

//********************************************
Double_t CoherentFitUtils::ABCDRDerivativeR(Double_t *x, Double_t *par){
//********************************************

  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t R = par[4];
  Double_t S = par[5];

  Double_t xx = x[0]+S;
  
  return A*B*pow(xx+R,B-1)*(log(xx)+C);
}

//********************************************
void CoherentFitUtils::SeparatePair(std::vector<std::pair<double,double>> p,
				    std::vector<double> &v1, std::vector<double> &v2){
//********************************************

  v1.clear();
  v2.clear();
  
  for(int i = 0; i < (int)p.size(); i++){
    v1.push_back(p[i].first);
    v2.push_back(p[i].second);
  }
}

//********************************************
void CoherentFitUtils::DeleteOutliers(std::vector<double> &x, std::vector<double> &x_error,
				      std::vector<double> &y, std::vector<double> &y_error){
//********************************************

  int i = 1;
  int n = std::min((int)x.size(),(int)y.size());
  while(i < n-1){
    double tc1 = abs((y[i]-y[i-1])/(x[i]-x[i-1]));
    double tc2 = abs((y[i]-y[i+1])/(x[i]-x[i+1]));
    if(tc2>100*tc1){
      x.erase(x.begin()+i+1);
      x_error.erase(x_error.begin()+i+1);
      y.erase(y.begin()+i+1);
      y_error.erase(y_error.begin()+i+1);
    }
    else i++;
  }
}

//********************************************
bool CoherentFitUtils::IsTrueSignal(AnaParticlePD* part){
//********************************************

  bool ItIs = false;
  
  if(!part)return ItIs;

  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return ItIs;

  if(truePart->PDG == 321 && truePart->ProcessEnd == 2)ItIs = true;

  return ItIs;
}

//********************************************
TGraphErrors* CoherentFitUtils::GetGraph(std::vector<std::pair<double,double>> x, std::vector<std::pair<double,double>> y){
//********************************************

  std::vector<double> xv, xev, yv, yev;
  SeparatePair(x,xv,xev);
  SeparatePair(y,yv,yev);

  return new TGraphErrors(xv.size(),&xv[0],&yv[0],&xev[0],&yev[0]);
}
