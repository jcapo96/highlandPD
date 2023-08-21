void DrawSyst_dEdx_dau(){

  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_nobranches_2023-04-18.root";
  //std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_all_dauhits_2023-04-19.root";
  std::string fmc   = dir+"systematics/syst_noBeamMomNorm.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

  //load drawing tools
  DrawingTools draw(fmc.c_str());

  //set ProtoDUNE style
  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  gStyle->SetErrorX(0.5);
  TGaxis::SetExponentOffset(0, 0.015, "y"); //exponent offset

  //binning for the plot
  const int nbins = 25;
  double edges[nbins+1] = {0,1,2,3,4,5,6,7,8,9,10};

  //get MC nominal
  TH1F* h_mc_nominal = new TH1F("h_mc_nominal","h_mc_nominal",nbins,0,10);
  syst->Project("h_mc_nominal","bestcandidate_dau_hit_dedx_toy","(accum_level>2 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50)*weight_syst_total*toy_weight");
  for(int i = 0; i < h_mc_nominal->GetNbinsX(); i++)h_mc_nominal->SetBinError(i+1,sqrt(h_mc_nominal->GetBinContent(i+1)));
  h_mc_nominal->Draw("e");gPad->Update();//gPad->WaitPrimitive();

  //clone nominal to add systematic
  TH1F* h_mc_syst = (TH1F*)h_mc_nominal->Clone("h_mc_syst");

  //get systematic error
  draw.DrawRelativeErrors(syst,"bestcandidate_dau_hit_dedx_toy",nbins,0,10,"accum_level>2 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50","","SYS");
  TH1D* h_mc_syst_rel = (TH1D*)gPad->GetPrimitive("alle_3");

  //loop over bins and add quadratically the systematic error
  for(int i = 0; i < nbins; i++){
    double st_error = h_mc_syst->GetBinError(i+1);
    double syst_error = h_mc_syst_rel->GetBinContent(i+1)*h_mc_syst->GetBinContent(i+1);
    double tot_error = sqrt(pow(st_error,2)+pow(syst_error,2));
    h_mc_syst->SetBinError(i+1,tot_error);
  }

  h_mc_syst->SetFillColor(38);
  h_mc_nominal->SetFillColor(9);
  h_mc_syst->SetMarkerSize(0);
  h_mc_nominal->SetMarkerSize(0);
  h_mc_syst->SetLineWidth(0);
  h_mc_nominal->SetLineWidth(0);

  //data
  TH1F* h_data = new TH1F("h_data","h_data",nbins,0,10);
  h_data->SetMarkerStyle(20);
  d->Project("h_data","bestcandidate_dau_hit_dedx","accum_level[0]>2 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<25");
  h_data->Draw();

  double integral_data = h_data->GetSumOfWeights();
  double integral_mc   = h_mc_nominal->GetSumOfWeights();

  h_data->SetTitle("");
  h_mc_nominal->SetTitle("");
  h_mc_syst->SetTitle("");

  h_mc_nominal->Sumw2(true);
  h_mc_nominal->Scale(integral_data/integral_mc);

  h_mc_syst->Sumw2(true);
  h_mc_syst->Scale(integral_data/integral_mc);

  h_mc_syst->GetXaxis()->SetTitle("Candidate dE/dx [MeV/cm]");
  h_mc_syst->GetXaxis()->CenterTitle();

  gStyle->SetOptStat("i");
  h_mc_syst->Draw("e2");
  gPad->Update();
  TPaveStats* ps = (TPaveStats*)h_mc_syst->FindObject("stats");
  ps->SetTextAlign(31);
  ps->SetX1NDC(0.7);
  ps->SetX2NDC(0.9);
  ps->SetY1NDC(0.925);
  ps->SetY2NDC(0.97);
  ps->SetBorderSize(1);
  gPad->Update();
  h_mc_nominal->Draw("e2 same");
  h_data->Draw("e same");

  //draw legend now
  TLegend* lg = new TLegend(0.52,0.52,0.85,0.85);
  lg->AddEntry(h_data,"Data","ple");
  lg->AddEntry(h_mc_nominal,"MC stat","f");
  lg->AddEntry(h_mc_syst,"MC stat + syst","f");
  lg->Draw("same");

  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis(); //gPad->Update();//always redraw axis
  gPad->Print("SYST_final_muon_dedx.pdf");
}
