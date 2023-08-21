void DrawSyst_length_kaon(){

  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-04-04.root";
  //std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_all_0__2023-04-15.root";
  //std::string fmc   = dir+"systematics/syst_noWeights.root";
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

  //binning for the plot
  // const int nbins = 8;
  // double edges[nbins+1] = {0,30,60,100,150,200,300,500,700};
  const int nbins = 19;
  double edges[nbins+1] = {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,500,700};

  //get MC nominal
  TH1F* h_mc_nominal = new TH1F("h_mc_nominal","h_mc_nominal",nbins,edges);
  syst->Project("h_mc_nominal","bestcandidate_length_toy","(accum_level>2 && bestcandidate_chi2_kaon_perndf_toy<50)*weight_syst_total*toy_weight");
  //set statistical error → sqrt(counts)
  for(int i = 0; i < h_mc_nominal->GetNbinsX(); i++){
    h_mc_nominal->SetBinError(i+1,sqrt(h_mc_nominal->GetBinContent(i+1)));
    std::cout << h_mc_nominal->GetBinContent(i+1) << " " << h_mc_nominal->GetBinError(i+1) << std::endl;
  }
  h_mc_nominal->Draw("e");

  //clone nominal to add systematic
  TH1F* h_mc_syst = (TH1F*)h_mc_nominal->Clone("h_mc_syst");

  //get relative error plot
  draw.DrawRelativeErrors(syst,"bestcandidate_length_toy",nbins,edges,"accum_level>2 && bestcandidate_chi2_kaon_perndf_toy<50","","SYS");
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
  TH1F* h_data = new TH1F("h_data","h_data",nbins,edges);
  h_data->SetMarkerStyle(20);
  d->Project("h_data","bestcandidate_length","accum_level[0][]>6 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50");

  double integral_data = h_data->GetSumOfWeights();
  double integral_mc   = h_mc_nominal->GetSumOfWeights();

  h_data->SetTitle("");
  h_mc_nominal->SetTitle("");
  h_mc_syst->SetTitle("");

  h_mc_nominal->Sumw2(true);
  h_mc_nominal->Scale(integral_data/integral_mc);

  h_mc_syst->Sumw2(true);
  h_mc_syst->Scale(integral_data/integral_mc);

  h_data->GetXaxis()->SetTitle("Candidate's Length [cm]");
  h_mc_syst->GetXaxis()->SetTitle("Candidate's Length [cm]");
  h_mc_syst->GetXaxis()->CenterTitle();
  h_data->GetXaxis()->CenterTitle();

  gStyle->SetOptStat("i");
  h_data->Draw("e");
  gPad->Update();
  TPaveStats* ps = (TPaveStats*)h_data->FindObject("stats");
  ps->SetTextAlign(31);
  ps->SetX1NDC(0.7);
  ps->SetX2NDC(0.9);
  ps->SetY1NDC(0.925);
  ps->SetY2NDC(0.97);
  ps->SetBorderSize(1);
  gPad->Update();
  h_mc_syst->Draw("e2 same");
  h_mc_nominal->Draw("e2 same");
  h_data->Draw("esame");

  //draw legend now
  TLegend* lg = new TLegend(0.52,0.52,0.85,0.85);
  lg->AddEntry(h_data,"Data","ple");
  lg->AddEntry(h_mc_nominal,"MC stat","f");
  lg->AddEntry(h_mc_syst,"MC stat + syst","f");
  lg->Draw("same");

  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis(); gPad->Update();//always redraw axis
  gPad->Print("SYST_final_kaon_length.pdf");
}
