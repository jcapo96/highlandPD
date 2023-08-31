void dEdx_plot_2D(){

  gROOT->ProcessLine(".L ./protoDUNEStyle.C");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);
  
  //get trees
  std::string dir   = "/dune/app/users/miagarc/technical_note/files/";
  std::string fdata = dir+"data/data.root";
  std::string fmc   = dir+"mc/mc.root";

  TFile* d_file = TFile::Open(fdata.c_str(),"r");
  TFile* mc_file = TFile::Open(fmc.c_str(),"r");

  if(!d_file || !mc_file){
    std::cout << "no file found" << std::endl;
    std::exit(1);
  }

  TTree* d  = (TTree*)d_file->Get("ana");
  TTree* mc = (TTree*)mc_file->Get("ana");
  
  const int nhits = 300;

  TH2D* d_h2  = new TH2D("d_h2" ,"",60,0,100,60,0,12);
  TH2D* mc_h2 = new TH2D("mc_h2","",60,0,100,60,0,12);
  TH2D* hdummy;

  for(int i = 1; i < nhits-1; i++){
    std::stringstream ssi;
    ssi << i;

    //data
    d->Draw(("bestcandidate_hit_dedx["+ssi.str()+"]:bestcandidate_hit_resrange["+ssi.str()+"]>>h(60,0,100,60,0,12)").c_str(),
	    "accum_level[0][]>6","colz");
    hdummy = (TH2D*)gPad->GetPrimitive("h");
    d_h2->Add(hdummy);
    hdummy->Reset();
    
    //MC
    mc->Draw(("bestcandidate_hit_dedx["+ssi.str()+"]:bestcandidate_hit_resrange["+ssi.str()+"]>>h(60,0,100,60,0,12)").c_str(),
    	     "accum_level[0][]>6","colz");
    hdummy = (TH2D*)gPad->GetPrimitive("h");
    mc_h2->Add(hdummy);
    hdummy->Reset();
  }

  //draw histograms
  TLatex tt1;
  tt1.SetNDC();

  TLatex tt2;
  tt2.SetNDC();
  tt2.SetTextAlign(31);
  
  //draw data
  double max = d_h2->GetMaximum();
  for(int ibinx = 0; ibinx < d_h2->GetNbinsX(); ibinx++)
    for(int ibiny = 0; ibiny < d_h2->GetNbinsY(); ibiny++)d_h2->SetBinContent(ibinx+1,ibiny+1,d_h2->GetBinContent(ibinx+1,ibiny+1)/max);

  d_h2->GetXaxis()->SetTitle("Residual Range [cm]");
  d_h2->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  d_h2->GetZaxis()->SetTitle("Relative Frequency");
  d_h2->GetXaxis()->CenterTitle();
  d_h2->GetYaxis()->CenterTitle();
  d_h2->GetZaxis()->CenterTitle();
  
  d_h2->Draw("colz");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.94,"Data");
  gPad->RedrawAxis();
  //gPad->SetLogz();

  gPad->Update();//gPad->WaitPrimitive();
  gPad->Print("plots/CF_dedx_2d_data.pdf");
  std::cout << "data printed" << std::endl;

  //draw MC
  //normalize to largest value
  max = mc_h2->GetMaximum();
  for(int ibinx = 0; ibinx < mc_h2->GetNbinsX(); ibinx++)
    for(int ibiny = 0; ibiny < mc_h2->GetNbinsY(); ibiny++)mc_h2->SetBinContent(ibinx+1,ibiny+1,mc_h2->GetBinContent(ibinx+1,ibiny+1)/max);

  mc_h2->GetXaxis()->SetTitle("Residual Range [cm]");
  mc_h2->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  mc_h2->GetZaxis()->SetTitle("Relative Frequency");
  mc_h2->GetXaxis()->CenterTitle();
  mc_h2->GetYaxis()->CenterTitle();
  mc_h2->GetZaxis()->CenterTitle();
  
  mc_h2->Draw("colz");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.94,"MC");
  gPad->RedrawAxis();

  gPad->Print("plots/CF_dedx_2d_mc.pdf");
}
