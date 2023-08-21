void trueeff(){
  
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/BrokenTracks_data_1GeV_new.root";
  std::string fmc   = dir+"mc/prod4a/BrokenTracks_trueff.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  // TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

  //drawing tools
  DrawingTools draw(fmc.c_str());
  draw.ChangeToProtoDUNEStyle(); //set default protodune style

  //Style specific stuff
  gStyle->SetErrorX(0.5);  //I like to have errors in X, at least for the selection
  draw.SetMarkerSize(0.9); //put correct marker size, I think protodune style is not working on this

  TGaxis::SetMaxDigits(3); //no more than two digits per axis. Too few?

  // gStyle->SetLabelOffset(0.009, "x"); //reset to nominal
  // gStyle->SetLabelOffset(0.009, "y"); //reset to nominal
  TGaxis::SetExponentOffset(0, 0.015, "y"); //exponent offset
  
  gStyle->SetOptStat(1000000);
  gStyle->SetStatBorderSize(1);
  draw.SetStatPos(-999,1-0.075+0.05); //top of the pad - margin for plot + width of statsbox + width of the frame
  draw.ChangeCategoryColor("broken",2,2); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("brokendaughter",2,2); //change kaon colors for candidates plots

  double default_legend_x_size = 0.3;

  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();

  draw.SetStatPos(0.9-0.2,1-0.075+0.05);
  TH1F* den = new TH1F("denominator","denominator",30,220,600);
  den->GetXaxis()->SetTitle("Candidate Start Z [cm]");
  den->GetXaxis()->CenterTitle();
  den->GetYaxis()->SetTitle("# Events/(12.67 cm)");
  den->GetYaxis()->CenterTitle();
    
  mc->Project("denominator","bestcandidate_pos[2]","bestcandidate_parent_pos[2]<220 && bestcandidate_parent_endpos[2]>220");
  den->SetMarkerColor(98);
  den->SetLineColor(98);
  den->Draw("e");gPad->Update();

  draw.SetStatPos(0.9,1-0.075+0.05);
  TH1F* num = new TH1F("numerator","numerator",30,220,600);
  mc->Project("numerator","bestcandidate_pos[2]","bestcandidate_parent_pos[2]<220 && bestcandidate_parent_endpos[2]>220 && bestcandidate_pos[2]<234 && abs(bestcandidate_truepos[2]-227)>7");
  num->SetMarkerColor(9);
  num->SetLineColor(9);
  num->Draw("sames e");

  gPad->BuildLegend(0.6,0.6,0.85,0.85,"","lep");

  tt1.DrawLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update();

  gPad->Print("SYST_brokentracks_trueeff.pdf");
}
