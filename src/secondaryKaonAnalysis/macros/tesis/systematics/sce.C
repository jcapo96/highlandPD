void sce(){

  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  //std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_sce_2023-04-12.root";
  ///std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_sceAllVoxels_2023-05-05.root";
  std::string fmc   = dir+"systematics/sce.root";

  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* syst = (TTree*)_file1->Get("all_syst");
  TTree* truth = (TTree*)_file1->Get("truth");

    //drawing tools
  DrawingTools draw(fmc.c_str());
  draw.ChangeToProtoDUNEStyle(); //set default protodune style

  //Style specific stuff
  gStyle->SetErrorX(0.5);  //I like to have errors in X, at least for the selection
  draw.SetMarkerSize(0.9); //put correct marker size, I think protodune style is not working on this
  TGaxis::SetMaxDigits(4); //no more than two digits per axis. Too few?

  // gStyle->SetLabelOffset(0.009, "x"); //reset to nominal
  // gStyle->SetLabelOffset(0.009, "y"); //reset to nominal
  TGaxis::SetExponentOffset(0, 0.015, "y"); //exponent offset
  
  gStyle->SetOptStat(1000000);
  gStyle->SetStatBorderSize(1);
  draw.SetStatPos(-999,1-0.075+0.05); //top of the pad - margin for plot + width of statsbox + width of the frame
  draw.ChangeCategoryColor("bestcandidateparticle",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidateparticle",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidateparticlereduced",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("particle",321,217);          //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidatedaughter",321,217);     //change kaon colors for candidates plots

  //text to be drawn
  TLatex tt1,tt2;
  tt1.SetNDC();
  tt2.SetNDC();
  tt2.SetTextAlign(31);  
  //binning
  double bins[] = {0,10,20,40,60,70,80,90,100,150,250};

  //plot
  draw.CenterTitles();
  draw.SetTitleX("Candidate #chi^{2}_{K}");
  draw.DrawRelativeErrors(syst,"bestcandidate_chi2_kaon_perndf_toy",10,bins,"accum_level>2","","SYS");

  TH1D* h = (TH1D*)gPad->GetPrimitive("alle_3");
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetNdivisions(503);

  h->SetLineWidth(3);
  h->Draw();

  tt1.DrawLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.94,"MC");
  gPad->RedrawAxis(); gPad->Update();//always redraw axis
  gPad->Print("syst_sce_rel.pdf");

  draw.SetTitleX("Candidate #chi^{2}_{K}");
  draw.DrawRelativeErrors(syst,"bestcandidate_chi2_kaon_perndf_toy",1,0,250,"accum_level>2","","SYS");

  // h = (TH1D*)gPad->GetPrimitive("alle_3");
  // h->GetXaxis()->CenterTitle();
  // h->GetXaxis()->SetTitleOffset(1.1);
  // h->GetYaxis()->CenterTitle();
  // h->GetYaxis()->SetTitleOffset(1.1);
  // h->GetYaxis()->SetNdivisions(503);

  // h->SetLineWidth(3);
  // h->Draw();

  // tt1.DrawLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // tt2.DrawLatex(0.9,0.94,"MC");
  // gPad->RedrawAxis(); gPad->Update();//always redraw axis
  // gPad->Print("syst_dqdx_rel.pdf");

  draw.CenterTitles();
  draw.SetTitleX("Candidate #chi^{2}_{K}");
  draw.DrawRelativeErrors(syst,"bestcandidate_Zstart_toy",30,0,600,"accum_level>2","","SYS");
}
