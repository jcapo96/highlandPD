void SystDiffError_AllSyst(){

  //open file
  std::string dir = "/dune/app/users/miagarc/technical_note/files/";
  std::string fmc   = dir+"systematics/syst_noBeamMomNorm.root";

  TFile *_file1 = TFile::Open(fmc.c_str());

  //get trees
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
  TGaxis::SetExponentOffset(0, 0.015, "y"); //exponent offset
  
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
  gPad->Print("plots/syst_rel_allsyst.pdf");

  draw.SetTitleX("Candidate #chi^{2}_{K}");
  draw.DrawRelativeErrors(syst,"bestcandidate_chi2_kaon_perndf_toy",1,0,250,"accum_level>2","","SYS");
}
