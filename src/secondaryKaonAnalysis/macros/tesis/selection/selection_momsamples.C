//#include "../../../../../../highland/src/highland2/highlandTools/src/SetProtoDUNEStyle.H"

bool save_plots = true;

void selection_momsamples(){

  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-04-04.root";
  std::string fmc   = dir+"mc/prod4a/6-7GeV_prod4a_microtree_2023-04-04.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
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
  draw.ChangeCategoryColor("bestcandidateparticle",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidateparticle",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidateparticlereduced",321,217); //change kaon colors for candidates plots
  draw.ChangeCategoryColor("particle",321,217);          //change kaon colors for candidates plots
  draw.ChangeCategoryColor("candidatedaughter",321,217);     //change kaon colors for candidates plots
  
  double default_legend_x_size = 0.3;
  
  //text to be drawn
  TLatex tt1;
  tt1.SetNDC();
  
  draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.Draw(d,mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"candidateparticlereduced","accum_level[0][]>6 && beam_nominal_mom==6","","area pur ignoreempty");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("selection_result_chi2_kaon_mom6.pdf");		       

  // draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  // draw.SetLegendPos("r");
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.Draw(d,mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"candidateparticlereduced","accum_level[0][]>6 && beam_nominal_mom==7","","area pur ignoreempty");
  // tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update(); //always redraw axis
  // if(save_plots)gPad->Print("selection_result_chi2_kaon_mom7.pdf");		       

  // draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  // draw.SetLegendPos("r");
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.Draw(d,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==6","","pur ignoreempty","6 GeV");
  // draw.Draw(d,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==7","same","pur ignoreempty","7 GeV");
  // tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update(); //always redraw axis
  // if(save_plots)gPad->Print("selection_result_chi2_kaon_mom_data.pdf");		       

  // draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  // draw.SetLegendPos("r");
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.Draw(mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==6","","pur ignoreempty","6 GeV");
  // draw.Draw(mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==7","same","pur ignoreempty","7 GeV");
  // tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update(); //always redraw axis
  // if(save_plots)gPad->Print("selection_result_chi2_kaon_mom_mc.pdf");		       

  draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.Draw(d,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==6","","pur ignoreempty area1","6 GeV");
  draw.Draw(d,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==7","same","pur ignoreempty area1","7 GeV");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  TH1F* dummy = (TH1F*)gPad->GetPrimitive("_all_18_19");
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("selection_result_chi2_kaon_mom_data_1.pdf");		       
  
  // draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  // draw.SetLegendPos("r");
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.Draw(mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==6","","pur ignoreempty area1","6 GeV");
  // draw.Draw(mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"all","accum_level[0][]>6 && beam_nominal_mom==7","same","pur ignoreempty area1","7 GeV");
  // tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // TH1F* dummy = (TH1F*)gPad->GetPrimitive("_all_18_19");
  // dummy->GetXaxis()->CenterTitle();
  // dummy->GetYaxis()->CenterTitle();
  // gPad->RedrawAxis();gPad->Update(); //always redraw axis
  // if(save_plots)gPad->Print("selection_result_chi2_kaon_mom_mc_1.pdf");		       
 }
