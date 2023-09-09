//#include "../../../../../../highland/src/highland2/highlandTools/src/SetProtoDUNEStyle.H"

bool save_plots = true;

void selection(){

  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/dune/app/users/miagarc/technical_note/files/";

  std::string fdata = dir+"data/data.root";
  std::string fmc   = dir+"mc/mc.root";

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
  
  //**************************************
  //selection plots
  //**************************************
  
  //beam momentum
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.SetTitleX("BI momentum [GeV]");
  draw.Draw(d,mc,"beam_mom",40,4,9,"beamparticlereduced","accum_level[0][0]>0 && beam_nominal_mom==6","","area pur ignoreempty");
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis("fG"); gPad->Update();//always redraw axis
  if(save_plots)gPad->Print("plots/selection_BI_mom.pdf");
  
  //candidates per event
  draw.SetTitleX("Candidates per event");
  draw.Draw(d,mc,"ncandidates",10,0,10,"beamparticlereduced","accum_level[0][]>0","","area pur ignoreempty");
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_ncandidates.pdf");
  
  //daughter type
  draw.SetTitleX("Daughter Type");
  draw.SetLegendPos("l");
  draw.SetLegendSize(0.24,-999);
  draw.Draw(d,mc,"candidates_dau_type",3,0,3,"candidatedauparticlereduced","accum_level[0][]>1","","area pur ignoreempty");
  draw.DrawCutLinesVertical(2,3,true);
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  gPad->Update();
  if(save_plots)gPad->Print("plots/selection_daughter_type.pdf");
  
  //daughter chi2
  draw.SetTitleX("Daughter #chi^{2}_{muon}");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.SetLegendSize(default_legend_x_size,-999);
  draw.Draw(d,mc,"candidates_dau_chi2_muon/candidates_dau_chi2_ndf",20,0,20,"candidatedaumuonreduced","accum_level[0][]>2 && candidates_dau_chi2_muon>0 && (abs(candidates_endpos[][2]-230)>20 && abs(candidates_endpos[][2]-460)>20) && (abs(candidates_dau_pos[][2]-230)>20 && abs(candidates_dau_pos[][2]-460)>20)","","area pur ignoreempty nostat");
  draw.DrawCutLineVertical(6,true,"l");
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_daughter_chi2.pdf");
  
  //daughter mom range
  draw.SetTitleX("Daughter Momentum by Range [MeV]");
  draw.Draw(d,mc,"candidates_dau_mom_muon",30,0,0.6,"candidatedaumuonreduced","accum_level[0][]>3 && (abs(candidates_endpos[][2]-230)>20 && abs(candidates_endpos[][2]-460)>20) && (abs(candidates_dau_pos[][2]-230)>20 && abs(candidates_dau_pos[][2]-460)>20)","","area ignoreempty pur");
  draw.DrawCutLinesVertical(0.221,0.245,false);
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_daughter_mom.pdf");
  
  //candidate daughter angle
  draw.SetTitleX("cos(Candidate-Daughter)");
  draw.SetMinY(5);
  draw.SetLogY();
  draw.SetLegendPos("l");
  draw.Draw(d,mc,"candidates_cos_dau",50,-1,1,"candidateparticlereduced","accum_level[0][]>4","","area pur ignoreempty");
  draw.DrawCutLineVertical(0.64,true,"l");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_candidate_cos.pdf");

  draw.SetTitleX("cos(Candidate-Daughter)");
  draw.SetMinY(5);
  draw.SetLogY();
  draw.SetLegendPos("l");
  draw.Draw(d,mc,"candidates_cos_dau",50,-1,1,"candidateendprocess","accum_level[0][]>4","","area pur ignoreempty");
  draw.DrawCutLineVertical(0.64,true,"l");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_candidate_cos_endprocess.pdf");
  
  //candidate daughter distance
  draw.SetTitleX("Distance Candidate-Daughter [cm]");
  draw.SetMinY(1);
  draw.SetLogY();
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.Draw(d,mc,"candidates_distance_dau",50,0,100,"candidateparticlereduced","accum_level[0][]>5","","area pur ignoreempty");
  draw.DrawCutLineVertical(10,true,"l");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->SetLogx();
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  gPad->Update();
  if(save_plots)gPad->Print("plots/selection_candidate_dis.pdf");
  
  //**************************************
  //selection result
  //**************************************
  gPad->SetLogx(0);
  draw.SetLogY(0);
  draw.SetMinY(0);
  draw.SetTitleX("Candidate #chi^{2}_{K}  [cm]");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.Draw(d,mc,"candidates_chi2_kaon/candidates_chi2_ndf",20,0,150,"candidateparticlereduced","accum_level[0][]>6","","area pur ignoreempty");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update(); //always redraw axis
  if(save_plots)gPad->Print("plots/selection_result_chi2_kaon.pdf");
  
  //**************************************
  //efficiency and purity
  //**************************************
  gStyle->SetPadBottomMargin(0.17);
  draw.DrawPurVSCut(mc,"bestcandidate_truepdg==321");
  draw.DrawEffVSCut(truth,"truekaon_truepdg==321","",-1,-1,"same");

  
  //**************************************
  //other stuff plots
  //**************************************
  draw.SetTitleX("Candidate Initial Momentum [GeV/c]");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88);
  draw.Draw(d,mc,"sqrt(pow(bestcandidate_calE,2)+2*bestcandidate_calE*473.7)/1000",12,0,2,"bestcandidateparticle","accum_level[0][0]>6 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50","","area pur ignoreempty");
  draw.DrawCutLineVertical(0.340,false,"l");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update();
  if(save_plots)gPad->Print("plots/endselection_candidate_mom.pdf");

  draw.SetTitleX("Candidate Initial Momentum [GeV/c]");
  draw.SetLegendPos("r");
  draw.SetLegendPos(0.57,0.88);
  draw.Draw(d,mc,"sqrt(pow(bestcandidate_calE,2)+2*bestcandidate_calE*473.7)/1000",12,0,2,"beamparticle","accum_level[0][0]>6 && bestcandidate_chi2_kaon/bestcandidate_chi2_ndf<50","","area pur ignoreempty catname");
  //draw.DrawCutLineVertical(0.340,false,"l");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update();
  if(save_plots)gPad->Print("plots/endselection_candidate_mom_beampart.pdf");		       
}
