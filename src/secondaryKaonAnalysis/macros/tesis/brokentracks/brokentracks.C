void brokentracks(){
  
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/BrokenTracks_data_1GeV_new.root";
  std::string fmc   = dir+"mc/prod4a/BrokenTracks_mc_1GeV_new.root";
  // std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-02-02.root";
  // std::string fmc   = dir+"systematics/6-7GeV_prod4a_microtree_syst_dQdxCal_2023-02-08.root";

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

  //**************************************
  //selection plots
  //**************************************
  
  //beam track end pos
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.SetTitleX("Beam Particle End Z [cm]");
  // draw.Draw(d,mc,"seltrk_endpos[2]",60,0,600,"beamparticle","accum_level>2","","area pur ignoreempty");
  // tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_controlsample.pdf");

  //beam track end pos
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.SetLegendEntryHeight(0.07);
  // draw.SetTitleX("Beam Particle End Z [cm]");
  // draw.Draw(d,mc,"seltrk_endpos[2]",50,210,260,"broken","accum_level>2","","area pur ignoreempty");
  // tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // draw.DrawCutLinesVertical(220,234,true);
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_aparegion.pdf");

  //cos
  // draw.SetLegendPos("l");
  // draw.SetMinY(10);
  // draw.SetLogY(1);
  // draw.SetTitleX("cos(Beam particle-daughter)");
  // draw.Draw(d,mc,"abs(seltrk_enddir[0]*seltrk_dau_dir[][0]+seltrk_enddir[1]*seltrk_dau_dir[][1]+seltrk_enddir[2]*seltrk_dau_dir[][2])",50,0,1,"broken","accum_level>4","","area pur ignoreempty");
  // tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // draw.DrawCutLineVertical(0.9,true,"r",0.8);
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_cos.pdf");

  //dis
  // draw.SetLegendPos("r");
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.SetLogY(0);
  // draw.SetMinY(0);
  // draw.SetTitleX("Distance Beam particle-daughter [cm]");
  // draw.Draw(d,mc,"sqrt(pow(seltrk_dau_pos[][0]-seltrk_endpos[0],2)+pow(seltrk_dau_pos[][1]-seltrk_endpos[1],2)+pow(seltrk_dau_pos[][2]-seltrk_endpos[2],2))",30,0,30,"broken","accum_level>4 && seltrk_dau_pos[][2]>200 && abs(seltrk_enddir[0]*seltrk_dau_dir[][0]+seltrk_enddir[1]*seltrk_dau_dir[][1]+seltrk_enddir[2]*seltrk_dau_dir[][2])>0.85","","area pur ignoreempty");
  // tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // draw.DrawCutLineVertical(13,true,"l",0.5);
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_dis.pdf");

  // //selection result
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.SetTitleX("Beam Particle End Z [cm]");
  // draw.Draw(d,mc,"seltrk_endpos[2]",50,210,260,"broken","accum_level>5","","area pur ignoreempty");
  // tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_endselection.pdf");

  //beam track end pos
  // draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  // draw.SetTitleX("#theta_{XZ} [#circ]");
  // draw.Draw(d,mc,"abs(atan2(seltrk_enddir[0],seltrk_enddir[2]))*180/TMath::Pi()",40,0,40,"broken","accum_level>5","","area pur ignoreempty");
  // tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  // gPad->RedrawAxis();gPad->Update();
  // gPad->Print("SYST_brokentracks_XZ.pdf");

  // //beam track end pos
  draw.SetLegendPos(0.57,0.88); //legend shouldnt be over the axis
  draw.SetTitleX("#theta_{YZ} [#circ]");
  draw.Draw(d,mc,"abs(atan2(seltrk_enddir[1],seltrk_enddir[2]))*180/TMath::Pi()",40,0,40,"broken","accum_level>5","","area pur ignoreempty");
  tt1.DrawLatex(0.18,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->RedrawAxis();gPad->Update();
  gPad->Print("SYST_brokentracks_YZ.pdf");
}
