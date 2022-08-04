{

  // load input files 
  gROOT->ProcessLine((".L " + string(gSystem->Getenv("PIONANALYSISROOT")) + "/macros/load.C").c_str());
  
  bool save_plots=false;
  bool wait_primi=true;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetLegendTextSize(0.);

  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1);

  
  /*
  gStyle->SetPaperSize(26,26);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.13);
  */
  
  //draw seltrk_length after initial cuts
  //first with no cuts
  draw.SetLegendSize(0.2,0.005);  
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",45,0,500,"beamparticle","","","AREA");

  if(save_plots)c1->Print("plot1.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //after BeamPionCut
  draw.SetTitle("After BeamPionCut");
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",45,0,500,"beamparticle","accum_level[0][0]>0","","AREA");

  if(save_plots)c1->Print("plot2.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //after BeamIsTrackCut
  draw.SetTitle("After BeamIsTrackCut");
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",45,0,500,"beamparticle","accum_level[0][0]>1","","AREA");
  draw.DrawCutLineVertical(226,true,"r");
  
  if(save_plots)c1->Print("plot3.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //after BeamGeometricCut

  draw.SetTitle("After BeamPionGeometricCut");
  draw.SetTitleX("Beam Particle End Z Position [cm]");
  draw.Draw(d,mc,"seltrk_endpos[2]",50,0,500,"beamparticle","accum_level[0][0]>2","","AREA");
  draw.DrawCutLineVertical(226,true,"r");
  
  if(save_plots)c1->Print("plot4.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //after BeamPionEndsAPA3
  draw.SetTitle("After PionEndAPA3Cut");
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",40,0,500,"beamparticle","accum_level[0][0]>3","","AREA");

  if(save_plots)c1->Print("plot5.pdf");
  if(wait_primi)c1->WaitPrimitive();


  //look at number of daughters
  draw.SetTitle("After PionEndAPA3Cut");
  draw.SetTitleX("Number of Daughters");
  draw.Draw(d,mc,"seltrk_ndau/2",10,0,10,"beamparticle","accum_level[0][0]>3","","AREA");

  if(save_plots)c1->Print("plot6.pdf");
  if(wait_primi)c1->WaitPrimitive();

  
  //look at daughters track CNN score (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetLegendPos("c");
  draw.SetTitleX("Daughters Track Score");
  draw.Draw(d,mc,"seltrk_dau_CNNscore[][0]",50,0,1,"dauparticle","accum_level[0][0]>3 && seltrk_dau_chi2_ndf!=-999","","AREA NOLEG");
  draw.DrawCutLineVertical(0.35,true,"r");
  
  if(save_plots)c1->Print("plot7.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at track-like daughters distance to vertex (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetLegendPos();
  draw.SetTitleX("Daughters Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",100,0,100,"dauparticle","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  draw.DrawCutLineVertical(10,true,"l");
  
  if(save_plots)c1->Print("plot8.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at track-like daughters distance to vertex logaritmic scale (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetMinY(0.1);
  draw.SetLogY();
  draw.SetTitleX("Daughters Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",100,0,100,"dauparticle","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  draw.DrawCutLineVertical(10,true,"l");
  
  if(save_plots)c1->Print("plot9.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at track-like daughters chi2proton distribution (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetMinY(0);
  draw.SetLogY(0);
  draw.SetTitleX("Daughters chi^{2}_{prot}/ndf");
  draw.Draw(d,mc,"seltrk_dau_chi2_prot/seltrk_dau_chi2_ndf",40,0,400,"dauparticle","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  draw.DrawCutLineVertical(50,true,"r");
  
  if(save_plots)c1->Print("plot10.pdf");
  if(wait_primi)c1->WaitPrimitive();





  //look at daughters track CNN score (only for shower hypothesis: seltrk_dau_chi2_ndf==-999)
  draw.SetTitle("After NoPionDaugtherCut");
  draw.SetLegendPos("r");
  draw.SetTitleX("Daughters Track Score");
  draw.Draw(d,mc,"seltrk_dau_CNNscore[][0]",50,0,1,"dauparticle","accum_level[0][0]>4 && seltrk_dau_chi2_ndf==-999","","AREA");
  draw.DrawCutLineVertical(0.35,true,"l");
  
  if(save_plots)c1->Print("plot11.pdf");
  if(wait_primi)c1->WaitPrimitive();



  draw.SetTitle("CEX branch");
  
  //look at shower-like daughters distance to vertex (only for shower hypothesis: seltrk_dau_chi2_ndf==-999)
  draw.SetTitleX("Dauther Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",50,0,200,"dauparticle","accum_level[0][0]>4 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(2,  true,"r");
  draw.DrawCutLineVertical(100,true,"l");
  
  if(save_plots)c1->Print("plot12.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daughters distance to vertex in log scale (only for shower hypothesis: seltrk_dau_chi2_ndf==-999)
  draw.SetMinY(0.1);
  draw.SetLogY();
  draw.SetTitleX("Dauther Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",175,0,350,"dauparticle","accum_level[0][0]>4 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(2,  true,"r");
  draw.DrawCutLineVertical(100,true,"l");

  
  if(save_plots)c1->Print("plot13.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daugthers number of hits
  draw.SetMinY(0);
  draw.SetLogY(0);
  draw.SetTitleX("Daughter Number of hits");
  draw.Draw(d,mc,"seltrk_dau_nhits",150,0,1500,"dauparticle","accum_level[0][0]>4 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_vtxdistance>2 && seltrk_dau_vtxdistance<100 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(12,  true,"r");
  draw.DrawCutLineVertical(1000,true,"l");

  
  if(save_plots)c1->Print("plot14.pdf");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daugthes number of hits in log scale
  draw.SetMinY(0.1);
  draw.SetLogY();
  draw.SetTitleX("Dauther Number of hits");
  draw.Draw(d,mc,"seltrk_dau_nhits",200,0,2000,"dauparticle","accum_level[0][0]>4 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_vtxdistance>2 && seltrk_dau_vtxdistance<100 & seltrk_dau_chi2_ndf!=-999","","AREA NOLEG");

  
  if(save_plots)c1->Print("plot15.pdf");
  if(wait_primi)c1->WaitPrimitive();



  //------------ custom category plots --------------
  draw.SetStackFillStyle(3001);

  draw.SetStackFillStyle(1001);
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",40,0,500,"pionana","accum_level[0][0]>3","","AREA");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetLegendSize(0.3,0.001);
  draw.SetLegendPos("c");  
  draw.SetTitleX("Track Score");
  draw.Draw(d,mc,"seltrk_dau_CNNscore[][0]",50,0,1,"daupionana","accum_level[0][0]>3 && seltrk_dau_chi2_ndf!=-999","","AREA");
  if(wait_primi)c1->WaitPrimitive();
  
}
