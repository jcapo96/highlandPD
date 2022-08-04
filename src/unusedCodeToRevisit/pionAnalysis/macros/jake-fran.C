{

  // reproduce plots in: https://docs.google.com/presentation/d/1g_v7DDgnkpo2PUtflCcyxU_fKrGMjAPJyxcuja8WrDE/edit#slide=id.g6de773c846_2_0

  
  // load input files 
  gROOT->ProcessLine((".L " + string(gSystem->Getenv("PIONANALYSISROOT")) + "/macros/load.C").c_str());
  
  bool save_plots=false;
  bool wait_primi=true;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetLegendTextSize(0.);

  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1);

  


  //------------ custom category plots --------------
  //  draw.SetStackFillStyle(3001);
  draw.SetStackFillStyle(1001);

  /*
  draw.SetTitle("Page 26: Problem in data");  
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",40,0,500,"pionana2","accum_level[0][0]>1","","AREA");
  if(wait_primi)c1->WaitPrimitive();
  */

  draw.SetTitle("Page 27: OK");  
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",40,0,500,"pionana2","accum_level[0][0]>2","","AREA");
  if(wait_primi)c1->WaitPrimitive();
    
  draw.SetTitle("Page 28: OK");  
  draw.SetTitleX("Beam Particle Track Length [cm]");
  draw.Draw(d,mc,"seltrk_length",40,0,500,"pionana2","accum_level[0][0]>3","","AREA");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetLegendSize(0.3,0.001);
  draw.SetLegendPos("c");
  draw.SetTitle("Page 30: OK");  
  draw.SetTitleX("Track Score");
  draw.Draw(d,mc,"seltrk_dau_CNNscore[][0]",50,0,1,"daupionana","accum_level[0][0]>3 && seltrk_dau_chi2_ndf!=-999","","AREA");
  if(wait_primi)c1->WaitPrimitive();

  //look at track-like daughters distance to vertex (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetLegendPos();
  draw.SetTitle("Page 31-Top: Some bins differ");  
  draw.SetTitleX("Daughters Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",10,0,20,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  draw.DrawCutLineVertical(10,true,"l");
  if(wait_primi)c1->WaitPrimitive();

  //look at track-like daughters distance to vertex logaritmic scale (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetMinY(1.3);
  draw.SetLogY();
  draw.SetTitle("Page 31-Bottom: Most bins differ");  
  draw.SetTitleX("Daughters Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",50,0,100,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  draw.DrawCutLineVertical(10,true,"r");
  if(wait_primi)c1->WaitPrimitive();

  //look at proton-like chi2/ndf for daughters (only for track hypothesis: seltrk_dau_chi2_ndf!=-999)
  draw.SetLegendPos("r");
  draw.SetMinY(0);
  draw.SetLogY(false);
  draw.SetTitle("Page 32: Small difference in data and MC");  
  draw.SetTitleX("Daughters chi^{2}_{prot}/ndf");
  draw.Draw(d,mc,"seltrk_dau_chi2_prot/seltrk_dau_chi2_ndf",40,0,400,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA");
  draw.DrawCutLineVertical(50,true,"r");
  if(wait_primi)c1->WaitPrimitive();
  

  
  //look at shower-like daughters distance to vertex (only for shower hypothesis: seltrk_dau_chi2_ndf==-999)
  draw.SetTitle("Page 34-Top: OK ");
  draw.SetTitleX("Dauther Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",50,0,100,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(2,  true,"r");
  draw.DrawCutLineVertical(100,true,"l");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daughters distance to vertex in log scale (only for shower hypothesis: seltrk_dau_chi2_ndf==-999)
  draw.SetMinY(0.5);
  draw.SetLogY();
  draw.SetTitle("Page 34-Bottom: very small differences in MC ");
  draw.SetTitleX("Dauther Distance To Vertex (cm)");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",175,0,350,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(2,  true,"r");
  draw.DrawCutLineVertical(100,true,"l");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daugthers number of hits
  draw.SetMinY(0);
  draw.SetLogY(0);
  draw.SetTitle("Page 35-Top: Not yet understood");
  draw.SetTitleX("Daughter Number of hits");
  draw.Draw(d,mc,"seltrk_dau_nhits",50,0,400,"daupionana","accum_level[0][0]>3 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_vtxdistance>2 && seltrk_dau_vtxdistance<100 && seltrk_dau_chi2_ndf==-999","","AREA NOLEG");
  draw.DrawCutLineVertical(12,  true,"r");
  draw.DrawCutLineVertical(1000,true,"l");
  if(wait_primi)c1->WaitPrimitive();

  //look at shower-like daugthes number of hits in log scale
  draw.SetMinY(0.1);
  draw.SetLogY();
  draw.SetTitle("Page 35-Bottom: Not yet understood");
  draw.SetTitleX("Daugther Number of hits");
  draw.Draw(d,mc,"seltrk_dau_nhits",200,0,2000,"daupionana","accum_level[0][0]>5 && seltrk_dau_CNNscore[][0]<0.35 && seltrk_dau_vtxdistance>2 && seltrk_dau_vtxdistance<100 & seltrk_dau_chi2_ndf!=-999","","AREA NOLEG");
  if(wait_primi)c1->WaitPrimitive();
  
}
