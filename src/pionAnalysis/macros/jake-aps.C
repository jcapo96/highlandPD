{

  // reproduce plots in: https://indico.fnal.gov/event/23989/contribution/0/material/slides/0.pdf

  
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

  //slide 9. We have less events than he. We are missing cosmics in the first bins. Maybe he is plotting pandora tracks before applying beam cuts?
  draw.SetTitle("Slide 9: NOT OK");  
  draw.SetTitleX("Track End Z [cm]");
  draw.Draw(d,mc,"seltrk_endpos[2]",40,0,500,"pionana2","accum_level[0][0]>0","","AREA");
  if(wait_primi)c1->WaitPrimitive();
    
  //slide 10, center. Less events than he and we are also missing cosmics.
  draw.SetLegendPos("l");
  draw.SetTitle("Slide 10, center: NOT OK");  
  draw.SetTitleX("#DeltaZ [cm]");
  draw.Draw(d,mc,"seltrk_pos[2]-beam_endpos[2]",40,0,50,"pionana2","accum_level[0][0]>0","","AREA");
  if(wait_primi)c1->WaitPrimitive();
  
  //slide 10, right.
  draw.SetLegendPos("r");
  draw.SetTitle("Slide 10, right: OK");  
  draw.SetTitleX("#Chi^{2}/ndf");
  draw.Draw(d,mc,"seltrk_chi2_prot/seltrk_chi2_ndf",100,0,400,"pionana2","accum_level[0][0]>3","","AREA");
  if(wait_primi)c1->WaitPrimitive();

  //slide 11, clean sample. Equal in statistics. Difference in legend,he has joint pi+- from interactions.
  draw.SetTitle("Slide 11: OK");  
  draw.SetTitleX("Track End Z [cm]");
  draw.Draw(d,mc,"seltrk_endpos[2]",100,0,400,"pionana2","accum_level[0][0]>4","","AREA");
  if(wait_primi)c1->WaitPrimitive();

  //slide 12, left. We don't have first bin. Statistics are not equal. Legend has changed
  draw.SetTitle("Slide 12, left: NOT OK");  
  draw.SetTitleX("CNN track-like score");
  draw.Draw(d,mc,"seltrk_dau_CNNscore[][0]",50,0,1,"daupionana","accum_level[0][0]>4 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  if(wait_primi)c1->WaitPrimitive();

  //slide 12, right. Statistics are not equal. Legend has changed
  draw.SetTitle("Slide 12, right: NOT OK");  
  draw.SetTitleX("#Chi^{2}/ndf");
  draw.Draw(d,mc,"seltrk_dau_chi2_prot/seltrk_dau_chi2_ndf",40,0,400,"daupionana","accum_level[0][0]>4 && seltrk_dau_chi2_CNNscore[][0]>0.35 && seltrk_dau_chi2_ndf!=-999","","AREA  NOLEG");
  if(wait_primi)c1->WaitPrimitive();

  //slide 13, left. Statistics are not ok
  draw.SetTitle("Slide 13, left: NOT OK");  
  draw.SetTitleX("Distance To Vertex [cm]");
  draw.Draw(d,mc,"seltrk_dau_vtxdistance",50,0,100,"daupionana","accum_level[][0]>4 && seltrk_dau_chi2_ndf!=-999","","area");
  if(wait_primi)c1->WaitPrimitive();

  //slide 13, right. Statistics are not ok
  draw.SetTitle("Slide 13, left: NOT OK");  
  draw.SetTitleX("Number of hits");
  draw.Draw(d,mc,"seltrk_dau_nhits2",48,0,400,"daupionana","accum_level[][0]>5","","area");
  if(wait_primi)c1->WaitPrimitive();
}
