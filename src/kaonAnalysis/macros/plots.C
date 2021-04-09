{

  // load input files 
  gROOT->ProcessLine((".L " + string(gSystem->Getenv("KAONANALYSISROOT")) + "/macros/load.C").c_str());
  
  bool save_plots=false;
  bool wait_primi=true;

  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetLegendTextSize(0.);

  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.3);


  //------------ custom category plots --------------
  //  draw.SetStackFillStyle(3001);
  draw.SetStackFillStyle(1001);
  draw.SetMarkerSize(0.6);


  //multiplicity and efficiency plots
  /*  draw.SetTitleX("True secondary kaons");
  draw.Draw(mc,"seltrk_truedaukaons",4,0,4,"particle","seltrk_truepdg==211");
  if(wait_primi)c1->WaitPrimitive();

  draw.Draw(mc,"seltrk_truedaukaons",3,1,4,"particle","seltrk_truepdg==211");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("Reco daughter true PDG");
  draw.Draw(mc,"seltrk_dau_truepdg",1,321,322,"daughter","seltrk_truepdg==211");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("True terciary muons");
  draw.Draw(mc,"seltrk_truedaukaon_nmu",3,1,4,"particle","seltrk_truepdg==211");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetTitleX("Reco gdaughter true PDG");
  draw.Draw(mc,"seltrk_gdau_truepdg",1,-13,-12,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg==321");
  if(wait_primi)c1->WaitPrimitive();

  
  //secondary kaons true info
  draw.SetTitleX("true mom (GeV)");
  draw.Draw(mc,"seltrk_dau_truemom",100,0,6,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg==321");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("true end mom (GeV)");
  draw.Draw(mc,"seltrk_dau_trueendmom",100,0,2,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg==321");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("true end process");
  draw.Draw(mc,"seltrk_dau_trueendproc",10,0,10,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg==321");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("true end mom (GeV)");
  draw.SetLineWidth(3);
  draw.Draw(mc,"seltrk_dau_trueendmom",100,0,2,"all","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==2","hist","","decaying kaons");
  draw.Draw(mc,"seltrk_dau_trueendmom",100,0,2,"all","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==3","same hist","","interacting kaons");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("Length (cm)");
  draw.SetLineWidth(1);
  draw.Draw(mc,"seltrk_dau_length",40,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg==321");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("Length (cm)");
  draw.SetLineWidth(3);
  draw.Draw(mc,"seltrk_dau_length",40,0,500,"all","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==2","hist","","decaying kaons");
  draw.Draw(mc,"seltrk_dau_length",40,0,500,"all","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==3","hist same","","interacting kaons");
  if(wait_primi)c1->WaitPrimitive();
  
  //terciary muons true info
  draw.SetTitleX("true mom (GeV)");
  draw.SetLineWidth(1);
  draw.Draw(mc,"seltrk_gdau_truemom",100,0,1,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_gdau_truepdg==-13");
  if(wait_primi)c1->WaitPrimitive();
  */
  draw.SetTitleX("true end mom (GeV)");
  draw.SetLegendPos("l");
  draw.Draw(mc,"seltrk_gdau_trueendmom",100,0,0.3,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_gdau_truepdg==-13");
  if(wait_primi)c1->WaitPrimitive();
  /*
  draw.SetTitleX("Length (cm)");
  draw.SetLegendPos("r");
  draw.Draw(mc,"seltrk_gdau_length",100,0,400,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_gdau_truepdg==-13");
  if(wait_primi)c1->WaitPrimitive();

  
  //daughter distribution
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",100,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22");
  if(wait_primi)c1->WaitPrimitive();

  //gd distribution
  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",100,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22");
  if(wait_primi)c1->WaitPrimitive();

  
  //first cut, kaon daughters
  draw.SetTitleX("Daughter daughters");
  draw.SetLineWidth(3);
  draw.Draw(mc,"seltrk_dau_ndau",5,0,5,"all","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==2","hist","","decaying kaons");
  draw.DrawCutLinesVertical(1,2);
  if(wait_primi)c1->WaitPrimitive();

  //distributions after cuts
  draw.SetLineWidth(1);
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",100,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",100,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1");
  if(wait_primi)c1->WaitPrimitive();
  

  //second cut, muon daughters
  draw.SetTitleX("GD daughters");
  draw.Draw(mc,"seltrk_gdau_ndau",5,0,5,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_gdau_truepdg==-13");
  draw.DrawCutLineVertical(3,true);
  if(wait_primi)c1->WaitPrimitive();
  
  //distributions after cuts
  draw.SetLineWidth(1);
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",100,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2");
  if(wait_primi)c1->WaitPrimitive();

  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",100,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2");
  if(wait_primi)c1->WaitPrimitive();
  

  //third cut, muons are track-like objects
  draw.SetTitleX("GD reco object type");
  draw.SetLegendPos("l");
  draw.Draw(mc,"seltrk_gdau_type",3,0,3,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2");
  draw.DrawCutLineVertical(2,true,"r");
  if(wait_primi)c1->WaitPrimitive();

  //distributions after cuts
  draw.SetLineWidth(1);
  draw.SetLegendPos("r");
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",100,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",100,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2");
  if(wait_primi)c1->WaitPrimitive();
  

  //fourth cut, CNN score
  draw.SetTitleX("GD CNN track score");
  draw.SetLegendPos("c");
  draw.Draw(mc,"seltrk_gdau_CNNscore[][][0]",40,0,1,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2");
  draw.DrawCutLineVertical(0.6,true,"r");
  if(wait_primi)c1->WaitPrimitive();

  //distributions after cuts
  draw.SetLineWidth(1);
  draw.SetLegendPos("r");
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",100,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2 && seltrk_gdau_CNNscore[][][0]>0.6");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",100,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2 && seltrk_gdau_CNNscore[][][0]>0.6");
  if(wait_primi)c1->WaitPrimitive();
  

  //fifth cut, chi2_muon
  draw.SetTitleX("GD #chi^{2}_{#mu}/ndf");
  draw.SetLegendPos("r");
  draw.Draw(mc,"seltrk_gdau_chi2_muon/seltrk_gdau_chi2_ndf",80,0,80,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2 && seltrk_gdau_CNNscore[][][0]>0.6");
  draw.DrawCutLineVertical(6,true,"l");
  if(wait_primi)c1->WaitPrimitive();

    //distributions after cuts
  draw.SetLineWidth(1);
  draw.SetLegendPos("r");
  draw.SetTitleX("Daughter length (cm)");
  draw.Draw(mc,"seltrk_dau_length",40,0,500,"daughter","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2 && seltrk_gdau_CNNscore[][][0]>0.6 && seltrk_gdau_chi2_muon/seltrk_gdau_chi2_ndf<6");
  if(wait_primi)c1->WaitPrimitive();
  
  draw.SetTitleX("GD Length (cm)");
  draw.Draw(mc,"seltrk_gdau_length",40,0,500,"gdaughtermuon","seltrk_truepdg==211 && seltrk_dau_truepdg!=22 && seltrk_dau_ndau==1 && seltrk_gdau_ndau<=2 && seltrk_gdau_type==2 && seltrk_gdau_CNNscore[][][0]>0.6 && seltrk_gdau_chi2_muon/seltrk_gdau_chi2_ndf<6");
  if(wait_primi)c1->WaitPrimitive();*/
}
