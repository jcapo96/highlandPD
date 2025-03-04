{

  // load input files 
  gROOT->ProcessLine((".L " + string(gSystem->Getenv("SECONDARYKAONANALYSISROOT")) + "/macros/tesis/selection/load.C").c_str());

  gStyle->SetOptStat(0);
  
  const int nhits = 300;

  TH2D* h[nhits];// 

  for(int i = 0; i < nhits; i++){
    std::stringstream ssi;
    ssi << i;
  
    h[i] = new TH2D(("h"+ssi.str()+"").c_str(),("h"+ssi.str()+"").c_str(),50,0,50,50,0,10);
    mc->Draw(("bestcandidate_hit_dedx["+ssi.str()+"]:bestcandidate_hit_resrange["+ssi.str()+"]>>h"+ssi.str()+"(50,0,50,50,0,10)").c_str(),
	     "accum_level[0][]>8","col z");
    h[i] = (TH2D*)gPad->GetPrimitive(("h"+ssi.str()+"").c_str());
    std::cout << "HOLA" << std::endl;
    if(i > 0)h[0]->Add(h[i]);
    std::cout << i << std::endl;
  }
  h[0]->GetXaxis()->SetTitle("Residual range [cm]");
  h[0]->GetYaxis()->SetTitle("dEdx [MeV/cm]");
  //h[0]->GetZaxis()->SetRangeUser(5,60);
  h[0]->Draw("colz");

}
