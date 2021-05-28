void rms(){
  
  //the input file name
  std::string filename= "/data4/DUNE/migue/analysis/6GeV_prod4_3.root";

  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* mc    = (TTree*)rfile->Get("ana");

  //set seltrk_ndau variable
  int seltrk_ndau;
  mc->SetBranchAddress("seltrk_ndau", &seltrk_ndau);

  //two canvas
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(3,1);

  int nentries = mc->GetEntries();

  for(int i = 0; i < 10000; i++){
    mc->GetEntry(i);
    c1->cd();
    mc->Draw("seltrk_truedaukaons","seltrk_truepdg==211 && seltrk_truedaukaons>0","",1,i);
    if(mc->GetSelectedRows()>0){
      for(int j = 0; j < seltrk_ndau; j++){
	std::stringstream ssj;
	ssj << j;
	/*mc->Draw(("seltrk_dau_hit_z["+ssj.str()+"]:seltrk_dau_hit_y["+ssj.str()+"]:seltrk_dau_hit_x["+ssj.str()+"]").c_str(),
		 ("seltrk_dau_hit_x["+ssj.str()+"]!=-999 && seltrk_dau_truepdg["+ssj.str()+"]==321 && seltrk_dau_trueendproc["+ssj.str()+"]==2 && seltrk_dau_ndau["+ssj.str()+"]>0 && seltrk_dau_trueendmom["+ssj.str()+"]<0.4 && seltrk_dau_hit_resrange["+ssj.str()+"]<14").c_str(),
		 "",1,i);*/
	mc->Draw(("seltrk_dau_hit_z["+ssj.str()+"]:seltrk_dau_hit_y["+ssj.str()+"]:seltrk_dau_hit_x["+ssj.str()+"]").c_str(),
	  ("seltrk_dau_hit_x["+ssj.str()+"]!=-999 && seltrk_dau_truepdg["+ssj.str()+"]==321 && seltrk_dau_trueendproc["+ssj.str()+"]!=2  && seltrk_dau_hit_resrange["+ssj.str()+"]<14").c_str(),"",1,i);
	if(mc->GetSelectedRows()>0){
	  std::cout << i << " " << j << std::endl;
	  TGraph2D* tg = new TGraph2D(mc->GetSelectedRows(),
				      mc->GetV3(),mc->GetV2(),mc->GetV1());
	  
	  tg->SetMarkerStyle(20);
	  tg->GetXaxis()->SetTitle("X");
	  tg->GetYaxis()->SetTitle("Y");
	  tg->GetZaxis()->SetTitle("Z");

	  c2->cd(1);
	  tg->Draw("lp");
	  gPad->SetTheta(90);
	  gPad->SetPhi(0);

	  c2->cd(2);
	  tg->Draw("lp");
	  gPad->SetTheta(0);
	  gPad->SetPhi(0);
	  
	  c2->cd(3);
	  tg->Draw("lp");
	  gPad->SetTheta(0);
	  gPad->SetPhi(-90);
	  
	  c2->Update();
	  c1->cd();c1->WaitPrimitive();
	}
      }
    }
    if(i%10==0)std::cout << i << " entries read" << std::endl;
  }
  
}
