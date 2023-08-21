void broken_tracks(){

  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-02-02.root";
  std::string fmc   = dir+"mc/prod4a/6-7GeV_prod4a_microtree_2023-02-02.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");

  double theta0 = 0;
  double thetaf = 180;
  double step = 180;
  int nsteps = (thetaf-theta0)/step;
  double thetamin = theta0;
  double thetamax = theta0+step;

  TH1F* hmc = new TH1F("hmc","hmc",nsteps,theta0,thetaf);
  TH1F* hd  = new TH1F("hd" ,"hd" ,nsteps,theta0,thetaf);

  std::string tot_pos_cut = "seltrk_endpos[2]>220 && seltrk_endpos[2]<800";
  std::string sub_pos_cut = "seltrk_endpos[2]>220 && seltrk_endpos[2]<230";

  for(int i = 0; i < nsteps; i++){
    std::stringstream sstmin, sstmax;
    sstmin << thetamin;
    sstmax << thetamax;
    std::string angle_cut = "abs(TMath::ATan2(seltrk_dir[0],seltrk_dir[2])*180/TMath::Pi())>"+sstmin.str()+" && abs(TMath::ATan2(seltrk_dir[0],seltrk_dir[2])*180/TMath::Pi())<"+sstmax.str();

    //MC
    double tot = mc->Draw("seltrk_endpos[2]",(""+tot_pos_cut+" && "+angle_cut+"").c_str());
    double sub = mc->Draw("seltrk_endpos[2]",(""+sub_pos_cut+" && "+angle_cut+"").c_str());
    std::cout << sub/tot << " +/- " << sub/tot*sqrt(pow(sqrt(sub)/sub,2)+pow(sqrt(tot)/tot,2)) << std::endl;
    hmc->SetBinContent(i+1,sub/tot);

    //data
    tot = d->Draw("seltrk_endpos[2]",(""+tot_pos_cut+" && "+angle_cut+"").c_str());
    sub = d->Draw("seltrk_endpos[2]",(""+sub_pos_cut+" && "+angle_cut+"").c_str());
    hd->SetBinContent(i+1,sub/tot);
    std::cout << sub/tot << " +/- " << sub/tot*sqrt(pow(sqrt(sub)/sub,2)+pow(sqrt(tot)/tot,2)) << std::endl;
    thetamin += step;
    thetamax += step;
  }

  hmc->SetLineWidth(3);
  hmc->SetLineColor(2);
  hmc->Draw("");
  hd->SetLineWidth(3);
  hd->SetLineColor(1);
  hd->Draw("same");
  
  gPad->Clear();
  TRatioPlot* r = new TRatioPlot(hd,hmc);
  r->Draw();
  
}
