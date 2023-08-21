void binning(){

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

  //MC calculation
  std::string cutxz[2] =    {"abs(atan2(seltrk_enddir[0],seltrk_enddir[2]))*180/TMath::Pi()<15","abs(atan2(seltrk_enddir[0],seltrk_enddir[2]))*180/TMath::Pi()>15"};
  std::string cutyz_mc[2] = {"abs(atan2(seltrk_enddir[1],seltrk_enddir[2]))*180/TMath::Pi()<17","abs(atan2(seltrk_enddir[1],seltrk_enddir[2]))*180/TMath::Pi()>17"};

  double den,num;
  std::vector<double> eff_mc, eff_mc_error;
  for(int ix = 0; ix < 2; ix++){
    for(int iy = 0; iy < 2; iy++){
      den = mc->Draw("seltrk_endpos[2]",("accum_level>3 && "+cutxz[ix]+" && "+cutyz_mc[iy]+"").c_str(),"goff");
      num = mc->Draw("seltrk_endpos[2]",("accum_level>5 && "+cutxz[ix]+" && "+cutyz_mc[iy]+"").c_str(),"goff");
      eff_mc.push_back(num/den);
      eff_mc_error.push_back(TEfficiency::ClopperPearson(den,num,0.683,true)-num/den);
      std::cout << num/den << " +/- " << num/den-TEfficiency::ClopperPearson(den,num,0.683,false) << " +/- " << TEfficiency::ClopperPearson(den,num,0.683,true)-num/den << std::endl;
    }
  }
  // std::cout << "MC EFF " << TMath::Mean(eff_mc.size(),&eff_mc[0]) << " +/- " << TMath::RMS(eff_mc.size(),&eff_mc[0]) << std::endl;

  //DATA
  std::string cutyz_d[2] = {"abs(atan2(seltrk_enddir[1],seltrk_enddir[2]))*180/TMath::Pi()<23","abs(atan2(seltrk_enddir[1],seltrk_enddir[2]))*180/TMath::Pi()>23"};

  std::vector<double> eff_d,eff_d_error;
  for(int ix = 0; ix < 2; ix++){
    for(int iy = 0; iy < 2; iy++){
      den = d->Draw("seltrk_endpos[2]",("accum_level>3 && "+cutxz[ix]+" && "+cutyz_d[iy]+"").c_str(),"goff");
      num = d->Draw("seltrk_endpos[2]",("accum_level>5 && "+cutxz[ix]+" && "+cutyz_d[iy]+"").c_str(),"goff");
      eff_d_error.push_back(TEfficiency::ClopperPearson(den,num,0.683,true)-num/den);
      std::cout << num/den << " +/- " << num/den-TEfficiency::ClopperPearson(den,num,0.683,false) << " +/- " << TEfficiency::ClopperPearson(den,num,0.683,true)-num/den << std::endl;
      eff_d.push_back(num/den);
    }
  }
  std::cout << "DATA EFF " << TMath::Mean(eff_d.size(),&eff_d[0]) << " +/- " << TMath::RMS(eff_d.size(),&eff_d[0]) << std::endl;
  
  //RATIO
  std::vector<double> r,r_error;
  for(int i = 0; i < 4; i++){
    r.push_back(eff_d[i]/eff_mc[i]);
    r_error.push_back(r[i]*sqrt(pow(eff_mc_error[i]/eff_mc[i],2)+pow(eff_d_error[i]/eff_d[i],2)));
    std::cout << r[i] << " +/- " << r_error[i] << std::endl;
  }

  std::cout << "RATIO " << TMath::Mean(r.size(),&r[0]) << " +/- " << TMath::RMS(r.size(),&r[0]) << std::endl;

  double mean = 0;
  double std_dev = 0;

  for(int i = 0; i < 4; i++)mean += r[i];
  mean /= 4;

  for(int i = 0; i < 4; i++)std_dev = pow(r[i]-mean,2);
  std_dev /= 3;
  std_dev = sqrt(std_dev);

  std::cout << mean << " +/- " << std_dev << " (" << std_dev/mean*100 << "%)" << std::endl;

}
