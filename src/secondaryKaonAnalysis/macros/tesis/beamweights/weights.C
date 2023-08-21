void weights(){
  
  //load trees
  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/6-7GeV_prod4a_reco2_microtree_2023-04-04.root";
  std::string fmc   = dir+"mc/prod4a/6-7GeV_prod4a_microtree_2023-04-04.root";

  TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());

  TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  TTree* truth = (TTree*)_file1->Get("truth");

  //data beam particles
  int n_total_6_data = d->Draw("0.5","beam_nominal_mom==6 && beam_npdgs>0","goff");
  int n_total_7_data = d->Draw("0.5","beam_nominal_mom==7 && beam_npdgs>0","goff");

  int n_unknown_6_data = d->Draw("0.5","beam_nominal_mom==6 && beam_npdgs==0","goff");
  int n_unknown_7_data = d->Draw("0.5","beam_nominal_mom==7 && beam_npdgs==0","goff");

  int n_pion_6_data = d->Draw("0.5","beam_nominal_mom==6 && (beam_pdgs[0]!=2212 && beam_pdgs[0]!=321)","goff");
  int n_pion_7_data = d->Draw("0.5","beam_nominal_mom==7 && (beam_pdgs[0]!=2212 && beam_pdgs[0]!=321)","goff");

  int n_proton_6_data = d->Draw("0.5","beam_nominal_mom==6 && beam_pdgs==2212","goff");
  int n_proton_7_data = d->Draw("0.5","beam_nominal_mom==7 && beam_pdgs==2212","goff");

  int n_kaon_6_data = d->Draw("0.5","beam_nominal_mom==6 && beam_pdgs==321","goff");
  int n_kaon_7_data = d->Draw("0.5","beam_nominal_mom==7 && beam_pdgs==321","goff");
  
  //std::cout << n_pion_6_data << " " << n_proton_6_data << " " << n_kaon_6_data << " " << n_unknown_6_data << " " << n_total_6_data << std::endl;
  //std::cout << n_pion_7_data << " " << n_proton_7_data << " " << n_kaon_7_data << " " << n_unknown_7_data << " " << n_total_7_data << std::endl;

  double r_pion_6_data = (double)n_pion_6_data/n_total_6_data;
  double r_pion_7_data = (double)n_pion_7_data/n_total_7_data;

  double r_proton_6_data = (double)n_proton_6_data/n_total_6_data;
  double r_proton_7_data = (double)n_proton_7_data/n_total_7_data;

  double r_kaon_6_data = (double)n_kaon_6_data/n_total_6_data;
  double r_kaon_7_data = (double)n_kaon_7_data/n_total_7_data;

  double r_pion_6_data_error = sqrt(1./n_pion_6_data+1./n_total_6_data);
  double r_pion_7_data_error = sqrt(1./n_pion_7_data+1./n_total_7_data);

  double r_proton_6_data_error = sqrt(1./n_proton_6_data+1./n_total_6_data);
  double r_proton_7_data_error = sqrt(1./n_proton_7_data+1./n_total_7_data);
  
  double r_kaon_6_data_error = sqrt(1./n_kaon_6_data+1./n_total_6_data);
  double r_kaon_7_data_error = sqrt(1./n_kaon_7_data+1./n_total_7_data);

  std::cout << "PERCENTAGES DATA" << std::endl;
  std::cout << "pion 6 GeV " << r_pion_6_data << " +- " << r_pion_6_data*r_pion_6_data_error << std::endl;
  std::cout << "pion 7 GeV " << r_pion_7_data << " +- " << r_pion_7_data*r_pion_7_data_error << std::endl;
  std::cout << "prot 6 GeV " << r_proton_6_data << " +- " << r_proton_6_data*r_proton_6_data_error << std::endl;
  std::cout << "prot 7 GeV " << r_proton_7_data << " +- " << r_proton_7_data*r_proton_7_data_error << std::endl;
  std::cout << "kaon 6 GeV " << r_kaon_6_data << " +- " << r_kaon_6_data*r_kaon_6_data_error << std::endl;
  std::cout << "kaon 7 GeV " << r_kaon_7_data << " +- " << r_kaon_7_data*r_kaon_7_data_error << std::endl;
  
  //mc beam particles
  int n_total_6_mc = mc->Draw("0.5","beam_nominal_mom==6","goff");
  int n_total_7_mc = mc->Draw("0.5","beam_nominal_mom==7","goff");

  int n_pion_6_mc = mc->Draw("0.5","beam_nominal_mom==6 && (abs(beam_truepdg)==211 || abs(beam_truepdg)==13 || abs(beam_truepdg)==11)","goff");
  int n_pion_7_mc = mc->Draw("0.5","beam_nominal_mom==7 && (abs(beam_truepdg)==211 || abs(beam_truepdg)==13 || abs(beam_truepdg)==11)","goff");

  int n_proton_6_mc = mc->Draw("0.5","beam_nominal_mom==6 && abs(beam_truepdg)==2212","goff");
  int n_proton_7_mc = mc->Draw("0.5","beam_nominal_mom==7 && abs(beam_truepdg)==2212","goff");

  int n_kaon_6_mc = mc->Draw("0.5","beam_nominal_mom==6 && abs(beam_truepdg)==321","goff");
  int n_kaon_7_mc = mc->Draw("0.5","beam_nominal_mom==7 && abs(beam_truepdg)==321","goff");

  double r_pion_6_mc = (double)n_pion_6_mc/n_total_6_mc;
  double r_pion_7_mc = (double)n_pion_7_mc/n_total_7_mc;

  double r_proton_6_mc = (double)n_proton_6_mc/n_total_6_mc;
  double r_proton_7_mc = (double)n_proton_7_mc/n_total_7_mc;

  double r_kaon_6_mc = (double)n_kaon_6_mc/n_total_6_mc;
  double r_kaon_7_mc = (double)n_kaon_7_mc/n_total_7_mc;

  double r_pion_6_mc_error = sqrt(1./n_pion_6_mc+1./n_total_6_mc);
  double r_pion_7_mc_error = sqrt(1./n_pion_7_mc+1./n_total_7_mc);

  double r_proton_6_mc_error = sqrt(1./n_proton_6_mc+1./n_total_6_mc);
  double r_proton_7_mc_error = sqrt(1./n_proton_7_mc+1./n_total_7_mc);
  
  double r_kaon_6_mc_error = sqrt(1./n_kaon_6_mc+1./n_total_6_mc);
  double r_kaon_7_mc_error = sqrt(1./n_kaon_7_mc+1./n_total_7_mc);

  std::cout << "PERCENTAGES MC" << std::endl;
  std::cout << "pion 6 GeV " << r_pion_6_mc << " +- " << r_pion_6_mc*r_pion_6_mc_error << std::endl;
  std::cout << "pion 7 GeV " << r_pion_7_mc << " +- " << r_pion_7_mc*r_pion_7_mc_error << std::endl;
  std::cout << "prot 6 GeV " << r_proton_6_mc << " +- " << r_proton_6_mc*r_proton_6_mc_error << std::endl;
  std::cout << "prot 7 GeV " << r_proton_7_mc << " +- " << r_proton_7_mc*r_proton_7_mc_error << std::endl;
  std::cout << "kaon 6 GeV " << r_kaon_6_mc << " +- " << r_kaon_6_mc*r_kaon_6_mc_error << std::endl;
  std::cout << "kaon 7 GeV " << r_kaon_7_mc << " +- " << r_kaon_7_mc*r_kaon_7_mc_error << std::endl;
  
  //ratios
  double r_pion_6_data_mc = r_pion_6_data/r_pion_6_mc;
  double r_pion_7_data_mc = r_pion_7_data/r_pion_7_mc;

  double r_proton_6_data_mc = r_proton_6_data/r_proton_6_mc;
  double r_proton_7_data_mc = r_proton_7_data/r_proton_7_mc;

  double r_kaon_6_data_mc = r_kaon_6_data/r_kaon_6_mc;
  double r_kaon_7_data_mc = r_kaon_7_data/r_kaon_7_mc;

  double r_pion_6_data_mc_error = r_pion_6_data_mc*sqrt(pow(r_pion_6_data_error,2)+pow(r_pion_6_mc_error,2));
  double r_pion_7_data_mc_error = r_pion_7_data_mc*sqrt(pow(r_pion_7_data_error,2)+pow(r_pion_7_mc_error,2));

  double r_proton_6_data_mc_error = r_proton_6_data_mc*sqrt(pow(r_proton_6_data_error,2)+pow(r_proton_6_mc_error,2));
  double r_proton_7_data_mc_error = r_proton_7_data_mc*sqrt(pow(r_proton_7_data_error,2)+pow(r_proton_7_mc_error,2));

  double r_kaon_6_data_mc_error = r_kaon_6_data_mc*sqrt(pow(r_kaon_6_data_error,2)+pow(r_kaon_6_mc_error,2));
  double r_kaon_7_data_mc_error = r_kaon_7_data_mc*sqrt(pow(r_kaon_7_data_error,2)+pow(r_kaon_7_mc_error,2));

  std::cout << "RATIOS DATA/MC" << std::endl;
  std::cout << "pion 6 GeV " << r_pion_6_data_mc << " +- " << r_pion_6_data_mc_error << std::endl;
  std::cout << "pion 7 GeV " << r_pion_7_data_mc << " +- " << r_pion_7_data_mc_error << std::endl;
  std::cout << "prot 6 GeV " << r_proton_6_data_mc << " +- " << r_proton_6_data_mc_error << std::endl;
  std::cout << "prot 7 GeV " << r_proton_7_data_mc << " +- " << r_proton_7_data_mc_error << std::endl;
  std::cout << "kaon 6 GeV " << r_kaon_6_data_mc << " +- " << r_kaon_6_data_mc_error << std::endl;
  std::cout << "kaon 7 GeV " << r_kaon_7_data_mc << " +- " << r_kaon_7_data_mc_error << std::endl;


  //total beam weights
  double r_7_6_data = (double)n_total_7_data/n_total_6_data;
  double r_7_6_mc   = (double)n_total_7_mc/n_total_6_mc;
  
  double r_7_6_data_error = sqrt(1./n_total_7_data+1./n_total_6_data);
  double r_7_6_mc_error   = sqrt(1./n_total_7_mc+1./n_total_7_mc);

  double r_7_6_data_mc = r_7_6_data/r_7_6_mc;
  double r_7_6_data_mc_error = r_7_6_data_mc*sqrt(pow(r_7_6_data_error,2)+pow(r_7_6_mc_error,2));
  std::cout << "---------------" << std::endl;
  std::cout << "DATA TOTAL" << std::endl;
  std::cout << "6 GeV " << n_total_6_data << std::endl;
  std::cout << "7 GeV " << n_total_7_data << std::endl;
  std::cout << "ratio 7/6 " << r_7_6_data << " +/- " << r_7_6_data*r_7_6_data_error << std::endl; 
  std::cout << "MC TOTAL" << std::endl;
  std::cout << "6 GeV " << n_total_6_mc << std::endl;
  std::cout << "7 GeV " << n_total_7_mc << std::endl;
  std::cout << "ratio 7/6 " << r_7_6_mc << " +/- " << r_7_6_mc*r_7_6_mc_error << std::endl; 
  std::cout << "RATIO DATA MC " << std::endl;
  std::cout << r_7_6_data_mc << " +/- " << r_7_6_data_mc_error << std::endl;
}
