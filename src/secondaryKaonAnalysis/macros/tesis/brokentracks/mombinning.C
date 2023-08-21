void mombinning(){

  std::string hostname = gSystem->HostName();

  std::string dir = "/data4/DUNE/migue/analysis/files/";

  std::string fdata = dir+"data/BrokenTracks_data_1GeV.root";
  std::string fmc   = dir+"mc/prod4a/BrokenTracks_mc_1GeV.root";
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
  std::string mom_down[3] = {"0","200","400"};
  std::string mom_up[3]   = {"200","400","600"};

  double den,num;
  for(int i = 0; i < 3; i++){
    std::string cutmom = "sqrt(pow(sqrt(beam_mom*beam_mom+0.105*0.105)-seltrk_deposited_energy/1000,2)-0.105)*1000>"+mom_down[i]+" && sqrt(pow(sqrt(beam_mom*beam_mom+0.105*0.105)-seltrk_deposited_energy/1000,2)-0.105)*1000<"+mom_up[i];
    den = mc->Draw("seltrk_endpos[2]",("accum_level>3 && "+cutmom+"").c_str(),"goff");
    num = mc->Draw("seltrk_endpos[2]",("accum_level>5 && "+cutmom+"").c_str(),"goff");
    std::cout << num/den << std::endl;;
  }
}
