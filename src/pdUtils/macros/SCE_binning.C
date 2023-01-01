//short macro to produce plaint text files with the information of the SCE maps


void CreateFile(TFile* rfile, std::string filename, std::string &histoname){

  TH3F* h = (TH3F*)rfile->Get(histoname.c_str());

  int nbinsx = h->GetNbinsX();
  int nbinsy = h->GetNbinsY();
  int nbinsz = h->GetNbinsZ();

  double xmin = h->GetXaxis()->GetXmin();
  double ymin = h->GetYaxis()->GetXmin();
  double zmin = h->GetZaxis()->GetXmin();

  double xmax = h->GetXaxis()->GetXmax();
  double ymax = h->GetYaxis()->GetXmax();
  double zmax = h->GetZaxis()->GetXmax();

  filename = filename.substr(filename.find("SCE"),filename.find(".root")-filename.find("SCE"));

  //open plain text file
  ofstream myfile;
  myfile.open(filename+"_"+histoname+".dat");

  //write binning
  myfile << nbinsx << "\t" << xmin << "\t" << xmax << "\n";
  myfile << nbinsy << "\t" << ymin << "\t" << ymax << "\n";
  myfile << nbinsz << "\t" << zmin << "\t" << zmax << "\n";

  //write bin content
  for(int x = 1; x < nbinsx+1; x++)
    for(int y = 1; y < nbinsy+1; y++)
      for(int z = 1; z < nbinsz+1; z++)
	myfile << x << "\t" << y << "\t" << z << "\t" << h->GetBinContent(x,y,z) << "\n";

  //close output file
  myfile.close();

}

void SCE_binning(){

  //open root file 
  std::string filename = "../data/SCE_DataDriven_180kV_v4.root";
  TFile* rfile = TFile::Open(filename.c_str());
  
  if(!rfile){
    std::cout << "cannot open file " << filename << std::endl;
    std::exit(1);
  }
  
  std::string histoname[18] = {"RecoFwd_Displacement_X_Neg",
			       "RecoFwd_Displacement_Y_Neg",
			       "RecoFwd_Displacement_Z_Neg",
			       "RecoBkwd_Displacement_X_Neg",
			       "RecoBkwd_Displacement_Y_Neg",
			       "RecoBkwd_Displacement_Z_Neg",
			       "Reco_ElecField_X_Neg",
			       "Reco_ElecField_Y_Neg",
			       "Reco_ElecField_Z_Neg",
			       "RecoFwd_Displacement_X_Pos",
			       "RecoFwd_Displacement_Y_Pos",
			       "RecoFwd_Displacement_Z_Pos",
			       "RecoBkwd_Displacement_X_Pos",
			       "RecoBkwd_Displacement_Y_Pos",
			       "RecoBkwd_Displacement_Z_Pos",
			       "Reco_ElecField_X_Pos",
			       "Reco_ElecField_Y_Pos",
			       "Reco_ElecField_Z_Pos"};
 
  for(int i = 0; i < 18; i++)CreateFile(rfile,filename,histoname[i]);

}
