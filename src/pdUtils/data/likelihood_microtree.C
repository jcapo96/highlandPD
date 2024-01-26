// == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf 
const double rho = 1.39; // [g/cm3], density of LAr
const double K = 0.307075; // [MeV cm2 / mol]
const double A = 39.948; // [g / mol], atomic mass of Ar 
const double I = 188.0e-6; // [MeV], mean excitation energy
const double Me = 0.511; // [Mev], mass of electron
// == Parameters for the density correction
const double density_C = 5.2146; 
const double density_y0 = 0.2; 
const double density_y1 = 3.0; 
const double density_a = 0.19559;
const double density_k = 3.0; 

double mass;

ROOT::Math::VavilovAccurate vav;
TGraph* tg_fcn;
TGraph* tg_rr_ke;
TGraph* tg_dedx_rr;

double Density_Correction(double beta, double gamma){
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma); 
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C; 
  }
  else if (density_y < density_y0){
    this_delta = 0.; 
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }
  
  return this_delta; 
} 

double dEdx_Bethe_Bloch(double KE, double mass){ 
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = Density_Correction(beta, gamma);
  
  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2); 
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm] 
  
  return this_dEdx;
}

double Get_Wmax(double KE, double mass){ 
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  
  return Wmax; 
} 

double Get_Landau_xi(double KE, double dx, double mass){ 
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2); 
  return xi; 
}    

double dEdx_PDF_setting(double *x, double *par){ 
  
  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width} 
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3]; 
  double y = (x[0] - b) / a; 
  
  double this_vav = 0.;
  
  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y); 
    this_vav =this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1])); 
    this_vav =TMath::Gaus(y, mu, sigma); 
  }
  else{ // == Vavilov
    this_vav =vav.Pdf(y, par[0], par[1]);
    this_vav =this_vav / a;
  }
  
  return this_vav;
}


double LikelihoodFit(TGraph* tg){

  double width = 0.65;

  TF1* pdf = new TF1("pdf", dEdx_PDF_setting, -10., 20., 5);
  
  double L  = 0;
  double Lf = 20;
  double step = 0.2;
  std::vector<double> L_v,Likelihood_v;
  L_v.clear();
  Likelihood_v.clear();

  while(L<Lf){
    double likelihood = 0;
    for(int i = 0; i < tg->GetN(); i++){
      double range = tg->GetPointX(i)+L;
      double dEdx = tg->GetPointY(i);
      double ke = tg_rr_ke->Eval(range);
      double gamma = (ke/mass)+1.0; 
      double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
      double xi = Get_Landau_xi(ke, width, mass); 
      double Wmax = Get_Wmax(ke, mass); 
      double kappa = xi / Wmax; 
      double dEdx_BB = dEdx_Bethe_Bloch(ke, mass);
      double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
      pdf->SetParameters(par);
      if(pdf->Eval(dEdx) == 0)continue;
      else
     	likelihood += log(pdf->Eval(dEdx));
    }
    L_v.push_back(L);
    Likelihood_v.push_back(likelihood);
    L += step;
  }

  // TGraph* ti = new TGraph(L_v.size(),&L_v[0],&Likelihood_v[0]);
  // ti->Draw("al*");gPad->Update();gPad->WaitPrimitive();
  
  auto it = std::max_element(Likelihood_v.begin(), Likelihood_v.end());
  int index = std::distance(Likelihood_v.begin(), it);

  return Likelihood_v[index];
}


void likelihood_microtree(){

  gStyle->SetOptFit(1);
  
  //get minitree
  TFile* file = TFile::Open("/home/migue/Documents/PhD/tesis/files/microtrees/xs/xs_03-05-2023.root");
  TTree* tree = (TTree*)file->Get("ana");

  //get kaon pid
  TFile* kfile = TFile::Open("/home/migue/Documents/PhD/tesis/macros/analysis/xs/dedx_profiles.root");
  TFile* rfile = TFile::Open("/home/migue/Documents/PhD/tesis/macros/analysis/xs/ke_vs_range.root");

  TCanvas* c1 = new TCanvas("c1","c1",700,500);
  TCanvas* c2 = new TCanvas("c2","c2",700,500);
  
  TH1F* hk = new TH1F("hkaon","hkaon",20,-5,5);
  TH1F* hp = new TH1F("hproton","hproton",20,-5,5);
  TH1F* ho = new TH1F("hother","hother",20,-5,5);
  hk->SetLineColor(2);
  hk->SetLineWidth(3);
  ho->SetLineColor(9);
  ho->SetLineWidth(3);
  hp->SetLineWidth(3);
  
  int counter = 0;
  int counter_k = 0;
  int counter_p = 0;
  
  int nentries = tree->GetEntries();
  //loop over entries;
  for(int ientry = 0; ientry < nentries; ientry++){
    if(ientry%10==0)std::cout << ientry << "/" << nentries << std::endl;

    c1->cd();
    tree->GetEntry(ientry);

    //check if kaon or proton
    bool kaon = false;
    bool proton = false;
    int nkaon = tree->Draw("bestcandidate_truepdg","bestcandidate_truepdg==321","goff",1,ientry);
    int nproton = tree->Draw("bestcandidate_truepdg","bestcandidate_truepdg==2212","goff",1,ientry);
    //if(nkaon == 0 && nproton == 0)continue;
    //else if(nkaon > 0 && nproton == 0)kaon = true;
    if(nkaon > 0 && nproton == 0)kaon = true;
    else if(nkaon == 0 && nproton>0)proton = true;

    int n = tree->Draw("bestcandidate_hit_dedx:bestcandidate_hit_resrange","bestcandidate_hit_dedx>1.5 && bestcandidate_hit_dedx<100 && bestcandidate_hit_resrange!=-999","",1,ientry);
    TGraph* tg = new TGraph(n,tree->GetV2(),tree->GetV1());
    tg->RemovePoint(0);
    tg->RemovePoint(tg->GetN()-1);

    //kaon computation
    tg_rr_ke = (TGraph*)rfile->Get("kaon");
    tg_dedx_rr = (TGraph*)kfile->Get("kaon");
    mass = 493.677;
    double lk_kaon = LikelihoodFit(tg);

    //proton computation
    tg_rr_ke = (TGraph*)rfile->Get("proton");
    tg_dedx_rr = (TGraph*)kfile->Get("proton");
    mass = 938.272;
    double lk_proton = LikelihoodFit(tg);
    
    double ratio = -2*log(lk_kaon/lk_proton);
    
    if(kaon){
      hk->Fill(ratio);
      if(ratio>0){
	counter++;
	counter_k++;
      }
      
    }
    else if(proton){
      hp->Fill(ratio);
      if(ratio>0){
	counter++;
	counter_p++;
      }
    }
    else if(!kaon && !proton){
      ho->Fill(ratio);
      if(ratio>0)
	counter++;
    }

    c2->cd();
    hp->Draw();
    hk->Draw("same");
    ho->Draw("same");
    gPad->Update();
  }

  
  c2->cd();
  hp->Draw();
  hk->Draw("same");
  ho->Draw("same");

  std::cout << counter << " " << counter_k << " " << counter_p << std::endl;
}
    
