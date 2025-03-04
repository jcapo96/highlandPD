TH1F* hm;

void Minimization(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  double lambda1 = par[0];
  double lambda2 = par[1];
  double lambda3 = par[2];
  double alpha   = par[3];
  double beta    = par[4];
  double gamma   = par[5];

  double lambda_av = alpha*lambda1+beta*lambda2+gamma*lambda3;
  double B = hm->Integral("width")*lambda_av/(exp(lambda_av*200)-exp(lambda_av*10));
  std::cout << B << std::endl;

  double chi = 0;
  for(int ibin = 0; ibin < hm->GetNbinsX(); ibin++){
    double x = hm->GetBinCenter(ibin+1);
    double hvalue = hm->GetBinContent(ibin+1);
    double fvalue = B*exp(x*(alpha*lambda1+beta*lambda2+gamma*lambda3));
    chi += pow(hvalue-fvalue,2)/fvalue;
  }

  double reg = abs(1-alpha-beta-gamma);

  f=chi+reg;
}


Double_t ExpoCombination(double *x, double *par){

  double f = 0;
  double expo = par[0]*exp(par[4]*par[1]*x[0]+par[5]*par[2]*x[0]+par[6]*par[3]*x[0]);
  double reg = abs(1-par[4]-par[5]-par[6]);

  f = expo+reg;
  return f;
}

void unknown_estimation(){
  
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

  TF1* fexpo = new TF1("fexpo","[0]*exp([1]*x)",10,200);
  fexpo->SetParameters(0.010,-0.008);
  
  //pions fit
  d->Draw("seltrk_length>>h_pion(190,10,200)","beam_pdgs==211");
  TH1F* h_pion = (TH1F*)gPad->GetPrimitive("h_pion");
  h_pion->Sumw2();
  h_pion->Scale(1/h_pion->Integral());
  h_pion->Fit("fexpo");
  double c_pion = fexpo->GetParameter(0);
  double s_pion = fexpo->GetParameter(1);
  
  d->Draw("seltrk_length>>h_proton(190,10,200)","beam_pdgs==2212");
  TH1F* h_proton = (TH1F*)gPad->GetPrimitive("h_proton");
  h_proton->Sumw2();
  h_proton->Scale(1/h_proton->Integral());
  h_proton->Fit("fexpo");
  double c_proton = fexpo->GetParameter(0);
  double s_proton = fexpo->GetParameter(1);

  d->Draw("seltrk_length>>h_kaon(190,10,200)","beam_pdgs==321");
  TH1F* h_kaon = (TH1F*)gPad->GetPrimitive("h_kaon");
  h_kaon->Sumw2();
  h_kaon->Scale(1/h_kaon->Integral());
  h_kaon->Fit("fexpo");
  double c_kaon = fexpo->GetParameter(0);
  double s_kaon = fexpo->GetParameter(1);

  h_pion->Draw();
  h_proton->Draw("same");
  h_kaon->Draw("same");
  
  //pions fit
  d->Draw("seltrk_length>>h_u(190,10,200)","beam_npdgs==0");
  TH1F* h_u = (TH1F*)gPad->GetPrimitive("h_u");
  //h_u->Scale(1/h_u->Integral());

  // TF1* f = new TF1("f",ExpoCombination,10,200,7);
  // f->FixParameter(1,s_pion);
  // f->FixParameter(2,s_proton);
  // f->FixParameter(3,s_kaon);
  // f->SetParLimits(4,0.8,0.84);
  // f->SetParLimits(5,0.05,0.14);
  // f->SetParLimits(6,0.02,0.1);
  // f->SetParameter(4,0.83);
  // f->SetParameter(5,0.11);
  // f->SetParameter(6,0.05);
  // h_u->Fit("f");


  TMinuit* t = new TMinuit(6);
  t->SetFCN(Minimization);

  hm = (TH1F*)h_u->Clone("hm");
  std::cout << hm->Integral() << std::endl;
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  // arglist[0] = 1;
  // t->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[6] = {s_pion,s_proton,s_kaon,0.83,0.12,0.05};
  Double_t step[6]   = {0.01  , 0.01   , 0.01 ,0.0001,0.0001,0.0001};
    
  t->mnparm(0,  "lambda1",  vstart[0],  step[0],  0, 0, ierflg);
  t->mnparm(1,  "lambda2",  vstart[1],  step[1],  0, 0, ierflg);
  t->mnparm(2,  "lambda3",  vstart[2],  step[2],  0, 0, ierflg);
  t->mnparm(3,  "alpha"  ,  vstart[3],  step[3],  0.79, 0.89, ierflg);
  t->mnparm(4,  "beta"   ,  vstart[4],  step[4],  0.08, 0.16, ierflg);
  t->mnparm(5,  "gamma"  ,  vstart[5],  step[5],  0.02, 0.09, ierflg);
  t->FixParameter(0);
  t->FixParameter(1);
  t->FixParameter(2);
  //t->FixParameter(5);

  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  t->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  t->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  t->mnprin(3,amin);

  double A,alpha,beta,gamma,dummy;
  t->GetParameter(3,alpha,dummy);
  t->GetParameter(4,beta,dummy);
  t->GetParameter(5,gamma,dummy);

  double lambda_av = s_pion*alpha+s_proton*beta+s_kaon*gamma;
  A = hm->Integral("width")*lambda_av/(exp(lambda_av*200)-exp(lambda_av*10));

  TF1* f = new TF1("f",ExpoCombination,10,200,7);
  f->SetParameters(A,s_pion,s_proton,s_kaon,alpha,beta,gamma);
  f->Draw("same");
  
}
