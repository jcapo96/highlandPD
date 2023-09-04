//********************************************
Double_t Langaus(Double_t *x, Double_t *par) {
//********************************************
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  //Double_t np = 100.0;      // number of convolution steps
  Double_t np = 1000.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

//********************************************
Double_t DoubleLangaus(Double_t *x, Double_t *par) {
//********************************************

  double slw    = par[0];
  double smpv   = par[1];
  double snorm  = par[2];
  double sgw    = par[3];
  double blw    = par[4];
  double bmpv   = par[5];
  double bnorm  = par[6];
  double bgw    = par[7];

  double spar[] = {slw,smpv,snorm,sgw};
  double bpar[] = {blw,bmpv,bnorm,bgw};

  float S = Langaus(x,spar);
  float B = Langaus(x,bpar);

  return S+B;
}

//**********************************************
void SetStyle(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();

  //MC histo style
  sb_sample->UseCurrentStyle();
  sb_sample->SetHistogramLineWidth(3);
  sb_sample->SetHistogramColor(1);
  for(int i = 0; i < sb_sample->GetSize(); i++){
    sb_sample->GetHistVector()[i]->SetTitle("");
    sb_sample->GetHistVector()[i]->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    sb_sample->GetHistVector()[i]->GetYaxis()->SetTitle("Probability density");
    sb_sample->GetHistVector()[i]->GetXaxis()->CenterTitle();
    sb_sample->GetHistVector()[i]->GetYaxis()->CenterTitle();
    sb_sample->GetHistVector()[i]->GetYaxis()->SetTitleSize(20);
    sb_sample->GetHistVector()[i]->GetYaxis()->SetTitleFont(43);
    sb_sample->GetHistVector()[i]->GetYaxis()->SetLabelFont(43); 
    sb_sample->GetHistVector()[i]->GetYaxis()->SetLabelSize(20);
    sb_sample->GetHistVector()[i]->GetYaxis()->SetTitleOffset(1.55);
    sb_sample->GetHistVector()[i]->GetXaxis()->SetTitleSize(20);
    sb_sample->GetHistVector()[i]->GetXaxis()->SetTitleFont(43);
    sb_sample->GetHistVector()[i]->GetXaxis()->SetLabelFont(43); 
    sb_sample->GetHistVector()[i]->GetXaxis()->SetLabelSize(20);
  }

  //MC fit style
  sb_sample->GetSignal()->SetCFitStyle(2);
  sb_sample->GetBackground()->SetCFitStyle(2);
  sb_sample->GetBackground()->GetSemiBackground()->SetCFitStyle(2);
  sb_sample->GetSignal()->SetCFitColor(1);
  sb_sample->GetBackground()->SetCFitColor(1);
  sb_sample->GetBackground()->GetSemiBackground()->SetCFitColor(1);

  //data histo style
  sbd_sample->UseCurrentStyle();
  sbd_sample->SetHistogramLineWidth(3);
  sbd_sample->SetHistogramColor(1);
  for(int i = 0; i < sb_sample->GetSize(); i++){
    sbd_sample->GetHistVector()[i]->SetTitle("");
    sbd_sample->GetHistVector()[i]->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetTitle("Probability density");
    sbd_sample->GetHistVector()[i]->GetXaxis()->CenterTitle();
    sbd_sample->GetHistVector()[i]->GetYaxis()->CenterTitle();
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetTitleSize(20);
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetTitleFont(43);
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetLabelFont(43); 
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetLabelSize(20);
    sbd_sample->GetHistVector()[i]->GetYaxis()->SetTitleOffset(1.55);
    sbd_sample->GetHistVector()[i]->GetXaxis()->SetTitleSize(20);
    sbd_sample->GetHistVector()[i]->GetXaxis()->SetTitleFont(43);
    sbd_sample->GetHistVector()[i]->GetXaxis()->SetLabelFont(43); 
    sbd_sample->GetHistVector()[i]->GetXaxis()->SetLabelSize(20);
  }

  //data fit style
  sbd_sample->GetSignal()->SetCFitStyle(2);
  sbd_sample->GetBackground()->SetCFitStyle(2);
  sbd_sample->GetBackground()->GetSemiBackground()->SetCFitStyle(2);
  sbd_sample->GetSignal()->SetCFitColor(1);
  sbd_sample->GetBackground()->SetCFitColor(1);
  sbd_sample->GetBackground()->GetSemiBackground()->SetCFitColor(1);
  
}
  
//**********************************************
void GetHistoPlots(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();

  //canvas for drawing
  TCanvas* c = new TCanvas("c","c",700,500);
  
  TLatex tt1;
  tt1.SetNDC();
  TLatex tt2;
  tt2.SetNDC();
  //tt2.SetTextSize(0.04);
  tt2.SetTextAlign(31);
  
  //print RR slices
  gStyle->SetOptStat(0);
  for(int i = 0; i < sb_sample->GetSize(); i++){
    //Draw Histograms
    sb_sample->GetHistVector()[i]->Draw("");
    gPad->RedrawAxis();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"MC");

    // //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/mc_slice_"+ssh.str()+".pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
  }

  for(int i = 0; i < sbd_sample->GetSize(); i++){
    //Draw Histograms
    sbd_sample->GetHistVector()[i]->Draw("");
    gPad->RedrawAxis();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"Data");

    // //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/d_slice_"+ssh.str()+".pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
  }

  delete c;

  //canvas for drawing
  c = new TCanvas("c","c",700,700);
  
  //Print RR slices comparing data and MC
  //histo style
  for(int i = 0; i < sb_sample->GetSize(); i++){
    //look for higher minimum
    double min1 = sb_sample->GetHistVector()[i]->GetBinCenter(1)-sb_sample->GetHistVector()[i]->GetBinWidth(1)/2;
    double min2 = sbd_sample->GetHistVector()[i]->GetBinCenter(1)-sbd_sample->GetHistVector()[i]->GetBinWidth(1)/2;
    double min = (min1<min2) ? min2 : min1;

    //look for lower minimum
    int nbins1 = sb_sample->GetHistVector()[i]->GetNbinsX();
    double max1 = sb_sample->GetHistVector()[i]->GetBinCenter(nbins1)+sb_sample->GetHistVector()[i]->GetBinWidth(nbins1)/2;
    int nbins2 = sbd_sample->GetHistVector()[i]->GetNbinsX();
    double max2 = sbd_sample->GetHistVector()[i]->GetBinCenter(nbins2)+sbd_sample->GetHistVector()[i]->GetBinWidth(nbins2)/2;
    double max = (max1<max2) ? max1 : max2;

    //get new array of bins edges
    std::vector<double> edges;
    edges.push_back(min);
    double width = sb_sample->GetHistVector()[i]->GetBinWidth(1);
    double lastedge = min+width;
    int nbins = 0;
    while(lastedge<=max){
      edges.push_back(lastedge);
      lastedge += width;
      nbins++;
    }

    //rebin both histograms
    TH1F* h1 = (TH1F*)sb_sample->GetHistVector()[i]->Rebin(nbins,"MC",&edges[0]);
    TH1F* h2 = (TH1F*)sbd_sample->GetHistVector()[i]->Rebin(nbins,"Data",&edges[0]);

    //draw ratio
    c->Clear();
    c->cd();
    TRatioPlot* r = new TRatioPlot(h1,h2,"diffsig");
    r->Draw();

    //ratio style
    r->GetLowerRefYaxis()->SetTitle("(MC-Data)/error_{MC}");
    r->GetLowerRefYaxis()->CenterTitle();
    r->GetLowerRefXaxis()->CenterTitle();
    // r->GetLowerRefYaxis()->SetRangeUser(-3,3);
    r->GetLowerRefYaxis()->SetTitleSize(20);
    r->GetLowerRefYaxis()->SetTitleFont(43);
    r->GetLowerRefYaxis()->SetLabelFont(43); 
    r->GetLowerRefYaxis()->SetLabelSize(20);
    r->GetLowerRefYaxis()->SetTitleOffset(1.55);
    r->GetLowerRefXaxis()->SetTitleSize(20);
    r->GetLowerRefXaxis()->SetTitleFont(43);
    r->GetLowerRefXaxis()->SetLabelFont(43);
    r->GetLowerRefXaxis()->SetLabelSize(20);
    r->GetLowerRefXaxis()->SetTitleOffset(1.55);

    r->GetUpperPad()->cd();

    //draw the three functions
    gPad->RedrawAxis();
    
    //write ProtoDUNE-SP
    tt1.DrawLatex(0.10,0.91,"#bf{DUNE:ProtoDUNE-SP}");
    //Write RR in the plot
    std::stringstream srmin, srmax;
    double rmin = sb_sample->GetRRVector()[i].first - 1;
    double rmax = sb_sample->GetRRVector()[i].first + 1;
    srmin << rmin;
    srmax << rmax;
    tt2.DrawLatex(0.9,0.91,(""+srmin.str()+" < RR [cm] < "+srmax.str()+"").c_str());

    //legend
    TLegend* lg = new TLegend(0.7,0.7,0.85,0.85);
    lg->AddEntry(h1,"MC","l");
    lg->AddEntry(h2,"Data","lp");
    lg->Draw("same");
    
    gPad->Update();//gPad->WaitPrimitive();

    //print plots
    std::stringstream ssh;
    ssh << i;
    c->Print(("plots/mc_d_slice_"+ssh.str()+".pdf").c_str());
  }

  delete c;
}

//**********************************************
void GetHistoPlotsWithIncoherentFits(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();

  //canvas for drawing
  TCanvas* c = new TCanvas("c","c",700,500);

  //style for plot results
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.4);
  gStyle->SetOptFit(001);
  gStyle->SetFitFormat("3.2f");
  gStyle->SetOptStat(0);

  //texts
  TLatex tt1;
  tt1.SetNDC();
  TLatex tt2;
  tt2.SetNDC();
  //tt2.SetTextSize(0.04);
  tt2.SetTextAlign(31);

  //function to be fitted
  TF1* f = new TF1("f",DoubleLangaus,0,30,8);
  f->SetParameters(0.1,1.7,0.5,0.05,0.1,8,0.5,0.1);
  f->SetParLimits(4,0.01,5);
  f->SetParNames("#sigma_{L,B}","#mu_{B}","N_{B}","#sigma_{G,B}","#sigma_{L,S}","#mu_{S}","N_{S}","#sigma_{G,S}");

  TH1F* dummy;
  //MC
  for(int i = 0; i < sb_sample->GetSize(); i++){
    //Draw Histograms
    dummy = (TH1F*)sb_sample->GetHistVector()[i]->Clone();
    dummy->Draw();
    dummy->Fit("f","L");
    // sb_sample->GetHistVector()[i]->Draw("");
    // sb_sample->GetHistVector()[i]->Fit("f","L");
    gPad->RedrawAxis();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"MC");

    // //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/CF_mc_slice_"+ssh.str()+"_IF.pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
    delete dummy;
  }

  f->SetParameters(0.1,1.7,0.5,0.05,0.1,8,0.5,0.1);
  for(int i = 0; i < sbd_sample->GetSize(); i++){
    //Draw Histograms
    dummy = (TH1F*)sbd_sample->GetHistVector()[i]->Clone();
    dummy->Draw();
    dummy->Fit("f","L");
    // sbd_sample->GetHistVector()[i]->Draw("");
    // sbd_sample->GetHistVector()[i]->Fit("f","L");
    gPad->RedrawAxis();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"Data");

    // //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/CF_d_slice_"+ssh.str()+"_IF.pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
    delete dummy;
  }

  delete c;
  
}

//**********************************************
void GetMCTrueHistoPlotsWithIncoherentFits(CoherentFit* cf_mc){
//**********************************************

  CoherentSample* ts = cf_mc->GetTrueSignalSample();
  CoherentSample* tb = cf_mc->GetTrueBackgroundSample();
  
  TLatex tt1;
  tt1.SetNDC();
  TLatex tt2;
  tt2.SetNDC();
  //tt2.SetTextSize(0.04);
  tt2.SetTextAlign(31);
  
  //style
  gStyle->SetFitFormat("3.4f");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  ts->UseCurrentStyle();
  ts->SetHistogramLineWidth(3);
  tb->UseCurrentStyle();
  tb->SetHistogramLineWidth(3);

  //signal
  for(int i = 0; i < ts->GetSize(); i++){
    //Draw Histograms
    ts->GetHistVector()[i]->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    ts->GetHistVector()[i]->GetYaxis()->SetTitle("Probability Density");
    ts->GetHistVector()[i]->GetXaxis()->CenterTitle();
    ts->GetHistVector()[i]->GetYaxis()->CenterTitle();
    ts->GetHistVector()[i]->SetTitle(""); //we don't want title displayed in the top of the canvas
    ts->GetHistVector()[i]->Draw("");
    ts->GetIFitVector()[i]->Draw("same");

    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

    //Write RR in the plot
    std::stringstream srmin, srmax;
    double rmin = ts->GetRRVector()[i].first - ts->GetRRVector()[i].second;
    double rmax = ts->GetRRVector()[i].first + ts->GetRRVector()[i].second;
    srmin << rmin;
    srmax << rmax;
    tt2.DrawLatex(0.9,0.94,(""+srmin.str()+" < RR [cm] < "+srmax.str()+"").c_str());

    gPad->RedrawAxis();
    	  
    //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/CF_IF_mc_ts_slice_"+ssh.str()+".pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
  }

  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.3);
  //background
  for(int i = 0; i < tb->GetSize(); i++){
    //Draw Histograms
    tb->GetHistVector()[i]->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    tb->GetHistVector()[i]->GetYaxis()->SetTitle("Probability Density");
    tb->GetHistVector()[i]->GetXaxis()->CenterTitle();
    tb->GetHistVector()[i]->GetYaxis()->CenterTitle();
    tb->GetHistVector()[i]->SetTitle(""); //we don't want title displayed in the top of the canvas
    tb->GetHistVector()[i]->Draw("");
    tb->GetIFitVector()[i]->Draw("same");

    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

    //Write RR in the plot
    std::stringstream srmin, srmax;
    double rmin = tb->GetRRVector()[i].first - tb->GetRRVector()[i].second;
    double rmax = tb->GetRRVector()[i].first + tb->GetRRVector()[i].second;
    srmin << rmin;
    srmax << rmax;
    tt2.DrawLatex(0.9,0.94,(""+srmin.str()+" < RR [cm] < "+srmax.str()+"").c_str());

    gPad->RedrawAxis();
    gPad->Update();
    	  
    //Print plots
    std::stringstream ssh;
    ssh << i;
    gPad->Print(("plots/CF_IF_mc_tb_slice_"+ssh.str()+".pdf").c_str());
    gPad->Update();//gPad->WaitPrimitive();
  }
  
}

//**********************************************
void GetParametrizationPlots(CoherentFit* cf_mc){
//**********************************************

  CoherentSample* ts = cf_mc->GetTrueSignalSample();
  CoherentSample* tb = cf_mc->GetTrueBackgroundSample();
  
  TLatex tt1;
  tt1.SetNDC();
   
  //style
  gStyle->SetFitFormat("3.4f");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  ts->UseCurrentStyle();
  ts->SetHistogramLineWidth(3);
  tb->UseCurrentStyle();
  tb->SetHistogramLineWidth(3);
 
  //signal
  TGraphErrors* tg_ts_lw = ts->GetIlwGraph();
  TGraphErrors* tg_ts_mpv = ts->GetImpvGraph();
  TGraphErrors* tg_ts_norm = ts->GetInormGraph();
  TGraphErrors* tg_ts_gw = ts->GetIgwGraph();

  tg_ts_lw->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_ts_lw->GetYaxis()->SetTitle("#sigma_{L,S} [MeV/cm]");
  tg_ts_lw->GetXaxis()->CenterTitle();
  tg_ts_lw->GetYaxis()->CenterTitle();
  tg_ts_lw->SetTitle("");
  tg_ts_lw->Draw("ap");
  TF1* f_ts_lw = new TF1("f_ts_lw","([0]/x-1)/x+[1]",1,61);
  f_ts_lw->SetParNames("#alpha","#beta");
  f_ts_lw->SetParameters(4.2845,0.0947);
  tg_ts_lw->Fit("f_ts_lw","");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_ts_lw.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  tg_ts_mpv->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_ts_mpv->GetYaxis()->SetTitle("#mu_{S} [MeV/cm]");
  tg_ts_mpv->GetXaxis()->CenterTitle();
  tg_ts_mpv->GetYaxis()->CenterTitle();
  tg_ts_mpv->SetTitle("");
  tg_ts_mpv->Draw("ap");
  TF1* f_ts_mpv = new TF1("f_ts_mpv","[0]*(([1]*x-1)/([1]*x+1))+[2]",1,61);
  f_ts_mpv->SetParNames("#alpha","#beta","#gamma");
  tg_ts_mpv->Fit("f_ts_mpv","");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_ts_mpv_nofunction.pdf");
  gPad->Update();//gPad->WaitPrimitive();
  
  tg_ts_norm->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_ts_norm->GetYaxis()->SetTitle("N_{S}");
  tg_ts_norm->GetXaxis()->CenterTitle();
  tg_ts_norm->GetYaxis()->CenterTitle();
  tg_ts_norm->SetTitle("");
  tg_ts_norm->Draw("ap");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_ts_norm.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  tg_ts_gw->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_ts_gw->GetYaxis()->SetTitle("#sigma_{G,S} [MeV/cm]");
  tg_ts_gw->GetXaxis()->CenterTitle();
  tg_ts_gw->GetYaxis()->CenterTitle();
  tg_ts_gw->SetTitle("");
  tg_ts_gw->Draw("ap");
  TF1* f_ts_gw = new TF1("f_ts_gw","[0]/(x+1)+[1]",1,61);
  f_ts_gw->SetParNames("#alpha","#beta");
  tg_ts_gw->Fit("f_ts_gw","");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_ts_gw.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  //background
  TGraphErrors* tg_tb_lw = tb->GetIlwGraph();
  TGraphErrors* tg_tb_mpv = tb->GetImpvGraph();
  TGraphErrors* tg_tb_norm = tb->GetInormGraph();
  TGraphErrors* tg_tb_gw = tb->GetIgwGraph();
  TGraphErrors* tg_tb_const = tb->GetIconstGraph();

  tg_tb_lw->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_tb_lw->GetYaxis()->SetTitle("#sigma_{L,B} [MeV/cm]");
  tg_tb_lw->GetXaxis()->CenterTitle();
  tg_tb_lw->GetYaxis()->CenterTitle();
  tg_tb_lw->SetTitle("");
  tg_tb_lw->Draw("ap");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_tb_lw.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  tg_tb_mpv->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_tb_mpv->GetYaxis()->SetTitle("#mu_{B} [MeV/cm]");
  tg_tb_mpv->GetXaxis()->CenterTitle();
  tg_tb_mpv->GetYaxis()->CenterTitle();
  tg_tb_mpv->SetTitle("");
  tg_tb_mpv->Draw("ap");
  TF1* f_tb_mpv = new TF1("f_tb_mpv","[0]+[1]*x",1,61);
  f_tb_mpv->SetParNames("#alpha","#beta");
  tg_tb_mpv->Fit("f_tb_mpv","");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_tb_mpv.pdf");
  gPad->Update();//gPad->WaitPrimitive();
  
  tg_tb_norm->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_tb_norm->GetYaxis()->SetTitle("N_{B}");
  tg_tb_norm->GetXaxis()->CenterTitle();
  tg_tb_norm->GetYaxis()->CenterTitle();
  tg_tb_norm->SetTitle("");
  tg_tb_norm->Draw("ap");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_tb_norm.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  tg_tb_gw->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_tb_gw->GetYaxis()->SetTitle("#sigma_{G,B} [MeV/cm]");
  tg_tb_gw->GetXaxis()->CenterTitle();
  tg_tb_gw->GetYaxis()->CenterTitle();
  tg_tb_gw->SetTitle("");
  tg_tb_gw->Draw("ap");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_tb_gw.pdf");
  gPad->Update();//gPad->WaitPrimitive();

  tg_tb_const->GetXaxis()->SetTitle("Residual Range [cm]");
  tg_tb_const->GetYaxis()->SetTitle("C_{B} [MeV/cm]");
  tg_tb_const->GetXaxis()->CenterTitle();
  tg_tb_const->GetYaxis()->CenterTitle();
  tg_tb_const->SetTitle("");
  tg_tb_const->Draw("ap");
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  gPad->Print("plots/CF_IF_mc_tb_const.pdf");
  gPad->Update();//gPad->WaitPrimitive();

}

//**********************************************
void GetCoherentFitPlots(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();
  
  TH1F* h_pulls_mc = new TH1F("","",10,-5,5);
  TH1F* h_pulls_d = new TH1F("","",10,-5,5);
  TH1F* h_pulls_mc_total = new TH1F("","",25,-5,5);
  TH1F* h_pulls_d_total = new TH1F("","",25,-5,5);

  h_pulls_mc->GetXaxis()->SetTitle("Pulls");
  h_pulls_mc->GetXaxis()->CenterTitle();
  h_pulls_d->GetXaxis()->SetTitle("Pulls");
  h_pulls_d->GetXaxis()->CenterTitle();
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  
  //print slices
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //canvas for drawing
  TCanvas* c = new TCanvas("c","c",700,700);

  TLatex tt1;
  tt1.SetNDC();

  TLatex tt2;
  tt2.SetNDC();
  //tt2.SetTextSize(0.04);
  tt2.SetTextAlign(31);
  
  for(int i = 0; i < sb_sample->GetSize(); i++){
    c->cd();
    TRatioPlot* r = sb_sample->GetFitRatio(i);
    r->Draw();
    r->GetLowerRefYaxis()->SetTitle("Pulls");
    r->GetLowerRefYaxis()->CenterTitle();
    r->GetLowerRefXaxis()->CenterTitle();
    r->GetLowerRefYaxis()->SetRangeUser(-3,3);
    r->GetLowerRefYaxis()->SetTitleSize(20);
    r->GetLowerRefYaxis()->SetTitleFont(43);
    r->GetLowerRefYaxis()->SetLabelFont(43); 
    r->GetLowerRefYaxis()->SetLabelSize(20);
    r->GetLowerRefYaxis()->SetTitleOffset(1.55);
    r->GetLowerRefXaxis()->SetTitleSize(20);
    r->GetLowerRefXaxis()->SetTitleFont(43);
    r->GetLowerRefXaxis()->SetLabelFont(43);
    r->GetLowerRefXaxis()->SetLabelSize(20);
    r->GetLowerRefXaxis()->SetTitleOffset(1.55);

    r->GetUpperPad()->cd();

    //draw the three functions
    sb_sample->GetSignal()->GetCFitVector()[i]->Draw("same");
    sb_sample->GetBackground()->GetCFitVector()[i]->Draw("same");
    sb_sample->GetBackground()->GetSemiBackground()->GetCFitVector()[i]->Draw("same");

    //draw the three functions
    gPad->RedrawAxis();

    tt1.DrawLatex(0.10,0.91,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.91,"MC");
  
    gPad->Update();//gPad->WaitPrimitive();
  
    //print plots
    std::stringstream ssh;
    ssh << i;
    c->Print(("plots/CF_mc_slice_"+ssh.str()+".pdf").c_str());

    //fill pulls histogram
    c->cd();
    TGraph* pulls = r->GetCalculationOutputGraph();
    for(int j = 0; j < pulls->GetN(); j++){
      h_pulls_mc->Fill(pulls->GetPointY(j));
      h_pulls_mc_total->Fill(pulls->GetPointY(j));
    }
    gStyle->SetOptStat("emr");
    h_pulls_mc->Draw();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"MC");
    gPad->Update();//gPad->WaitPrimitive();
    c->Print(("plots//CF_mc_slice_"+ssh.str()+"_pulls.pdf").c_str());
    h_pulls_mc->Reset();
    gStyle->SetOptStat(0);
  }

  for(int i = 0; i < sbd_sample->GetSize(); i++){
    c->cd();
    TRatioPlot* r = sbd_sample->GetFitRatio(i);
    r->Draw();
    r->GetLowerRefYaxis()->SetTitle("Pulls");
    r->GetLowerRefYaxis()->CenterTitle();
    r->GetLowerRefXaxis()->CenterTitle();
    r->GetLowerRefYaxis()->SetRangeUser(-3,3);
    r->GetLowerRefYaxis()->SetTitleSize(20);
    r->GetLowerRefYaxis()->SetTitleFont(43);
    r->GetLowerRefYaxis()->SetLabelFont(43); 
    r->GetLowerRefYaxis()->SetLabelSize(20);
    r->GetLowerRefYaxis()->SetTitleOffset(1.55);
    r->GetLowerRefXaxis()->SetTitleSize(20);
    r->GetLowerRefXaxis()->SetTitleFont(43);
    r->GetLowerRefXaxis()->SetLabelFont(43);
    r->GetLowerRefXaxis()->SetLabelSize(20);
    r->GetLowerRefXaxis()->SetTitleOffset(1.55);

    r->GetUpperPad()->cd();

    //draw the three functions
    sbd_sample->GetSignal()->GetCFitVector()[i]->Draw("same");
    sbd_sample->GetBackground()->GetCFitVector()[i]->Draw("same");
    sbd_sample->GetBackground()->GetSemiBackground()->GetCFitVector()[i]->Draw("same");
  
    //draw the three functions
    gPad->RedrawAxis(); gPad->Update();

    tt1.DrawLatex(0.10,0.91,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.91,"Data");
  
    gPad->Update();//gPad->WaitPrimitive();
  
    //print plots
    std::stringstream ssh;
    ssh << i;
    c->Print(("plots//CF_d_slice_"+ssh.str()+".pdf").c_str());

    //fill pulls histogram
    c->cd();
    TGraph* pulls = r->GetCalculationOutputGraph();
    for(int j = 0; j < pulls->GetN(); j++){
      h_pulls_d->Fill(pulls->GetPointY(j));
      h_pulls_d_total->Fill(pulls->GetPointY(j));
    }

    gStyle->SetOptStat("emr");
    h_pulls_d->Draw();
    tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
    tt2.DrawLatex(0.9,0.94,"Data");
    gPad->Update();//gPad->WaitPrimitive();
    c->Print(("plots//CF_d_slice_"+ssh.str()+"_pulls.pdf").c_str());
    h_pulls_d->Reset();
    gStyle->SetOptStat(0);
  }

  gStyle->SetOptStat("emr");
  TCanvas* c2 = new TCanvas("c2","c2",700,500);
  c2->cd();
  h_pulls_mc_total->GetXaxis()->SetTitle("Pulls");
  h_pulls_mc_total->GetXaxis()->CenterTitle();
  h_pulls_mc_total->Draw();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.94,"MC");
  gPad->Update();//gPad->WaitPrimitive();
  c->Print("plots//CF_mc_pulls_total.pdf");

  h_pulls_d_total->GetXaxis()->SetTitle("Pulls");
  h_pulls_d_total->GetXaxis()->CenterTitle();
  h_pulls_d_total->Draw();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  tt2.DrawLatex(0.9,0.94,"Data");
  gPad->Update();//gPad->WaitPrimitive();
  c->Print("plots//CF_d_pulls_total.pdf");

  delete c;
  delete c2;
}

//**********************************************
void GetCorrelatedRandomNumbers3D(double (&crandom)[3], TRandom3* r, double (&cholesk)[3][3]){
//**********************************************
  
  double random[3] = {0};
  for(int j = 0; j < 3; j++)random[j] = r->Gaus(0,1);
  for(int j = 0; j < 3; j++){
    double product = 0;
    for(int k = 0; k < 3; k++)product += cholesk[j][k]*random[k];
    crandom[j] = product;
  }
}

//**********************************************
void GetCholeskyMatrix3D(double (&cholesk)[3][3], double c12, double c13, double c23){
//**********************************************
  
  double cmatrix[3][3] = {{1,c12,c13},{c12,1.000,c23},{c13,c23,1.000}};

  for(int i = 0; i < 3; i++){
    for(int j = 0; j <= i; j++){
      int sum = 0;
      
      if(j == i){
	for (int k = 0; k < j; k++)
	  sum += pow(cholesk[j][k], 2);
	cholesk[j][j] = sqrt(cmatrix[j][j]-sum);
      }
      else {
	for (int k = 0; k < j; k++)
	  sum += (cholesk[i][k] * cholesk[j][k]);
	cholesk[i][j] = (cmatrix[i][j] - sum) / cholesk[j][j];
      }
    }
  }
}

//**********************************************
TGraphErrors* GetMPVGraph(double* params, double* params_error, double* correlations){
//**********************************************
  
  double cholesk[3][3] = {0};
  GetCholeskyMatrix3D(cholesk,correlations[0],correlations[1],correlations[2]);

  TF1* f = new TF1("f","[0]*(([1]*x-1)/([1]*x+1)+[2])",1,61);
  TH2F* h = new TH2F("h","h",3000,1,61,8000,2,10);
  
  TRandom3* r = new TRandom3();
  r->SetSeed(1);

  double crandom[3] = {0};
  for(int i = 0; i < 10000; i++){

    GetCorrelatedRandomNumbers3D(crandom,r,cholesk);

    for(int j = 0; j < 3; j++)
      f->SetParameter(j,params[j]+params_error[j]*crandom[j]);

    for(int j = 0; j < 3000; j++){
      h->Fill(1+60./3000*j,f->Eval(1+60./3000*j));
    }
  }

  TH1D* dummy;
  std::vector<double> x,y,e;
  for(int i = 0; i < 3000; i++){
    x.push_back(h->GetXaxis()->GetBinCenter(i+1));
    
    dummy = h->ProjectionY("dummy",i+1,i+2);
    dummy->Fit("gaus","Q");

    y.push_back(dummy->GetFunction("gaus")->GetParameter(1));
    e.push_back(dummy->GetFunction("gaus")->GetParameter(2));
  }

  return new TGraphErrors(x.size(),&x[0],&y[0],0,&e[0]);
}

//**********************************************
void GetCorrelatedRandomNumbers2D(double (&crandom)[2], TRandom3* r, double (&cholesk)[2][2]){
//**********************************************
  
  double random[3] = {0};
  for(int j = 0; j < 2; j++)random[j] = r->Gaus(0,1);
  for(int j = 0; j < 2; j++){
    double product = 0;
    for(int k = 0; k < 2; k++)product += cholesk[j][k]*random[k];
    crandom[j] = product;
  }
}

//**********************************************
void GetCholeskyMatrix2D(double (&cholesk)[2][2], double c12){
//**********************************************
  
  double cmatrix[2][2] = {{1,c12},{c12,1.000}};

  for(int i = 0; i < 2; i++){
    for(int j = 0; j <= i; j++){
      int sum = 0;
      
      if(j == i){
	for (int k = 0; k < j; k++)
	  sum += pow(cholesk[j][k], 2);
	cholesk[j][j] = sqrt(cmatrix[j][j]-sum);
      }
      else {
	for (int k = 0; k < j; k++)
	  sum += (cholesk[i][k] * cholesk[j][k]);
	cholesk[i][j] = (cmatrix[i][j] - sum) / cholesk[j][j];
      }
    }
  }
}

//**********************************************
TGraphErrors* GetGWGraph(double* params, double* params_error, double* correlations){
//**********************************************
  
  double cholesk[2][2] = {0};
  GetCholeskyMatrix2D(cholesk,correlations[1]);

  TF1* f = new TF1("f","[0]/(x+1)+[1]",1,61);
  TH2F* h = new TH2F("h","h",2000,1,61,5000,0,5);
  
  TRandom3* r = new TRandom3();
  r->SetSeed(1);

  double crandom[2] = {0};
  for(int i = 0; i < 20000; i++){

    GetCorrelatedRandomNumbers2D(crandom,r,cholesk);

    for(int j = 0; j < 2; j++)
      f->SetParameter(j,params[j]+params_error[j]*crandom[j]);

    for(int j = 0; j < 2000; j++){
      h->Fill(1+60./2000*j,f->Eval(1+60./2000*j));
    }
  }

  // h->Draw("colz");gPad->Update();gPad->WaitPrimitive();
  
  TH1D* dummy;
  std::vector<double> x,y,e;
  for(int i = 0; i < 2000; i++){
    x.push_back(h->GetXaxis()->GetBinCenter(i+1));
    
    dummy = h->ProjectionY("dummy",i+1,i+2);
    dummy->Fit("gaus","Q");
    // gPad->Update();gPad->WaitPrimitive();

    y.push_back(dummy->GetFunction("gaus")->GetParameter(1));
    e.push_back(dummy->GetFunction("gaus")->GetParameter(2));
  }

  return new TGraphErrors(x.size(),&x[0],&y[0],0,&e[0]);
}

//**********************************************
TGraphErrors* GetLWGraph(double* params, double* params_error, double* correlations){
//**********************************************
  
  double cholesk[2][2] = {0};
  GetCholeskyMatrix2D(cholesk,correlations[1]);

  TF1* f = new TF1("f","([0]/x-1)/x+[1]",1,61);
  TH2F* h = new TH2F("h","h",2000,1,61,5000,0,5);
  
  TRandom3* r = new TRandom3();
  r->SetSeed(1);

  double crandom[2] = {0};
  for(int i = 0; i < 20000; i++){

    GetCorrelatedRandomNumbers2D(crandom,r,cholesk);

    for(int j = 0; j < 2; j++)
      f->SetParameter(j,params[j]+params_error[j]*crandom[j]);

    for(int j = 0; j < 2000; j++){
      h->Fill(1+60./2000*j,f->Eval(1+60./2000*j));
    }
  }

  // h->Draw("colz");gPad->Update();gPad->WaitPrimitive();
  
  TH1D* dummy;
  std::vector<double> x,y,e;
  for(int i = 0; i < 2000; i++){
    x.push_back(h->GetXaxis()->GetBinCenter(i+1));
    
    dummy = h->ProjectionY("dummy",i+1,i+2);
    dummy->Fit("gaus","Q");
    // gPad->Update();gPad->WaitPrimitive();

    y.push_back(dummy->GetFunction("gaus")->GetParameter(1));
    e.push_back(dummy->GetFunction("gaus")->GetParameter(2));
  }

  return new TGraphErrors(x.size(),&x[0],&y[0],0,&e[0]);
}

//**********************************************
void GetLWPlots(CoherentSample* sb_sample, CoherentSample* sbd_sample){
//**********************************************

  TCanvas* c = new TCanvas("c","c",700,500);
  
  //data
  double d_params[2]       = {sbd_sample->GetSignal()->GetClwA().first,  sbd_sample->GetSignal()->GetClwB().first };
  double d_params_error[2] = {sbd_sample->GetSignal()->GetClwA().second, sbd_sample->GetSignal()->GetClwB().second};
  double d_correlations[2] = {1          , sbd_sample->GetCorrelationMatrixElement(0,1)};
  
  TGraphErrors* tgd = GetLWGraph(d_params,d_params_error,d_correlations);
  tgd->SetTitle("Data");
  tgd->SetFillColorAlpha(9,0.5);

  //MC
  double mc_params[2]       = {sb_sample->GetSignal()->GetClwA().first,  sb_sample->GetSignal()->GetClwB().first};
  double mc_params_error[2] = {sb_sample->GetSignal()->GetClwA().second, sb_sample->GetSignal()->GetClwB().second};
  double mc_correlations[2] = {1           , sb_sample->GetCorrelationMatrixElement(0,1)};
  
  TGraphErrors* tgmc = GetLWGraph(mc_params,mc_params_error,mc_correlations);
  tgmc->SetTitle("MC");
  tgmc->SetFillColorAlpha(98,0.5);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(tgmc,"3");
  mg->Add(tgd,"3");

  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#sigma_{S,L} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();

  mg->GetYaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetRangeUser(0.02,4.8);
  mg->Draw("ap");
  gPad->BuildLegend(0.7,0.7,0.85,0.85,"","f");
  gPad->SetLogy();

  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update(); //gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_lw_comparative.pdf");

  //plot incoherent-coherent for MC
  TGraphErrors* tg_ts_lw = sb_sample->GetTrueSignal()->GetIlwGraph();
  tg_ts_lw->SetMarkerStyle(20);
  mg->RecursiveRemove(tgd);
  mg->Add(tg_ts_lw,"p");
  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#sigma_{S,L} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(0.9);
  //mg->GetYaxis()->SetRangeUser(0.02,4.8);
  mg->Draw("ap");

  TLegend* lg = new TLegend(0.6,0.6,0.86,0.86);
  lg->AddEntry(tgmc,"MC Coherent Fit", "f");
  lg->AddEntry(tg_ts_lw,"MC Incoherent Fit", "lp");
  lg->Draw();

  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update();//gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_lw_comparative_coherent-incoherent.pdf");

  delete c;
}

//**********************************************
void GetMPVPlots(CoherentSample* sb_sample, CoherentSample* sbd_sample){
//**********************************************

  TCanvas* c = new TCanvas("c","c",700,500);
  
  //data
  double d_params[3]       = {sbd_sample->GetSignal()->GetCmpvA().first, sbd_sample->GetSignal()->GetCmpvB().first, sbd_sample->GetSignal()->GetCmpvC().first };
  double d_params_error[3] = {sbd_sample->GetSignal()->GetCmpvA().second,sbd_sample->GetSignal()->GetCmpvB().second,sbd_sample->GetSignal()->GetCmpvC().second};
  double d_correlations[3] = {sbd_sample->GetCorrelationMatrixElement(2,3),sbd_sample->GetCorrelationMatrixElement(2,4),sbd_sample->GetCorrelationMatrixElement(3,4)};
  
  TGraphErrors* tgd = GetMPVGraph(d_params,d_params_error,d_correlations);
  tgd->SetTitle("Data");
  tgd->SetFillColorAlpha(9,0.5);

  //MC
  double mc_params[3]       = {sb_sample->GetSignal()->GetCmpvA().first, sb_sample->GetSignal()->GetCmpvB().first, sb_sample->GetSignal()->GetCmpvC().first };
  double mc_params_error[3] = {sb_sample->GetSignal()->GetCmpvA().second,sb_sample->GetSignal()->GetCmpvB().second,sb_sample->GetSignal()->GetCmpvC().second};
  double mc_correlations[3] = {sb_sample->GetCorrelationMatrixElement(2,3),sb_sample->GetCorrelationMatrixElement(2,4),sb_sample->GetCorrelationMatrixElement(3,4)};
  
  TGraphErrors* tgmc = GetMPVGraph(mc_params,mc_params_error,mc_correlations);
  tgmc->SetTitle("MC");
  tgmc->SetFillColorAlpha(98,0.5);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(tgmc,"3");
  mg->Add(tgd,"3");

  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#mu_{S} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();

  mg->GetYaxis()->SetTitleOffset(0.9);
  mg->Draw("ap");
  gPad->BuildLegend(0.7,0.7,0.85,0.85,"","f");

  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update();//gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_mpv_comparative.pdf");

  //plot incoherent-coherent for MC
  TGraphErrors* tg_ts_mpv = sb_sample->GetTrueSignal()->GetImpvGraph();
  tg_ts_mpv->SetMarkerStyle(20);
  mg->RecursiveRemove(tgd);
  mg->Add(tg_ts_mpv,"p");
  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#mu_{S} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(0.9);
  mg->Draw("ap");

  TLegend* lg = new TLegend(0.6,0.6,0.86,0.86);
  lg->AddEntry(tgmc,"MC Coherent Fit", "f");
  lg->AddEntry(tg_ts_mpv,"MC Incoherent Fit", "lp");
  lg->Draw();

  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update();//gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_mpv_comparative_coherent-incoherent.pdf");

  delete c;
}

//**********************************************
void GetGWPlots(CoherentSample* sb_sample, CoherentSample* sbd_sample){
//**********************************************

  TCanvas* c = new TCanvas("c","c",700,500);
  
  //data
  double d_params[2]       = {sbd_sample->GetSignal()->GetCgwA().first,  sbd_sample->GetSignal()->GetCgwB().first };
  double d_params_error[2] = {sbd_sample->GetSignal()->GetCgwA().second, sbd_sample->GetSignal()->GetCgwB().second};
  double d_correlations[2] = {1, sbd_sample->GetCorrelationMatrixElement(5,6)};
  
  TGraphErrors* tgd = GetGWGraph(d_params,d_params_error,d_correlations);
  tgd->SetTitle("Data");
  tgd->SetFillColorAlpha(9,0.5);

  //MC
  double mc_params[2]       = {sb_sample->GetSignal()->GetCgwA().first,  sb_sample->GetSignal()->GetCgwB().first };
  double mc_params_error[2] = {sb_sample->GetSignal()->GetCgwA().second, sb_sample->GetSignal()->GetCgwB().second};
  double mc_correlations[2] = {1, sb_sample->GetCorrelationMatrixElement(5,6)};
  
  TGraphErrors* tgmc = GetGWGraph(mc_params,mc_params_error,mc_correlations);
  tgmc->SetFillColorAlpha(98,0.5);
  tgmc->SetTitle("MC");

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(tgmc,"3");
  mg->Add(tgd,"3");

  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#sigma_{S,G} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();

  mg->GetYaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetRangeUser(0.07,4);
  mg->Draw("ap");
  gPad->BuildLegend(0.7,0.7,0.85,0.85,"","f");
  gPad->SetLogy();

  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update(); //gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_gw_comparative.pdf");

  //plot incoherent-coherent for MC
  TGraphErrors* tg_ts_gw = sb_sample->GetTrueSignal()->GetIgwGraph();
  tg_ts_gw->SetMarkerStyle(20);
  mg->RecursiveRemove(tgd);
  mg->Add(tg_ts_gw,"p");
  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("#sigma_{S,G} [MeV/cm]");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();

  mg->GetYaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetRangeUser(0.07,4);
  mg->Draw("ap");

  TLegend* lg = new TLegend(0.6,0.6,0.86,0.86);
  lg->AddEntry(tgmc,"MC Coherent Fit", "f");
  lg->AddEntry(tg_ts_gw,"MC Incoherent Fit", "lp");
  lg->Draw();

  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update();//gPad->WaitPrimitive();
  
  gPad->Print("plots/CF_gw_comparative_coherent-incoherent.pdf");

  delete c;
}

//**********************************************
void GetCoherentFunctions(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();

  GetLWPlots(sb_sample, sbd_sample);
  GetMPVPlots(sb_sample, sbd_sample);
  GetGWPlots(sb_sample, sbd_sample);
}

//**********************************************
void GetLikelihoodPlot(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();
  
  std::vector<double> mc_lk = sb_sample->GetLikelihoodVector();
  std::vector<double> d_lk  = sbd_sample->GetLikelihoodVector();
  std::vector<double> rr;
  for(int i = 0; i < sb_sample->GetSize(); i++)rr.push_back(sb_sample->GetRRVector()[i].first);

  TGraph* tg_mc_lk = new TGraph(rr.size(),&rr[0],&mc_lk[0]);
  tg_mc_lk->SetTitle("MC");
  tg_mc_lk->SetLineWidth(3);
  tg_mc_lk->SetLineColor(98);
  tg_mc_lk->SetMarkerColor(98);

  TGraph* tg_d_lk = new TGraph(rr.size(),&rr[0],&d_lk[0]);
  tg_d_lk->SetTitle("Data");
  tg_d_lk->SetLineWidth(3);
  tg_d_lk->SetLineColor(9);
  tg_d_lk->SetMarkerColor(9);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(tg_mc_lk,"l");
  mg->Add(tg_d_lk,"l");

  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("Likelihood Weight");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();

  mg->Draw("ap");
  gPad->BuildLegend(0.7,0.7,0.85,0.85,"","l");

  TLatex tt1;
  tt1.SetNDC();
  tt1.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();
  gPad->Update();
  
  gPad->Print("plots/CF_likelihood.pdf");

}
  
//**********************************************
void CoherentFitPlots(){        //Main function*
//**********************************************
  
  //set ProtoDUNE-SP style
  gROOT->ProcessLine(".L ../selection/protoDUNEStyle.C");

  //open files
  CoherentFit* cf_mc = new CoherentFit("/home/migue/Public/highland/files/test/mc_dedx.root",true);
  CoherentFit* cf_d  = new CoherentFit("/home/migue/Public/highland/files/test/data_dedx.root",false);

  //run coherent fit in MC
  cf_mc->CreateCoherentSamples(999);
  cf_mc->GenerateTrueMCHistograms(1,61,2,999,1,20); //parameters are Min RR, Max RR, RR slice size, chi2 cut, dEdx min histogram value, dEdx max histogram value
  cf_mc->SequentialCoherentFit(true);

  //run coherent fit in data
  cf_d->CreateCoherentSamples(999);
  cf_d->GenerateHistograms(1,61,2,999,1,20);
  cf_d->DataCoherentFit(cf_mc);

  SetStyle(cf_mc,cf_d);
  GetHistoPlots(cf_mc,cf_d);
  GetHistoPlotsWithIncoherentFits(cf_mc,cf_d);
  GetMCTrueHistoPlotsWithIncoherentFits(cf_mc);
  GetParametrizationPlots(cf_mc);
  GetCoherentFitPlots(cf_mc,cf_d);
  GetCoherentFunctions(cf_mc,cf_d);
  GetLikelihoodPlot(cf_mc,cf_d);
}
