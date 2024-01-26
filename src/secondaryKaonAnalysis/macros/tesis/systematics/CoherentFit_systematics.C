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
void PlotMPVWithSystematicErrors(CoherentFit* cf_mc, CoherentFit* cf_d){
//**********************************************

  //get signal+background samples
  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* sbd_sample = cf_d->GetSignalPlusBackgroundSample();
  
  //data
  double d_params[3]       = {sbd_sample->GetSignal()->GetCmpvA().first, sbd_sample->GetSignal()->GetCmpvB().first, sbd_sample->GetSignal()->GetCmpvC().first };
  double d_params_error[3] = {sbd_sample->GetSignal()->GetCmpvA().second,sbd_sample->GetSignal()->GetCmpvB().second,sbd_sample->GetSignal()->GetCmpvC().second};
  double d_correlations[3] = {sbd_sample->GetCorrelationMatrixElement(2,3),sbd_sample->GetCorrelationMatrixElement(2,4),sbd_sample->GetCorrelationMatrixElement(3,4)};
  
  TGraphErrors* tgd = GetMPVGraph(d_params,d_params_error,d_correlations);
  tgd->SetTitle("Data");
  tgd->SetFillColorAlpha(9,0.5);

  //MC
  double mean_A = cf_mc->GetmpvASystHisto()->GetMean();
  double mean_B = cf_mc->GetmpvBSystHisto()->GetMean();
  double mean_C = cf_mc->GetmpvCSystHisto()->GetMean();
  double syst_error_A = cf_mc->GetmpvASystHisto()->GetRMS()/mean_A;
  double syst_error_B = cf_mc->GetmpvBSystHisto()->GetRMS()/mean_B;
  double syst_error_C = cf_mc->GetmpvCSystHisto()->GetRMS()/mean_C;
  double stat_error_A = sb_sample->GetSignal()->GetCmpvA().second/sb_sample->GetSignal()->GetCmpvA().first;
  double stat_error_B = sb_sample->GetSignal()->GetCmpvB().second/sb_sample->GetSignal()->GetCmpvB().first;
  double stat_error_C = sb_sample->GetSignal()->GetCmpvC().second/sb_sample->GetSignal()->GetCmpvC().first;
  double total_error_A = sqrt(pow(stat_error_A,2)+pow(syst_error_A,2))*mean_A;
  double total_error_B = sqrt(pow(stat_error_B,2)+pow(syst_error_B,2))*mean_B;
  double total_error_C = sqrt(pow(stat_error_C,2)+pow(syst_error_C,2))*mean_C;

  //stat errors
  double mc_params[3]       = {sb_sample->GetSignal()->GetCmpvA().first, sb_sample->GetSignal()->GetCmpvB().first, sb_sample->GetSignal()->GetCmpvC().first };
  double mc_params_error[3] = {sb_sample->GetSignal()->GetCmpvA().second,sb_sample->GetSignal()->GetCmpvB().second,sb_sample->GetSignal()->GetCmpvC().second};
  double mc_correlations[3] = {sb_sample->GetCorrelationMatrixElement(2,3),sb_sample->GetCorrelationMatrixElement(2,4),sb_sample->GetCorrelationMatrixElement(3,4)};
  
  TGraphErrors* tgmc_stat = GetMPVGraph(mc_params,mc_params_error,mc_correlations);
  tgmc_stat->SetTitle("MC");
  tgmc_stat->SetFillColorAlpha(98,0.5);

  //syst errors
  double syst_mc_params[3]       = {mean_A,mean_B,mean_C};
  double syst_mc_params_error[3] = {total_error_A,total_error_B,total_error_C};
  
  TGraphErrors* tgmc_syst = GetMPVGraph(syst_mc_params,syst_mc_params_error,mc_correlations);
  tgmc_syst->SetTitle("MC stat+syst");
  tgmc_syst->SetFillColorAlpha(98,0.5);

  TMultiGraph* mg = new TMultiGraph();
  mg->Add(tgmc_syst,"3");
  mg->Add(tgd,"3");

  mg->GetXaxis()->SetTitle("Residual Range [cm]");
  mg->GetYaxis()->SetTitle("MPV_{S} [MeV/cm]");
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
  
  gPad->Print("plots/CF_mpv_syst.pdf");

  //data mean substraction plot
  TCanvas* c2 = new TCanvas("c2","c2",700,500);
  c2->SetRightMargin(0.08);
  c2->SetLeftMargin(0.12);

  TGraphErrors* tg_diffdata_data = new TGraphErrors();
  TGraphErrors* tg_diffdata_mc_st = new TGraphErrors();
  TGraphErrors* tg_diffdata_mc_syst = new TGraphErrors();

  for(int i = 0; i < tgmc_syst->GetN(); i++){
    double x   = tgmc_syst->GetPointX(i);
    double y_d  = tgd->GetPointY(i)-tgd->GetPointY(i);
    double y_mc = tgmc_syst->GetPointY(i)-tgd->GetPointY(i);
    tg_diffdata_data->SetPoint(i,x,y_d);
    tg_diffdata_mc_st->SetPoint(i,x,y_mc);
    tg_diffdata_mc_syst->SetPoint(i,x,y_mc);
    tg_diffdata_data->SetPointError(i,0,tgd->GetErrorY(i));
    tg_diffdata_mc_st->SetPointError(i,0,tgmc_stat->GetErrorY(i));
    tg_diffdata_mc_syst->SetPointError(i,0,tgmc_syst->GetErrorY(i));
  }
  tg_diffdata_mc_syst->SetFillColorAlpha(98,0.25);
  tg_diffdata_mc_st->SetFillColorAlpha(98,0.5);
  tg_diffdata_data->SetFillColorAlpha(9,0.5);

  TMultiGraph* mg2 = new TMultiGraph();
  mg2->Add(tg_diffdata_mc_syst,"3");
  mg2->Add(tg_diffdata_mc_st,"3");
  mg2->Add(tg_diffdata_data,"3");

  mg2->GetXaxis()->SetTitle("Residual Range [cm]");
  mg2->GetYaxis()->SetTitle("(MPV_{S}-MPV_{S}^{data}) [MeV/cm]");
  mg2->GetXaxis()->CenterTitle();
  mg2->GetYaxis()->CenterTitle();
  mg2->GetYaxis()->SetTitleOffset(1.2);
  
  mg2->Draw("a");

  tg_diffdata_mc_syst->SetTitle("MC stat+syst");
  tg_diffdata_data->SetTitle("Data");
  tg_diffdata_mc_st->SetTitle("MC stat");
  
  gPad->BuildLegend(0.5,0.2,0.85,0.45,"","f");
  tt1.DrawLatex(0.12,0.94,"#bf{DUNE:ProtoDUNE-SP}");

  gPad->RedrawAxis();gPad->Update();
  gPad->Print("plots/CF_mpv_syst_diff.pdf");
}

//**********************************************
void CoherentFit_systematics(){
//**********************************************
  
  //set ProtoDUNE-SP style
  gROOT->ProcessLine(".L ../selection/protoDUNEStyle.C");

  //open files
  CoherentFit* cf_mc = new CoherentFit("/dune/app/users/miagarc/technical_note/files/systematics/all_systematics_dedx/all_systematics_dedx_merged.root",true);
  CoherentFit* cf_d  = new CoherentFit("/dune/app/users/miagarc/technical_note/files/data/data_dedx.root",false);

  //run coherent fit in MC
  cf_mc->CreateCoherentSamples(999);
  cf_mc->GenerateTrueMCHistograms(1,61,2,999,1,20); //parameters are Min RR, Max RR, RR slice size, chi2 cut, dEdx min histogram value, dEdx max histogram value
  cf_mc->SequentialCoherentFit(true);

  //run coherent fit in data
  cf_d->CreateCoherentSamples(999);
  cf_d->GenerateHistograms(1,61,2,999,1,20);
  cf_d->DataCoherentFit(cf_mc);

  //propagate systematic errors
  cf_mc->PropagateSystematicErrors("CoherentFit_systematic_results.root",true,true); //output filename, propagate weight systematics, propagate variation systematics

  //make plots
  PlotMPVWithSystematicErrors(cf_mc,cf_d);
}
