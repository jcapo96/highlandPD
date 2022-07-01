#include "CoherentSample.hxx"
#include "CoherentFitUtils.hxx"

#include "TVector.h"
#include "TParameter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "Math/Math.h"

#include <iostream>

static CoherentSample* CSMinuit;

//********************************************************************
CoherentSample::CoherentSample(){
//********************************************************************

  fSignal = NULL;
  fBackground = NULL;
  fTrueSignal = NULL;
  fTrueBackground = NULL;

  fType = SampleTypeEnum::kUnassigned;
  fBackgroundModel = BackgroundModelEnum::kFitUnassigned;

  fChi2Cut = 999;
  
  fMinuit = NULL;
  
  fh.clear();
  fSystHist.clear();

  fIntegral.clear();

  fIFit.clear();
  fCFit.clear();

  fClwFit  = NULL;
  fCmpvFit = NULL;
  fCgwFit  = NULL;

  fRR.clear();
  
  fIlw.clear();
  fImpv.clear();
  fInorm.clear();
  fIgw.clear();
  fIfw.clear();

  fClwA = std::make_pair(-999.,-999.);
  fClwB = std::make_pair(-999.,-999.);
  fClwC = std::make_pair(-999.,-999.);
  
  fCmpvA = std::make_pair(-999.,-999.);
  fCmpvB = std::make_pair(-999.,-999.);
  fCmpvC = std::make_pair(-999.,-999.);
  fCmpvD = std::make_pair(-999.,-999.);
  fCmpvR = std::make_pair(-999.,-999.);

  fCgwA = std::make_pair(-999.,-999.);
  fCgwB = std::make_pair(-999.,-999.);
  fCgwC = std::make_pair(-999.,-999.);

  fClwQa  = std::make_pair(0,0);
  fCshift = std::make_pair(0,0);
  fCgwQa  = std::make_pair(0,0);

  fCnorm.clear();
}

//********************************************************************
CoherentSample::CoherentSample(CoherentSample::SampleTypeEnum type){
//********************************************************************

  CoherentSample();
  fType = type;
}

//********************************************************************
CoherentSample::CoherentSample(const CoherentSample &c){
//********************************************************************

  fSignal = c.GetSignal();
  fBackground = c.GetBackground();
  fTrueSignal = c.GetTrueSignal();
  fTrueBackground = c.GetTrueBackground();

  fType = c.GetSampleType();
  fBackgroundModel = c.GetBackgroundModel();

  fChi2Cut = c.GetChi2Cut();
  
  fMinuit = NULL;
  
  fh = c.GetHistVector();
  fSystHist = c.GetSystHistVector();

  fIntegral = c.GetIntegralVector();

  fIFit = c.GetIFitVector();
  fCFit = c.GetCFitVector();

  fClwFit  = c.GetClwFit();
  fCmpvFit = c.GetCmpvFit();
  fCgwFit  = c.GetCgwFit();

  fRR = c.GetRRVector();
  
  fIlw = c.GetIlwVector();
  fImpv = c.GetImpvVector();
  fInorm = c.GetInormVector();
  fIgw = c.GetIgwVector();
  fIfw = c.GetIfwVector();

  fClwA = c.GetClwA();
  fClwB = c.GetClwB();
  fClwC = c.GetClwC();
  
  fCmpvA = c.GetCmpvA();
  fCmpvB = c.GetCmpvB();
  fCmpvC = c.GetCmpvC();
  fCmpvD = c.GetCmpvD();
  fCmpvR = c.GetCmpvR();

  fCgwA = c.GetCgwA();
  fCgwB = c.GetCgwB();
  fCgwC = c.GetCgwC();

  fClwQa  = c.GetClwQa();
  fCshift = c.GetCshift();
  fCgwQa  = c.GetCgwQa();

  fCnorm = c.GetCnorm();
  
}

//********************************************************************
void CoherentSample::WriteToRootFile(const std::string& filename){
//********************************************************************

  TFile* rfile = new TFile(filename.c_str(),"NEW");

  TParameter<double> par("Chi2Cut",fChi2Cut);
  rfile->WriteObject(&par,"Chi2Cut");
  
  rfile->WriteObject(&fh,"hist_vector");

  rfile->WriteObject(&fIntegral,"Int_vector");
  rfile->WriteObject(&fIIntegral,"IInt_vector");
  rfile->WriteObject(&fCIntegral,"CInt_vector");
  
  rfile->WriteObject(&fIFit,"IFit_vector");
  rfile->WriteObject(&fCFit,"CFit_vector");

  if(fClwFit)fClwFit->Write();
  if(fCmpvFit)fCmpvFit->Write();
  if(fCgwFit)fCgwFit->Write();
  
  rfile->WriteObject(&fRR,"RR_vector");

  rfile->WriteObject(&fIlw  ,"Ilw_vector");
  rfile->WriteObject(&fImpv ,"Impv_vector");
  rfile->WriteObject(&fInorm,"Inorm_vector");
  rfile->WriteObject(&fIgw  ,"Igw_vector");
  rfile->WriteObject(&fIfw  ,"Ifw_vector");

  rfile->WriteObject(&fClwA,"ClwA");
  rfile->WriteObject(&fClwB,"ClwB");
  rfile->WriteObject(&fClwC,"ClwC");

  rfile->WriteObject(&fCmpvA,"CmpvA");
  rfile->WriteObject(&fCmpvB,"CmpvB");
  rfile->WriteObject(&fCmpvC,"CmpvC");
  rfile->WriteObject(&fCmpvD,"CmpvD");
  rfile->WriteObject(&fCmpvR,"CmpvR");

  rfile->WriteObject(&fCshift,"Cshift");
  
  rfile->WriteObject(&fCgwA,"CgwA");
  rfile->WriteObject(&fCgwB,"CgwB");
  rfile->WriteObject(&fCgwC,"CgwC");
  
  rfile->WriteObject(&fCnorm,"Cnorm");
    
  rfile->Close();
}

//********************************************************************
void CoherentSample::ReadFromRootFile(const std::string& filename){
//********************************************************************

  TFile* rfile = TFile::Open(filename.c_str(),"READ");

  TParameter<double> *par;
  rfile->GetObject("Chi2Cut",par);
  if(par)fChi2Cut = par->GetVal();
  
  std::vector<TH1F*> *h;
  rfile->GetObject("hist_vector",h);
  if(h)fh = *h;
  
  std::vector<double> *d;
  rfile->GetObject("Int_vector",d);
  if(d)fIntegral = *d;
  rfile->GetObject("IInt_vector",d);
  if(d)fIIntegral = *d;
  rfile->GetObject("CInt_vector",d);
  if(d)fCIntegral = *d;
  
  std::vector<TF1*> *f;
  rfile->GetObject("IFit_vector",f);
  if(f)fIFit = *f;
  rfile->GetObject("CFit_vector",f);
  if(f)fCFit = *f;

  fClwFit  = (TF1*)rfile->Get("lw_par");
  fCmpvFit = (TF1*)rfile->Get("mpv_par");
  fCgwFit  = (TF1*)rfile->Get("gw_par");
  
  std::vector<std::pair<double,double>> *vp;
  rfile->GetObject("RR_vector",vp);
  if(vp)fRR = *vp;

  rfile->GetObject("Ilw_vector",vp);
  if(vp)fIlw = *vp;
  rfile->GetObject("Impv_vector",vp);
  if(vp)fImpv = *vp;
  rfile->GetObject("Inorm_vector",vp);
  if(vp)fInorm = *vp;
  rfile->GetObject("Igw_vector",vp);
  if(vp)fIgw = *vp;
  rfile->GetObject("Ifw_vector",vp);
  if(vp)fIfw = *vp;

  std::pair<double,double> *p;
  rfile->GetObject("ClwA",p);
  if(p)fClwA = *p;
  rfile->GetObject("ClwB",p);
  if(p)fClwB = *p;
  rfile->GetObject("ClwC",p);
  if(p)fClwC = *p;

  rfile->GetObject("CmpvA",p);
  if(p)fCmpvA = *p;
  rfile->GetObject("CmpvB",p);
  if(p)fCmpvB = *p;
  rfile->GetObject("CmpvC",p);
  if(p)fCmpvC = *p;
  rfile->GetObject("CmpvD",p);
  if(p)fCmpvD = *p;
  rfile->GetObject("CmpvR",p);
  if(p)fCmpvR = *p;

  rfile->GetObject("Cshift",p);
  if(p)fCshift = *p;

  rfile->GetObject("CgwA",p);
  if(p)fCgwA = *p;
  rfile->GetObject("CgwB",p);
  if(p)fCgwB = *p;
  rfile->GetObject("CgwC",p);
  if(p)fCgwC = *p;

  rfile->GetObject("Cnorm",vp);
  if(vp)fCnorm = *vp;
  
  rfile->Close();
}

//********************************************************************
void CoherentSample::NormalizeHistograms(){
//********************************************************************
  
  if(fh.empty()){
    std::cout << "no histogram stored yet" << std::endl;
    return;
  }

  if((int)fh[0]->Integral("width") <= 1){
    std::cout << "histograms already normalized" << std::endl;
    return;
  }
  
  if(fIntegral.empty())ComputeIntegral();
  
  for(int i = 0; i < (int)fh.size(); i++)
    fh[i]->Scale(1/fIntegral[i]);
}

//********************************************************************
void CoherentSample::NormalizeHistograms(std::vector<double> Integral){
//********************************************************************

  if((int)fh[0]->Integral("width") <= 1){
    std::cout << "histograms already normalized" << std::endl;
    return;
  }
  
  if(fh.size() > Integral.size()){
    std::cout << "new integral vector is smaller than histogram vector" << std::endl;
    std::cout << "cannot normalize it" << std::endl;
    std::exit(1);
  }
  SetIntegralVector(Integral);
  NormalizeHistograms();
  
}

//********************************************************************
void CoherentSample::ComputeIntegral(){
//********************************************************************

  if(fh.empty()){
    std::cout << "no histogram stored yet" << std::endl;
    return;
  }
  
  if((int)fh[0]->Integral("width") == 1){
    std::cout << "histograms already normalized" << std::endl;
    return;
  }

  for(int i = 0; i < (int)fh.size(); i++)
    fIntegral.push_back(fh[i]->Integral("width"));
}

//********************************************************************
void CoherentSample::StoreCoherentFits(){
//********************************************************************

  if(fType != SampleTypeEnum::kSignalPlusBackground){

    CoherentSample* satb = NULL;
    if(fType == SampleTypeEnum::kTrueBackground)satb = GetTrueSignal();
    if(fType == SampleTypeEnum::kBackground)    satb = GetSignal();
    
    if(fBackgroundModel == BackgroundModelEnum::kQuadraticWidths){
      fClwFit = new TF1("lw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fClwFit->SetParameters(satb->GetClwA().first,satb->GetClwB().first,satb->GetClwC().first, fClwA.first, fClwB.first, fClwC.first);
    }
    else{
      fClwFit = new TF1("lw_par" ,CoherentFitUtils::ABCParametrization,0,60,4);
      fClwFit->SetParameters(fClwA.first,fClwB.first,fClwC.first,fClwQa.first);
    }

    fCmpvFit = new TF1("mpv_par",CoherentFitUtils::ABCDRParametrization,0,60,6);
    fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,fCshift.first);
    
    if(fBackgroundModel == BackgroundModelEnum::kQuadraticWidths){
      fCgwFit = new TF1("gw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fCgwFit->SetParameters(satb->GetCgwA().first,satb->GetCgwB().first,satb->GetCgwC().first, fCgwA.first, fCgwB.first, fCgwC.first);
    }
    else{
      fCgwFit = new TF1("gw_par" ,CoherentFitUtils::ABCParametrization,0,60,4);
      fCgwFit->SetParameters(fCgwA.first,fCgwB.first,fCgwC.first,fCgwQa.first);
    }
    
    //store individual histogram fits
    for(int i = 0; i < (int)fRR.size(); i++){
      std::stringstream ssi;
      ssi << i;
      TF1* f = new TF1(("CF_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);

      double lw,mpv,norm,gw;
      
      lw   = fClwFit->Eval(fRR[i].first);
      mpv  = fCmpvFit->Eval(fRR[i].first);
      norm = fCnorm[0].first;
      gw   = fCgwFit->Eval(fRR[i].first);

      f->SetParameters(lw,mpv,norm,gw);		       
      fCFit.push_back(f);
      //fCIntegral.push_back(CoherentFitUtils::GetFunctionNormalizationInsideHistogramBoundaries(fh[i],f));
    }
  }
  else{
    //store parametric functions
    GetSignal()->StoreCoherentFits();
    GetBackground()->StoreCoherentFits();

    //store individual histogram fits
    for(int i = 0; i < (int)fh.size(); i++){
      std::stringstream ssi;
      ssi << i;
      TF1* f = new TF1(("CF_b"+ssi.str()+"").c_str(),CoherentFitUtils::DoubleLangaus,1,30,8);

      double slw, smpv, snorm, sgw, blw, bmpv, bnorm, bgw;
      
      slw   = GetSignal()->GetClwFit()->Eval(fRR[i].first);
      smpv  = GetSignal()->GetCmpvFit()->Eval(fRR[i].first);
      snorm = GetSignal()->GetCnorm()[0].first;
      sgw   = GetSignal()->GetCgwFit()->Eval(fRR[i].first);
      blw   = GetBackground()->GetClwFit()->Eval(fRR[i].first);
      bmpv  = GetBackground()->GetCmpvFit()->Eval(fRR[i].first);
      bnorm = GetBackground()->GetCnorm()[0].first;
      bgw   = GetBackground()->GetCgwFit()->Eval(fRR[i].first);
      
      f->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
      fCFit.push_back(f);
      fCIntegral.push_back(CoherentFitUtils::GetFunctionNormalizationInsideHistogramBoundaries(fh[i],f));
    }
  }
}

//********************************************************************
void CoherentSample::SetCFitParameters(const CoherentSample* sample){
//********************************************************************

  CoherentSample* copysample = NULL;
  if(sample->GetSampleType() == SampleTypeEnum::kTrueSignal)copysample = GetSignal();
  else if (sample->GetSampleType() == SampleTypeEnum::kTrueBackground)copysample = GetBackground();
  else if (sample->GetSampleType() == SampleTypeEnum::kSignal || sample->GetSampleType() == SampleTypeEnum::kBackground)copysample = this;
  
  copysample->SetClwA(sample->GetClwA());
  copysample->SetClwB(sample->GetClwB());
  copysample->SetClwC(sample->GetClwC());
  copysample->SetCmpvA(sample->GetCmpvA());
  copysample->SetCmpvB(sample->GetCmpvB());
  copysample->SetCmpvC(sample->GetCmpvC());
  copysample->SetCmpvD(sample->GetCmpvD());
  copysample->SetCmpvR(sample->GetCmpvR());
  copysample->SetCgwA(sample->GetCgwA());
  copysample->SetCgwB(sample->GetCgwB());
  copysample->SetCgwC(sample->GetCgwC());
  copysample->SetCnorm(sample->GetCnorm());
  copysample->SetCshift(sample->GetCshift());
}

//********************************************************************
void CoherentSample::SetCFitParametersWithVariations(const CoherentSample* sample, TRandom3* r, const double sigma,
						     const bool apply_all_var,
						     const bool apply_only_lw_var, 
						     const bool apply_only_mpv_var, 
						     const bool apply_only_gw_var,
						     const bool apply_only_shift_var,
						     const bool apply_only_norm_var){
//********************************************************************

  if(apply_all_var && 
     (apply_only_lw_var || apply_only_mpv_var || apply_only_gw_var || 
      apply_only_shift_var || apply_only_norm_var)){
    std::cout << "bad parameter configuration, check it!" << std::endl;
    std::exit(1);
  }
  
  //fill default values
  this->fClwA.first     = sample->GetClwA().first;
  this->fClwB.first     = sample->GetClwB().first;
  this->fClwC.first     = sample->GetClwC().first;
  this->fCmpvA.first    = sample->GetCmpvA().first;
  this->fCmpvB.first    = sample->GetCmpvB().first;
  this->fCmpvC.first    = sample->GetCmpvC().first;
  this->fCmpvD.first    = sample->GetCmpvD().first;
  this->fCmpvR.first    = sample->GetCmpvR().first;
  this->fCgwA.first     = sample->GetCgwA().first;
  this->fCgwB.first     = sample->GetCgwB().first;
  this->fCgwC.first     = sample->GetCgwC().first;
  this->fCshift.first   = sample->GetCshift().first;
  this->fCnorm[0].first = sample->GetCnorm()[0].first;
  
  //apply desired variations
  if(apply_all_var || apply_only_lw_var){
    this->fClwA.first     = r->Gaus(sample->GetClwA().first,sample->GetClwA().first*sigma);
    this->fClwB.first     = r->Gaus(sample->GetClwB().first,sample->GetClwB().first*sigma);
    this->fClwC.first     = r->Gaus(sample->GetClwC().first,sample->GetClwC().first*sigma);
  }
  if(apply_all_var || apply_only_mpv_var){
    this->fCmpvA.first    = r->Gaus(sample->GetCmpvA().first,sample->GetCmpvA().first*sigma);
    //this->fCmpvB.first    = r->Gaus(sample->GetCmpvB().first,sample->GetCmpvB().first*sigma);
    //this->fCmpvC.first    = r->Gaus(sample->GetCmpvC().first,sample->GetCmpvC().first*sigma);
    //this->fCmpvD.first    = r->Gaus(sample->GetCmpvD().first,sample->GetCmpvD().first*sigma);
    //this->fCmpvR.first    = r->Gaus(sample->GetCmpvR().first,sample->GetCmpvR().first*sigma);
  }
  if(apply_all_var || apply_only_gw_var){
    this->fCgwA.first     = r->Gaus(sample->GetCgwA().first,sample->GetCgwA().first*sigma);
    this->fCgwB.first     = r->Gaus(sample->GetCgwB().first,sample->GetCgwB().first*sigma);
    this->fCgwC.first     = r->Gaus(sample->GetCgwC().first,sample->GetCgwC().first*sigma);
  }
  if(apply_all_var || apply_only_shift_var)
    this->fCshift.first   = r->Gaus(sample->GetCshift().first,sample->GetCshift().first*sigma);
  if(apply_all_var || apply_only_norm_var)
    this->fCnorm[0].first = r->Gaus(sample->GetCnorm()[0].first,sample->GetCnorm()[0].first*sigma);

  std::cout << this->fClwA.first << " " << this->fCmpvA.first << std::endl;
}
  
//********************************************************************
void CoherentSample::SequentialCoherentFit(){
//********************************************************************

  if(fh.empty()){
    std::cout << "No histograms to fit" << std::endl;
    std::exit(1);
  }

  std::cout << "sequential coherent fit" << std::endl;

  //if(fType == SampleTypeEnum::kBackground){
  //  for(int i = 0; i < (int)fh.size(); i++)fh[i]->Rebin(2);
  //}

  IncoherentFit();
  GetInitialParValuesForCoherentFit();
  CoherentFit();
}

//********************************************************************
void CoherentSample::IncoherentFit(){
//********************************************************************

  std::cout << "incoherent fit" << std::endl;
  
  if(fType != SampleTypeEnum::kSignalPlusBackground){
    for(int i = 0; i < (int)fh.size(); i++){
      fIFit.push_back(CoherentFitUtils::LangausFit(fh[i]));
      //save parameters and errors
      fIlw.push_back(  std::make_pair(fIFit.back()->GetParameter(0),fIFit.back()->GetParError(0)));
      fImpv.push_back( std::make_pair(fIFit.back()->GetParameter(1),fIFit.back()->GetParError(1)));
      fInorm.push_back(std::make_pair(fIFit.back()->GetParameter(2),fIFit.back()->GetParError(2)));
      fIgw.push_back(  std::make_pair(fIFit.back()->GetParameter(3),fIFit.back()->GetParError(3)));
      fIfw.push_back(std::make_pair(sqrt(pow(fIlw.back().first ,2)+pow(fIgw.back().first ,2)),
				    sqrt(pow(fIlw.back().second,2)+pow(fIgw.back().second,2))));
      fIIntegral.push_back(CoherentFitUtils::GetFunctionNormalizationInsideHistogramBoundaries(fh[i],fIFit[i]));
    }
  }
  else if(fType == SampleTypeEnum::kSignalPlusBackground){
    return;
  }
}

//********************************************************************
void CoherentSample::GetInitialParValuesForCoherentFit(bool equal_weights, bool draw_fits){
//********************************************************************

  std::cout << "Initial par values" << std::endl;
  
  if(fType == SampleTypeEnum::kTrueSignal || fType == SampleTypeEnum::kTrueBackground){
    CoherentFitUtils::GetABCParametrization(fClwA.first,fClwB.first,fClwC.first,
					    fRR,fIlw,true,draw_fits);
    CoherentFitUtils::GetABCParametrization(fCgwA.first,fCgwB.first,fCgwC.first,
					    fRR,fIgw,equal_weights,draw_fits);
    CoherentFitUtils::GetABCDRParametrization(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,
					      fRR,fImpv,equal_weights,draw_fits);
  }
  else if(fType == SampleTypeEnum::kSignalPlusBackground){
    SetCFitParameters(GetTrueSignal());
    SetCFitParameters(GetTrueBackground());
  }
}

//********************************************************************
void CoherentSample::CoherentFit(){
//********************************************************************

  std::cout << "coherent fit" << std::endl;
  
  if(fType == SampleTypeEnum::kTrueSignal){
    CoherentFitSignal();
  }
  else if(fType == SampleTypeEnum::kTrueBackground){
    if(fBackgroundModel == BackgroundModelEnum::kShift)                CoherentFitBackgroundShift();
    else if (fBackgroundModel == BackgroundModelEnum::k3Par)           CoherentFitBackground3Par();
    else if (fBackgroundModel == BackgroundModelEnum::kQuadraticWidths)CoherentFitBackgroundQuadraticWidths();
    else                                                               CoherentFitBackgroundAllFree();
  }
  else if(fType == SampleTypeEnum::kSignalPlusBackground){
    if(fBackgroundModel == BackgroundModelEnum::kShift)CoherentFitSignalPlusBackgroundShift();
    else if (fBackgroundModel == BackgroundModelEnum::kQuadraticWidths)CoherentFitSignalPlusBackgroundQuadraticWidths();
    else CoherentFitSignalPlusBackgroundAllFree();
  }
}

//********************************************************************
TGraphErrors* CoherentSample::GetMPVErrorBand(){
//********************************************************************

  //get covariance matrix
  const int ndim = fMinuit->GetNumPars();
  double matrix[ndim][ndim] = {0};
  fMinuit->mnemat(&matrix[0][0],ndim);

  //get hessian
  TF1* J[5];
  J[0] = new TF1("derA" ,CoherentFitUtils::ABCDRDerivativeA,0,60,6);
  J[1] = new TF1("derB" ,CoherentFitUtils::ABCDRDerivativeB,0,60,6);
  J[2] = new TF1("derC" ,CoherentFitUtils::ABCDRDerivativeC,0,60,6);
  J[3] = new TF1("derD" ,CoherentFitUtils::ABCDRDerivativeD,0,60,6);
  J[4] = new TF1("derR" ,CoherentFitUtils::ABCDRDerivativeR,0,60,6);
  for(int i = 0; i < 5; i++)J[i]->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,0);
  TF1* mpvFit = NULL;
  if(fType==SampleTypeEnum::kSignalPlusBackground)mpvFit = fSignal->GetCmpvFit();
  else if(fType==SampleTypeEnum::kTrueSignal)mpvFit = fCmpvFit;
  else return NULL;
  
  //compute error at each RR
  std::vector<double> rr,rr_error,mpv,mpv_error;
  rr.clear();
  rr_error.clear();
  mpv.clear();
  mpv_error.clear();
  double error    = 0;
  double aux      = 0;
  double v_aux[5] = {0};
  for(int irr = 0; irr < fRR.size(); irr++){
    rr.push_back(fRR[irr].first);
    rr_error.push_back(fRR[irr].second);
    mpv.push_back(mpvFit->Eval(fRR[irr].first));
    for(int icolumn = 3; icolumn < 8; icolumn++){
      for(int irow = 3; irow < 8; irow++)
	aux += J[irow-3]->Eval(fRR[irr].first)*matrix[irow][icolumn];
      v_aux[icolumn-3] = aux;
      aux = 0;
    }
    for(int i = 0; i < 5; i++)error += v_aux[i]*J[i]->Eval(fRR[irr].first);
    mpv_error.push_back(error);
    std::cout << fRR[irr].first << " " << error << std::endl;
    error = 0;
  }

  return new TGraphErrors(mpv.size(),&rr[0],&mpv[0],&rr_error[0],&mpv_error[0]);
}
  
//********************************************************************
void CoherentSample::CoherentFitSignal(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 11+1;
  const int NBINS = fh.size();
  const int NPAR = CPAR;// + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignal);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {fClwA.first, fClwB.first, fClwC.first,
			   fCmpvA.first, fCmpvB.first, fCmpvC.first, fCmpvD.first, fCmpvR.first,
			   fCgwA.first, fCgwB.first, fCgwC.first};
  Double_t step[NPAR]   = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fInorm[ipar].first;
    step[ipar+CPAR]   = 10;
    }*/
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = mean_norm / fInorm.size();
    
  fMinuit->mnparm(0,  "lw A",  vstart[0],  step[0],  0, 0, ierflg);
  fMinuit->mnparm(1,  "lw B",  vstart[1],  step[1],  0, 0, ierflg);
  fMinuit->mnparm(2,  "lw C",  vstart[2],  step[2],  0, 0, ierflg);
  fMinuit->mnparm(3,  "mpv A", vstart[3],  step[3],  0, 0, ierflg);
  fMinuit->mnparm(4,  "mpv B", vstart[4],  step[4],  0, 0, ierflg);
  fMinuit->mnparm(5,  "mpv C", vstart[5],  step[5],  0, 0, ierflg);
  fMinuit->mnparm(6,  "mpv D", vstart[6],  step[6],  0, 0, ierflg);
  fMinuit->mnparm(7,  "mpv R", vstart[7],  step[7],  0, 0, ierflg);
  fMinuit->mnparm(8,  "gw A",  vstart[8],  step[8],  0, 0, ierflg);
  fMinuit->mnparm(9,  "gw B",  vstart[9],  step[9],  0, 0, ierflg);
  fMinuit->mnparm(10, "gw C",  vstart[10], step[10], 0, 0, ierflg);
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 0, ierflg);
		    }*/
  fMinuit->mnparm(11, "norm",  mean_norm, 0.01, 0, 1, ierflg);
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fClwA.first,  fClwA.second);
  fMinuit->GetParameter(1,  fClwB.first,  fClwB.second);
  fMinuit->GetParameter(2,  fClwC.first,  fClwC.second);
  fMinuit->GetParameter(3,  fCmpvA.first, fCmpvA.second);
  fMinuit->GetParameter(4,  fCmpvB.first, fCmpvB.second);
  fMinuit->GetParameter(5,  fCmpvC.first, fCmpvC.second);
  fMinuit->GetParameter(6,  fCmpvD.first, fCmpvD.second);
  fMinuit->GetParameter(7,  fCmpvR.first, fCmpvR.second);
  fMinuit->GetParameter(8,  fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(9,  fCgwB.first,  fCgwB.second);
  fMinuit->GetParameter(10, fCgwC.first,  fCgwC.second);
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCnorm.push_back(std::make_pair(norm, norm_error));
    }*/
  double norm, norm_error;
  fMinuit->GetParameter(11, norm, norm_error);
  fCnorm.push_back(std::make_pair(norm, norm_error));

  StoreCoherentFits();
  
}

//********************************************************************
void CoherentSample::fcnSignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 11;
  
  //set each parameter
  double lwA  = par[0];
  double lwB  = par[1];
  double lwC  = par[2];
  double mpvA = par[3];
  double mpvB = par[4];
  double mpvC = par[5];
  double mpvD = par[6];
  double mpvR = par[7];
  double gwA  = par[8];
  double gwB  = par[9];
  double gwC  = par[10];
  double norm = par[11];
  
  //define function to fit to each histogram
  double lw,mpv,gw;//norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double lw_par[]  = {lwA,lwB,lwC,0};
    double mpv_par[] = {mpvA,mpvB,mpvC,mpvD,mpvR,0};
    double gw_par[]  = {gwA,gwB,gwC,0};
    double RR_par[]  = {CSMinuit->fRR[ibin].first};
    
    lw   = CoherentFitUtils::ABCParametrization(RR_par,lw_par);
    mpv  = CoherentFitUtils::ABCDRParametrization(RR_par,mpv_par);
    gw   = CoherentFitUtils::ABCParametrization(RR_par,gw_par);
    //norm = par[CPAR+ibin];

    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//BACKGROUND MODELS

//********************************************************************
void CoherentSample::CoherentFitBackgroundAllFree(){
  //********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 11;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundAllFree);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {fClwA.first, fClwB.first, fClwC.first,
			   fCmpvA.first, fCmpvB.first, fCmpvC.first, fCmpvD.first, fCmpvR.first,
			   fCgwA.first, fCgwB.first, fCgwC.first};
  Double_t step[NPAR]   = {         0.1,          0.1,          0.1,           0.1,           0.1,           0.1,           0.1,           0.1,           0.1,          0.1,         0.1};
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fInorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
  }
  fMinuit->mnparm(0,  "lwA" ,  vstart[0] , step[0],  1, 20, ierflg);  ;//0, 0, ierflg);
  fMinuit->mnparm(1,  "lwB" ,  vstart[1] , step[1],  0, 0, ierflg);
  fMinuit->mnparm(2,  "lwC" ,  vstart[2] , step[2],  0.05, 1, ierflg);//0, 0, ierflg);
  fMinuit->mnparm(3,  "mpvA",  vstart[3] , step[3],  0, 0, ierflg);
  fMinuit->mnparm(4,  "mpvB",  vstart[4] , step[5],  0, 0, ierflg);
  fMinuit->mnparm(5,  "mpvC",  vstart[5] , step[5],  0, 0, ierflg);
  fMinuit->mnparm(6,  "mpvD",  vstart[6] , step[6],  0, 0, ierflg);
  fMinuit->mnparm(7,  "mpvR",  vstart[7] , step[7],  0, 0, ierflg);
  fMinuit->mnparm(8,  "gwA" ,  vstart[8] , step[8],  0, 0, ierflg);
  fMinuit->mnparm(9,  "gwB" ,  vstart[9] , step[9],  0, 0, ierflg);
  fMinuit->mnparm(10, "gwC" ,  vstart[10], step[10], 0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fClwA.first ,  fClwA.second);
  fMinuit->GetParameter(1,  fClwB.first ,  fClwB.second);
  fMinuit->GetParameter(2,  fClwC.first ,  fClwC.second);
  fMinuit->GetParameter(3,  fCmpvA.first,  fCmpvA.second);
  fMinuit->GetParameter(4,  fCmpvB.first,  fCmpvB.second);
  fMinuit->GetParameter(5,  fCmpvC.first,  fCmpvC.second);
  fMinuit->GetParameter(6,  fCmpvD.first,  fCmpvD.second);
  fMinuit->GetParameter(7,  fCmpvR.first,  fCmpvR.second);
  fMinuit->GetParameter(8,  fCgwA.first ,  fCgwA.second);
  fMinuit->GetParameter(9,  fCgwB.first ,  fCgwB.second);
  fMinuit->GetParameter(10, fCgwC.first ,  fCgwC.second);
  for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCnorm.push_back(std::make_pair(norm, norm_error));
  }

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundAllFree(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  //********************************************************************
  
  const int CPAR = 11;
  
  //set each parameter
  double lwA   = par[0];
  double lwB   = par[1];
  double lwC   = par[2];
  double mpvA  = par[3];
  double mpvB  = par[4];
  double mpvC  = par[5];
  double mpvD  = par[6];
  double mpvR  = par[7];
  double gwA   = par[8];
  double gwB   = par[9];
  double gwC   = par[10];
  
  //define function to fit to each histogram
  double lw,mpv,norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double lw_par[]  = {lwA,lwB,lwC,0};
    double mpv_par[] = {mpvA,mpvB,mpvC,mpvD,mpvR,0};
    double gw_par[]  = {gwA,gwB,gwC,0};
    double RR_par[]  = {CSMinuit->fRR[ibin].first};
    
    lw   = CoherentFitUtils::ABCParametrization(RR_par,lw_par);
    mpv  = CoherentFitUtils::ABCDRParametrization(RR_par,mpv_par);
    gw   = CoherentFitUtils::ABCParametrization(RR_par,gw_par);
    norm = par[CPAR+ibin];
    
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundAllFree(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 22;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundAllFree);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,smpvR,sgwA,sgwB,sgwC,blwA,blwB,blwC,bmpvA,bmpvB,bmpvC,bmpvD,bmpvR,bgwA,bgwB,bgwC,slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,smpvR_error,sgwA_error,sgwB_error,sgwC_error,blwA_error,blwB_error,blwC_error,bmpvA_error,bmpvB_error,bmpvC_error,bmpvD_error,bmpvR_error,bgwA_error,bgwB_error,bgwC_error;
  slwA = GetSignal()->GetClwA().first;
  slwB = GetSignal()->GetClwB().first;
  slwC = GetSignal()->GetClwC().first;
  smpvA = GetSignal()->GetCmpvA().first;
  smpvB = GetSignal()->GetCmpvB().first;
  smpvC = GetSignal()->GetCmpvC().first;
  smpvD = GetSignal()->GetCmpvD().first;
  smpvR = GetSignal()->GetCmpvR().first;
  sgwA = GetSignal()->GetCgwA().first;
  sgwB = GetSignal()->GetCgwB().first;
  sgwC = GetSignal()->GetCgwC().first;
  blwA = GetBackground()->GetClwA().first;
  blwB = GetBackground()->GetClwB().first;
  blwC = GetBackground()->GetClwC().first;
  bmpvA = GetBackground()->GetCmpvA().first;
  bmpvB = GetBackground()->GetCmpvB().first;
  bmpvC = GetBackground()->GetCmpvC().first;
  bmpvD = GetBackground()->GetCmpvD().first;
  bmpvR = GetBackground()->GetCmpvR().first;
  bgwA = GetBackground()->GetCgwA().first;
  bgwB = GetBackground()->GetCgwB().first;
  bgwC = GetBackground()->GetCgwC().first;
  Double_t vstart[NPAR] = {slwA, slwB, slwC, smpvA, smpvB, smpvC, smpvD, smpvR, sgwA, sgwB, sgwC,
			   blwA, blwB, blwC, bmpvA, bmpvB, bmpvC, bmpvD, bmpvR, bgwA, bgwB, bgwC};
  Double_t step[NPAR]   = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  std::vector<std::pair<double,double>> snorm = GetSignal()->GetCnorm();
  std::vector<std::pair<double,double>> bnorm = GetBackground()->GetCnorm();
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = (snorm[ipar].first+(1-bnorm[ipar].first))/2;
    step[ipar+CPAR]   = 0.01;
  }
  fMinuit->mnparm(0,  "s lw A" , vstart[0] ,  step[0], 0, 0, ierflg);
  fMinuit->mnparm(1,  "s lw B" , vstart[1] ,  step[1], 0, 0, ierflg);
  fMinuit->mnparm(2,  "s lw C" , vstart[2] ,  step[2], 0, 0, ierflg);
  fMinuit->mnparm(3,  "s mpv A", vstart[3] ,  step[3], 0, 0, ierflg);
  fMinuit->mnparm(4,  "s mpv B", vstart[4] ,  step[4], 0, 0, ierflg);
  fMinuit->mnparm(5,  "s mpv C", vstart[5] ,  step[5], 0, 0, ierflg);
  fMinuit->mnparm(6,  "s mpv D", vstart[6] ,  step[6], 0, 0, ierflg);
  fMinuit->mnparm(7,  "s mpv R", vstart[7] ,  step[7], 0, 0, ierflg);
  fMinuit->mnparm(8,  "s gw A" , vstart[8] ,  step[8], 0, 0, ierflg);
  fMinuit->mnparm(9,  "s gw B" , vstart[9] ,  step[9], 0, 0, ierflg);
  fMinuit->mnparm(10, "s gw C" , vstart[10], step[10], 0, 0, ierflg);
  fMinuit->mnparm(11, "b lw A" , vstart[11], step[11], 0, 0, ierflg);
  fMinuit->mnparm(12, "b lw B" , vstart[12], step[12], 0, 0, ierflg);
  fMinuit->mnparm(13, "b lw C" , vstart[13], step[13], 0, 0, ierflg);
  fMinuit->mnparm(14, "b mpv A", vstart[14], step[14], 0, 0, ierflg);
  fMinuit->mnparm(15, "b mpv B", vstart[15], step[15], 0, 0, ierflg);
  fMinuit->mnparm(16, "b mpv C", vstart[16], step[16], 0, 0, ierflg);
  fMinuit->mnparm(17, "b mpv D", vstart[17], step[17], 0, 0, ierflg);
  fMinuit->mnparm(18, "b mpv R", vstart[18], step[18], 0, 0, ierflg);
  fMinuit->mnparm(19, "b gw A" , vstart[19], step[19], 0, 0, ierflg);
  fMinuit->mnparm(20, "b gw B" , vstart[20], step[20], 0, 0, ierflg);
  fMinuit->mnparm(21, "b gw C" , vstart[21], step[21], 0, 0, ierflg);
  fMinuit->FixParameter(0);
  fMinuit->FixParameter(1);
  fMinuit->FixParameter(2);
  fMinuit->FixParameter(3);
  fMinuit->FixParameter(4);
  fMinuit->FixParameter(5);
  fMinuit->FixParameter(6);
  fMinuit->FixParameter(7);
  fMinuit->FixParameter(8);
  fMinuit->FixParameter(9);
  fMinuit->FixParameter(10);
  //fMinuit->FixParameter(11);
  //fMinuit->FixParameter(12);
  //fMinuit->FixParameter(13);
  //fMinuit->FixParameter(14);
  //fMinuit->FixParameter(15);
  //fMinuit->FixParameter(16);
  //fMinuit->FixParameter(17);
  //fMinuit->FixParameter(18);
  //fMinuit->FixParameter(19);
  //fMinuit->FixParameter(20);
  //fMinuit->FixParameter(21);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("s Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
    //fMinuit->FixParameter(ipar+CPAR);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  slwA , slwA_error);
  fMinuit->GetParameter(1,  slwB , slwB_error);
  fMinuit->GetParameter(2,  slwC , slwC_error);
  fMinuit->GetParameter(3,  smpvA, smpvA_error);
  fMinuit->GetParameter(4,  smpvB, smpvB_error);
  fMinuit->GetParameter(5,  smpvC, smpvC_error);
  fMinuit->GetParameter(6,  smpvD, smpvD_error);
  fMinuit->GetParameter(7,  smpvR, smpvR_error);
  fMinuit->GetParameter(8,  sgwA , sgwA_error);
  fMinuit->GetParameter(9,  sgwB , sgwB_error);
  fMinuit->GetParameter(10, sgwC , sgwC_error);
  fMinuit->GetParameter(11, blwA , blwA_error);
  fMinuit->GetParameter(12, blwB , blwB_error);
  fMinuit->GetParameter(13, blwC , blwC_error);
  fMinuit->GetParameter(14, bmpvA, bmpvA_error);
  fMinuit->GetParameter(15, bmpvB, bmpvB_error);
  fMinuit->GetParameter(16, bmpvC, bmpvC_error);
  fMinuit->GetParameter(17, bmpvD, bmpvD_error);
  fMinuit->GetParameter(18, bmpvR, bmpvR_error);
  fMinuit->GetParameter(19, bgwA , bgwA_error);
  fMinuit->GetParameter(20, bgwB , bgwB_error);
  fMinuit->GetParameter(21, bgwC , bgwC_error);
  std::vector<std::pair<double,double>> sn;
  std::vector<std::pair<double,double>> bn;
  for(int ipar = 0; ipar < NBINS; ipar++){
    double snorm, snorm_error;//, bnorm, bnorm_error;
    fMinuit->GetParameter(ipar+CPAR, snorm, snorm_error);
    sn.push_back(std::make_pair(snorm, snorm_error));
    //fMinuit->GetParameter(ipar+CPAR+NBINS, bnorm, bnorm_error);
    bn.push_back(std::make_pair(1-snorm, snorm_error));
  }

  GetSignal()->SetClwA(std::make_pair(slwA,slwA_error));
  GetSignal()->SetClwB(std::make_pair(slwB,slwB_error));
  GetSignal()->SetClwC(std::make_pair(slwC,slwC_error));
  GetSignal()->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  GetSignal()->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  GetSignal()->SetCmpvC(std::make_pair(smpvC,smpvR_error));
  GetSignal()->SetCmpvD(std::make_pair(smpvD,smpvC_error));
  GetSignal()->SetCmpvR(std::make_pair(smpvR,smpvD_error));
  GetSignal()->SetCgwA(std::make_pair(sgwA,sgwA_error));
  GetSignal()->SetCgwB(std::make_pair(sgwB,sgwB_error));
  GetSignal()->SetCgwC(std::make_pair(sgwC,sgwC_error));
  GetSignal()->SetCnorm(sn);

  GetBackground()->SetClwA(std::make_pair(blwA,blwA_error));
  GetBackground()->SetClwB(std::make_pair(blwB,blwB_error));
  GetBackground()->SetClwC(std::make_pair(blwC,blwC_error));
  GetBackground()->SetCmpvA(std::make_pair(bmpvA,bmpvA_error));
  GetBackground()->SetCmpvB(std::make_pair(bmpvB,bmpvB_error));
  GetBackground()->SetCmpvC(std::make_pair(bmpvC,bmpvR_error));
  GetBackground()->SetCmpvD(std::make_pair(bmpvD,bmpvC_error));
  GetBackground()->SetCmpvR(std::make_pair(bmpvR,bmpvD_error));
  GetBackground()->SetCgwA(std::make_pair(bgwA,bgwA_error));
  GetBackground()->SetCgwB(std::make_pair(bgwB,bgwB_error));
  GetBackground()->SetCgwC(std::make_pair(bgwC,bgwC_error));
  GetBackground()->SetCnorm(bn);

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundAllFree(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  const int CPAR = 22;
  //const int NBINS = CSMinuit->fh.size();
  
  //set each parameter
  double slwA  = par[0];
  double slwB  = par[1];
  double slwC  = par[2];
  double smpvA = par[3];
  double smpvB = par[4];
  double smpvC = par[5];
  double smpvD = par[6];
  double smpvR = par[7];
  double sgwA  = par[8];
  double sgwB  = par[9];
  double sgwC  = par[10];
  double blwA  = par[11];
  double blwB  = par[12];
  double blwC  = par[13];
  double bmpvA = par[14];
  double bmpvB = par[15];
  double bmpvC = par[16];
  double bmpvD = par[17];
  double bmpvR = par[18];
  double bgwA  = par[19];
  double bgwB  = par[20];
  double bgwC  = par[21];
  
  //define function to fit to each histogram
  double slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double slw_par[]  = {slwA,slwB,slwC,0};
    double smpv_par[] = {smpvA,smpvB,smpvC,smpvD,smpvR,0};
    double sgw_par[]  = {sgwA,sgwB,sgwC,0};
    double blw_par[]  = {blwA,blwB,blwC,0};
    double bmpv_par[] = {bmpvA,bmpvB,bmpvC,bmpvD,bmpvR,0};
    double bgw_par[]  = {bgwA,bgwB,bgwC,0};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    
    slw   = CoherentFitUtils::ABCParametrization(RR_par,slw_par);
    smpv  = CoherentFitUtils::ABCDRParametrization(RR_par,smpv_par);
    sgw   = CoherentFitUtils::ABCParametrization(RR_par,sgw_par);
    snorm = par[CPAR+ibin];
    blw   = CoherentFitUtils::ABCParametrization(RR_par,blw_par);
    bmpv  = CoherentFitUtils::ABCDRParametrization(RR_par,bmpv_par);
    bgw   = CoherentFitUtils::ABCParametrization(RR_par,bgw_par);
    bnorm = 1-snorm;//par[CPAR+ibin+NBINS];
    
    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    if(smpv>bmpv)
      Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
    else
      Likelihood = Likelihood + 0.1*CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//shift on mpv

//********************************************************************
void CoherentSample::CoherentFitBackgroundShift(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 7;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundShift);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {1, fClwA.first, fClwB.first, fClwC.first, fCgwA.first, fCgwB.first, fCgwC.first};
  Double_t step[NPAR]   = {0.1 ,       0.1,         0.1,         0.1,         0.1,         0.1,         0.1};
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fInorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
  }
  fMinuit->mnparm(0,  "shift",  vstart[0],  step[0],  0, 0, ierflg);
  fMinuit->mnparm(1,  "lwA"  ,  vstart[1],  step[1],  1, 20, ierflg);
  fMinuit->mnparm(2,  "lwB"  ,  vstart[2],  step[2],  0, 0, ierflg);
  fMinuit->mnparm(3,  "lwC"  ,  vstart[3],  step[3],  0.05, 1, ierflg);
  fMinuit->mnparm(4,  "gwA"  ,  vstart[4],  step[4],  0, 0, ierflg);
  fMinuit->mnparm(5,  "gwB"  ,  vstart[5],  step[5],  0, 0, ierflg);
  fMinuit->mnparm(6,  "gwC"  ,  vstart[6],  step[6],  0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fCshift.first,  fCshift.second);
  fMinuit->GetParameter(1,  fClwA.first  ,  fClwA.second);
  fMinuit->GetParameter(2,  fClwB.first  ,  fClwB.second);
  fMinuit->GetParameter(3,  fClwC.first  ,  fClwC.second);
  fMinuit->GetParameter(4,  fCgwA.first  ,  fCgwA.second);
  fMinuit->GetParameter(5,  fCgwB.first  ,  fCgwB.second);
  fMinuit->GetParameter(6,  fCgwC.first  ,  fCgwC.second);
  for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCnorm.push_back(std::make_pair(norm, norm_error));
  }

  fCmpvA = fTrueSignal->GetCmpvA();
  fCmpvB = fTrueSignal->GetCmpvB();
  fCmpvC = fTrueSignal->GetCmpvC();
  fCmpvD = fTrueSignal->GetCmpvD();
  fCmpvR = fTrueSignal->GetCmpvR();	  
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundShift(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 7;
  
  //set each parameter
  double shift = par[0];
  double lwA   = par[1];
  double lwB   = par[2];
  double lwC   = par[3];
  double mpvA  = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  double mpvB  = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  double mpvC  = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  double mpvD  = CSMinuit->GetTrueSignal()->GetCmpvD().first;
  double mpvR  = CSMinuit->GetTrueSignal()->GetCmpvR().first;
  double gwA   = par[4];
  double gwB   = par[5];
  double gwC   = par[6];

  //define function to fit to each histogram
  double lw,mpv,norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double lw_par[]   = {lwA,lwB,lwC,0};
    double mpv_par[]  = {mpvA,mpvB,mpvC,mpvD,mpvR,shift};
    double gw_par[]   = {gwA,gwB,gwC,0};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    lw   = CoherentFitUtils::ABCParametrization(RR_par,lw_par);
    mpv  = CoherentFitUtils::ABCDRParametrization(RR_par,mpv_par);
    gw   = CoherentFitUtils::ABCParametrization(RR_par,gw_par);
    norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundShift(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 18;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundShift);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,smpvR,sgwA,sgwB,sgwC,blwA,blwB,blwC,bshift,bgwA,bgwB,bgwC,slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,smpvR_error,sgwA_error,sgwB_error,sgwC_error,blwA_error,blwB_error,blwC_error,bshift_error,bgwA_error,bgwB_error,bgwC_error;
  slwA = GetSignal()->GetClwA().first;
  slwB = GetSignal()->GetClwB().first;
  slwC = GetSignal()->GetClwC().first;
  smpvA = GetSignal()->GetCmpvA().first;
  smpvB = GetSignal()->GetCmpvB().first;
  smpvC = GetSignal()->GetCmpvC().first;
  smpvD = GetSignal()->GetCmpvD().first;
  smpvR = GetSignal()->GetCmpvR().first;
  sgwA = GetSignal()->GetCgwA().first;
  sgwB = GetSignal()->GetCgwB().first;
  sgwC = GetSignal()->GetCgwC().first;
  blwA = GetBackground()->GetClwA().first;
  blwB = GetBackground()->GetClwB().first;
  blwC = GetBackground()->GetClwC().first;
  bshift = GetBackground()->GetCshift().first;
  bgwA = GetBackground()->GetCgwA().first;
  bgwB = GetBackground()->GetCgwB().first;
  bgwC = GetBackground()->GetCgwC().first;
  Double_t vstart[NPAR] = {slwA, slwB, slwC, smpvA, smpvB, smpvC, smpvD, smpvR, sgwA, sgwB, sgwC,
			   blwA, blwB, blwC, bshift,                            bgwA, bgwB, bgwC};
  Double_t step[NPAR]   = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  std::vector<std::pair<double,double>> snorm = GetSignal()->GetCnorm();
  std::vector<std::pair<double,double>> bnorm = GetBackground()->GetCnorm();
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = (snorm[ipar].first+(1-bnorm[ipar].first))/2;
    step[ipar+CPAR]   = 0.01;
  }
  fMinuit->mnparm(0,  "s lw A" , vstart[0] ,  step[0], 0, 0, ierflg);
  fMinuit->mnparm(1,  "s lw B" , vstart[1] ,  step[1], 0, 0, ierflg);
  fMinuit->mnparm(2,  "s lw C" , vstart[2] ,  step[2], 0, 0, ierflg);
  fMinuit->mnparm(3,  "s mpv A", vstart[3] ,  step[3], 0, 0, ierflg);
  fMinuit->mnparm(4,  "s mpv B", vstart[4] ,  step[4], 0, 0, ierflg);
  fMinuit->mnparm(5,  "s mpv C", vstart[5] ,  step[5], 0, 0, ierflg);
  fMinuit->mnparm(6,  "s mpv D", vstart[6] ,  step[6], 0, 0, ierflg);
  fMinuit->mnparm(7,  "s mpv R", vstart[7] ,  step[7], 0, 0, ierflg);
  fMinuit->mnparm(8,  "s gw A" , vstart[8] ,  step[8], 0, 0, ierflg);
  fMinuit->mnparm(9,  "s gw B" , vstart[9] ,  step[9], 0, 0, ierflg);
  fMinuit->mnparm(10, "s gw C" , vstart[10], step[10], 0, 0, ierflg);
  fMinuit->mnparm(11, "b lw A" , vstart[11], step[11], 0, 0, ierflg);
  fMinuit->mnparm(12, "b lw B" , vstart[12], step[12], 0, 0, ierflg);
  fMinuit->mnparm(13, "b lw C" , vstart[13], step[13], 0, 0, ierflg);
  fMinuit->mnparm(14, "b shift", vstart[14], step[14], 0, 0, ierflg);
  fMinuit->mnparm(15, "b gw A" , vstart[15], step[15], 0, 0, ierflg);
  fMinuit->mnparm(16, "b gw B" , vstart[16], step[16], 0, 0, ierflg);
  fMinuit->mnparm(17, "b gw C" , vstart[17], step[17], 0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("s Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 0, ierflg);//0.82, 0.90, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  slwA  , slwA_error);
  fMinuit->GetParameter(1,  slwB  , slwB_error);
  fMinuit->GetParameter(2,  slwC  , slwC_error);
  fMinuit->GetParameter(3,  smpvA , smpvA_error);
  fMinuit->GetParameter(4,  smpvB , smpvB_error);
  fMinuit->GetParameter(5,  smpvC , smpvC_error);
  fMinuit->GetParameter(6,  smpvD , smpvD_error);
  fMinuit->GetParameter(7,  smpvR , smpvR_error);
  fMinuit->GetParameter(8,  sgwA  , sgwA_error);
  fMinuit->GetParameter(9,  sgwB  , sgwB_error);
  fMinuit->GetParameter(10, sgwC  , sgwC_error);
  fMinuit->GetParameter(11, blwA  , blwA_error);
  fMinuit->GetParameter(12, blwB  , blwB_error);
  fMinuit->GetParameter(13, blwC  , blwC_error);
  fMinuit->GetParameter(14, bshift, bshift_error);
  fMinuit->GetParameter(15, bgwA  , bgwA_error);
  fMinuit->GetParameter(16, bgwB  , bgwB_error);
  fMinuit->GetParameter(17, bgwC  , bgwC_error);
  std::vector<std::pair<double,double>> sn;
  std::vector<std::pair<double,double>> bn;
  for(int ipar = 0; ipar < NBINS; ipar++){
    double snorm, snorm_error;//, bnorm, bnorm_error;
    fMinuit->GetParameter(ipar+CPAR, snorm, snorm_error);
    sn.push_back(std::make_pair(snorm, snorm_error));
    //fMinuit->GetParameter(ipar+CPAR+NBINS, bnorm, bnorm_error);
    bn.push_back(std::make_pair(1-snorm, snorm_error));
  }

  GetSignal()->SetClwA(std::make_pair(slwA,slwA_error));
  GetSignal()->SetClwB(std::make_pair(slwB,slwB_error));
  GetSignal()->SetClwC(std::make_pair(slwC,slwC_error));
  GetSignal()->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  GetSignal()->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  GetSignal()->SetCmpvC(std::make_pair(smpvC,smpvR_error));
  GetSignal()->SetCmpvD(std::make_pair(smpvD,smpvC_error));
  GetSignal()->SetCmpvR(std::make_pair(smpvR,smpvD_error));
  GetSignal()->SetCgwA(std::make_pair(sgwA,sgwA_error));
  GetSignal()->SetCgwB(std::make_pair(sgwB,sgwB_error));
  GetSignal()->SetCgwC(std::make_pair(sgwC,sgwC_error));
  GetSignal()->SetCnorm(sn);

  GetBackground()->SetClwA(std::make_pair(blwA,blwA_error));
  GetBackground()->SetClwB(std::make_pair(blwB,blwB_error));
  GetBackground()->SetClwC(std::make_pair(blwC,blwC_error));
  GetBackground()->SetCmpvA(GetSignal()->GetCmpvA());
  GetBackground()->SetCmpvB(GetSignal()->GetCmpvB());
  GetBackground()->SetCmpvC(GetSignal()->GetCmpvC());
  GetBackground()->SetCmpvD(GetSignal()->GetCmpvD());
  GetBackground()->SetCmpvR(GetSignal()->GetCmpvR());
  GetBackground()->SetCshift(std::make_pair(bshift,bshift_error));
  GetBackground()->SetCgwA(std::make_pair(bgwA,bgwA_error));
  GetBackground()->SetCgwB(std::make_pair(bgwB,bgwB_error));
  GetBackground()->SetCgwC(std::make_pair(bgwC,bgwC_error));
  GetBackground()->SetCnorm(bn);
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundShift(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  const int CPAR = 18;
  const int NBINS = CSMinuit->fh.size();
  
  //set each parameter
  double slwA  = par[0];
  double slwB  = par[1];
  double slwC  = par[2];
  double smpvA = par[3];
  double smpvB = par[4];
  double smpvC = par[5];
  double smpvD = par[6];
  double smpvR = par[7];
  double sgwA  = par[8];
  double sgwB  = par[9];
  double sgwC  = par[10];
  double blwA  = par[11];
  double blwB  = par[12];
  double blwC  = par[13];
  double shift = par[14];
  double bgwA  = par[15];
  double bgwB  = par[16];
  double bgwC  = par[17];
  
  //define function to fit to each histogram
  double slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < NBINS; ibin++){
    //get RR depending parameters
    double slw_par[]  = {slwA,slwB,slwC,0};
    double smpv_par[] = {smpvA,smpvB,smpvC,smpvD,smpvR,0};
    double sgw_par[]  = {sgwA,sgwB,sgwC,0};
    double blw_par[]  = {blwA,blwB,blwC,0};
    double bmpv_par[] = {smpvA,smpvB,smpvC,smpvD,smpvR,shift};
    double bgw_par[]  = {bgwA,bgwB,bgwC,0};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    
    slw   = CoherentFitUtils::ABCParametrization(RR_par,slw_par);
    smpv  = CoherentFitUtils::ABCDRParametrization(RR_par,smpv_par);
    sgw   = CoherentFitUtils::ABCParametrization(RR_par,sgw_par);
    snorm = par[CPAR+ibin];
    blw   = CoherentFitUtils::ABCParametrization(RR_par,blw_par);
    bmpv  = CoherentFitUtils::ABCDRParametrization(RR_par,bmpv_par);
    bgw   = CoherentFitUtils::ABCParametrization(RR_par,bgw_par);
    bnorm = 1-par[CPAR+ibin];
    
    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//shift and withs as cuadratic sums from signal widths
//********************************************************************
void CoherentSample::CoherentFitBackground3Par(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 3;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackground3Par);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {0.1,  -1, 0.1};
  Double_t step[NPAR]   = {0.1, 0.1, 0.1};
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fInorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
  }
  fMinuit->mnparm(0,  "lw_shift" ,  vstart[0],  step[0],  0, 0, ierflg);
  fMinuit->mnparm(1,  "mpv_shift",  vstart[1],  step[1],  0, 0, ierflg);
  fMinuit->mnparm(2,  "gw_shift" ,  vstart[2],  step[2],  0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fClwQa.first,  fClwQa.second);
  fMinuit->GetParameter(1,  fCshift.first, fCshift.second);
  fMinuit->GetParameter(2,  fCgwQa.first,  fCgwQa.second);

  for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCnorm.push_back(std::make_pair(norm, norm_error));
  }

  fClwA = fTrueSignal->GetClwA();
  fClwB = fTrueSignal->GetClwB();
  fClwC = fTrueSignal->GetClwC();
  fCmpvA = fTrueSignal->GetCmpvA();
  fCmpvB = fTrueSignal->GetCmpvB();
  fCmpvC = fTrueSignal->GetCmpvC();
  fCmpvD = fTrueSignal->GetCmpvD();
  fCmpvR = fTrueSignal->GetCmpvR();
  fCgwA = fTrueSignal->GetCgwA();
  fCgwB = fTrueSignal->GetCgwB();
  fCgwC = fTrueSignal->GetCgwC();
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackground3Par(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 3;
  
  //set each parameter
  double lwA   = CSMinuit->GetTrueSignal()->GetClwA().first;
  double lwB   = CSMinuit->GetTrueSignal()->GetClwB().first;
  double lwC   = CSMinuit->GetTrueSignal()->GetClwC().first;
  double mpvA  = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  double mpvB  = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  double mpvC  = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  double mpvD  = CSMinuit->GetTrueSignal()->GetCmpvD().first;
  double mpvR  = CSMinuit->GetTrueSignal()->GetCmpvR().first;
  double gwA   = CSMinuit->GetTrueSignal()->GetCgwA().first;
  double gwB   = CSMinuit->GetTrueSignal()->GetCgwB().first;
  double gwC   = CSMinuit->GetTrueSignal()->GetCgwC().first;
  double lw_s  = par[0];
  double mpv_s = par[1];
  double gw_s  = par[2];
  
  //define function to fit to each histogram
  double lw,mpv,norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double lw_par[]   = {lwA,lwB,lwC,lw_s};
    double mpv_par[]  = {mpvA,mpvB,mpvC,mpvD,mpvR,mpv_s};
    double gw_par[]   = {gwA,gwB,gwC,gw_s};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    lw   = CoherentFitUtils::ABCParametrization(RR_par,lw_par);
    mpv  = CoherentFitUtils::ABCDRParametrization(RR_par,mpv_par);
    gw   = CoherentFitUtils::ABCParametrization(RR_par,gw_par);
    norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//shift and widths as quadratic ABC addition to signal width
//shift and withs as cuadratic sums from signal widths
//********************************************************************
void CoherentSample::CoherentFitBackgroundQuadraticWidths(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 7+1;
  const int NBINS = fh.size();
  const int NPAR = CPAR;// + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundQuadraticWidths);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {2.6,  4.60, -1.98, 0.063, 36.9, -2.01, 0.09};
  Double_t step[NPAR]   = {0.1,   0.1,   0.1,   0.1,  0.1,  0.01, 0.01};
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fInorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
    }*/
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = mean_norm / fInorm.size();
  fMinuit->mnparm(0,  "shift", vstart[0], step[0], 0,  10, ierflg);
  fMinuit->mnparm(1,  "lwA"  , vstart[1], step[1], 0,  10, ierflg);
  fMinuit->mnparm(2,  "lwB"  , vstart[2], step[2], -10, 0, ierflg);
  fMinuit->mnparm(3,  "lwC"  , vstart[3], step[3], 0,   1, ierflg);
  fMinuit->mnparm(4,  "gwA"  , vstart[4], step[4], 0,  90, ierflg);
  fMinuit->mnparm(5,  "gwB"  , vstart[5], step[5], -10, 0, ierflg);
  fMinuit->mnparm(6,  "gwC"  , vstart[6], step[6], 0,   1, ierflg);
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
		    }*/
  fMinuit->mnparm(7,  "norm"  , mean_norm, 0.01, 0,   1, ierflg);
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0, fCshift.first, fCshift.second);
  fMinuit->GetParameter(1, fClwA.first,   fClwA.second);
  fMinuit->GetParameter(2, fClwB.first,   fClwB.second);
  fMinuit->GetParameter(3, fClwC.first,   fClwC.second);
  fMinuit->GetParameter(4, fCgwA.first,   fCgwA.second);
  fMinuit->GetParameter(5, fCgwB.first,   fCgwB.second);
  fMinuit->GetParameter(6, fCgwC.first,   fCgwC.second);

  /*for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCnorm.push_back(std::make_pair(norm, norm_error));
    }*/

  double norm, norm_error;
  fMinuit->GetParameter(7, norm, norm_error);
  fCnorm.push_back(std::make_pair(norm, norm_error));
  
  fCmpvA = fTrueSignal->GetCmpvA();
  fCmpvB = fTrueSignal->GetCmpvB();
  fCmpvC = fTrueSignal->GetCmpvC();
  fCmpvD = fTrueSignal->GetCmpvD();
  fCmpvR = fTrueSignal->GetCmpvR();
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundQuadraticWidths(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 7;
  
  //set each parameter
  double shift = par[0];
  double blwA  = par[1];
  double blwB  = par[2];
  double blwC  = par[3];
  double bgwA  = par[4];
  double bgwB  = par[5];
  double bgwC  = par[6];
  double norm  = par[7];
  double slwA  = CSMinuit->GetTrueSignal()->GetClwA().first;
  double slwB  = CSMinuit->GetTrueSignal()->GetClwB().first;
  double slwC  = CSMinuit->GetTrueSignal()->GetClwC().first;
  double smpvA = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  double smpvB = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  double smpvC = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  double smpvD = CSMinuit->GetTrueSignal()->GetCmpvD().first;
  double smpvR = CSMinuit->GetTrueSignal()->GetCmpvR().first;
  double sgwA  = CSMinuit->GetTrueSignal()->GetCgwA().first;
  double sgwB  = CSMinuit->GetTrueSignal()->GetCgwB().first;
  double sgwC  = CSMinuit->GetTrueSignal()->GetCgwC().first;
  
  //define function to fit to each histogram
  double lw,mpv,gw;//norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double lw_par[]   = {slwA,slwB,slwC,blwA,blwB,blwC};
    double mpv_par[]  = {smpvA,smpvB,smpvC,smpvD,smpvR,shift};
    double gw_par[]   = {sgwA,sgwB,sgwC,bgwA,bgwB,bgwC};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    lw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,lw_par);
    mpv  = CoherentFitUtils::ABCDRParametrization(RR_par,mpv_par);
    gw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,gw_par);
    //norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundQuadraticWidths(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 18+1;
  const int NBINS = fh.size();
  const int NPAR = CPAR;// + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundQuadraticWidths);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,smpvR,sgwA,sgwB,sgwC,blwA,blwB,blwC,bshift,bgwA,bgwB,bgwC,slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,smpvR_error,sgwA_error,sgwB_error,sgwC_error,blwA_error,blwB_error,blwC_error,bshift_error,bgwA_error,bgwB_error,bgwC_error;

  slwA   = GetSignal()->GetClwA().first;         //  4.70175e+00;
  slwB   = GetSignal()->GetClwB().first;	 // -3.16696e+00;
  slwC   = GetSignal()->GetClwC().first;	 //  5.16751e-02;
  smpvA  = GetSignal()->GetCmpvA().first;	 // -5.64484e-01;
  smpvB  = GetSignal()->GetCmpvB().first;	 // -1.11127e+00;
  smpvC  = GetSignal()->GetCmpvC().first;	 // -1.06019e+02;
  smpvD  = GetSignal()->GetCmpvD().first;	 //  1.98386e+00;
  smpvR  = GetSignal()->GetCmpvR().first;	 //  6.08834e+00;
  sgwA   = GetSignal()->GetCgwA().first;	 //  2.43725e+00;
  sgwB   = GetSignal()->GetCgwB().first;	 // -1.07204e+00;
  sgwC   = GetSignal()->GetCgwC().first;	 //  1.46027e-01;
  blwA   = GetBackground()->GetClwA().first;	 //  3.01259e-05;
  blwB   = GetBackground()->GetClwB().first;	 // -4.11839e+00;
  blwC   = GetBackground()->GetClwC().first;	 //  3.74797e-02;
  bshift = GetBackground()->GetCshift().first;   //  4.97744e+00;
  bgwA   = GetBackground()->GetCgwA().first;	 //  9.34709e-01;
  bgwB   = GetBackground()->GetCgwB().first;	 // -3.94219e-01;
  bgwC   = GetBackground()->GetCgwC().first;	 //  1.98209e-07;
  Double_t vstart[NPAR] = {slwA, slwB, slwC, smpvA, smpvB, smpvC, smpvD, smpvR, sgwA, sgwB, sgwC,
			   blwA, blwB, blwC, bshift,                            bgwA, bgwB, bgwC};
  Double_t step[NPAR]   = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  std::vector<std::pair<double,double>> snorm = GetSignal()->GetCnorm();
  std::vector<std::pair<double,double>> bnorm = GetBackground()->GetCnorm();
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = (snorm[ipar].first+(1-bnorm[ipar].first))/2;//dummy[ipar];
    step[ipar+CPAR]   = 0.01;
    }*/
  vstart[18] = (snorm[0].first+(1-bnorm[0].first))/2;//dummy[ipar];
  step[18]   = 0.01;
    
  fMinuit->mnparm(0,  "s lw A" , vstart[0] ,  step[0],   0,  0, ierflg);
  fMinuit->mnparm(1,  "s lw B" , vstart[1] ,  step[1],   0,  0, ierflg);
  fMinuit->mnparm(2,  "s lw C" , vstart[2] ,  step[2],   0,  0, ierflg);
  fMinuit->mnparm(3,  "s mpv A", vstart[3] ,  step[3],   0,  0, ierflg);
  fMinuit->mnparm(4,  "s mpv B", vstart[4] ,  step[4],   0,  0, ierflg);
  fMinuit->mnparm(5,  "s mpv C", vstart[5] ,  step[5],   0,  0, ierflg);
  fMinuit->mnparm(6,  "s mpv D", vstart[6] ,  step[6],   0,  0, ierflg);
  fMinuit->mnparm(7,  "s mpv R", vstart[7] ,  step[7],   0,  0, ierflg);
  fMinuit->mnparm(8,  "s gw A" , vstart[8] ,  step[8],   0,  0, ierflg);
  fMinuit->mnparm(9,  "s gw B" , vstart[9] ,  step[9],   0,  0, ierflg);
  fMinuit->mnparm(10, "s gw C" , vstart[10], step[10],   0,  0, ierflg);
  fMinuit->mnparm(11, "b lw A" , vstart[11], step[11],   0, 10, ierflg);
  fMinuit->mnparm(12, "b lw B" , vstart[12], step[12], -10,  0, ierflg);
  fMinuit->mnparm(13, "b lw C" , vstart[13], step[13],   0, 10, ierflg);
  fMinuit->mnparm(14, "b shift", vstart[14], step[14],   0,  0, ierflg);//vstart[14]*0.95, vstart[14]*1.05, ierflg);
  fMinuit->mnparm(15, "b gw A" , vstart[15], step[15],   0, 90, ierflg);
  fMinuit->mnparm(16, "b gw B" , vstart[16], step[16], -10,  0, ierflg);
  fMinuit->mnparm(17, "b gw C" , vstart[17], step[17],   0, 10, ierflg);
  //fMinuit->FixParameter(0);
  //fMinuit->FixParameter(1);
  //fMinuit->FixParameter(2);
  //fMinuit->FixParameter(3);
  //fMinuit->FixParameter(4);
  //fMinuit->FixParameter(5);
  //fMinuit->FixParameter(6);
  //fMinuit->FixParameter(7);
  //fMinuit->FixParameter(8);
  //fMinuit->FixParameter(9);
  //fMinuit->FixParameter(10);
  //fMinuit->FixParameter(11);
  //fMinuit->FixParameter(12);
  fMinuit->FixParameter(13);
  //fMinuit->FixParameter(14);
  //fMinuit->FixParameter(15);
  //fMinuit->FixParameter(16);
  fMinuit->FixParameter(17);
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    if(ipar == 0 || ipar == NBINS-1)
      fMinuit->mnparm(ipar+CPAR, ("s Ar "+ssi.str()+"").c_str(),
		      vstart[ipar+CPAR], step[ipar+CPAR], 0.975*vstart[ipar+CPAR], 1.025*vstart[ipar+CPAR], ierflg); //0.85, 0.88, ierflg);//
    else
      fMinuit->mnparm(ipar+CPAR, ("s Ar "+ssi.str()+"").c_str(),
		      vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg); //0.85, 0.88, ierflg);//
    //fMinuit->FixParameter(ipar+CPAR);
    }*/
  //fMinuit->mnparm(18, "norm" , vstart[18], step[18], 0.95*vstart[18], 1.05*vstart[18], ierflg);
  fMinuit->mnparm(18, "norm" , vstart[18], step[18], 0.9*vstart[18], 1.1*vstart[18], ierflg);
  
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  slwA  , slwA_error);
  fMinuit->GetParameter(1,  slwB  , slwB_error);
  fMinuit->GetParameter(2,  slwC  , slwC_error);
  fMinuit->GetParameter(3,  smpvA , smpvA_error);
  fMinuit->GetParameter(4,  smpvB , smpvB_error);
  fMinuit->GetParameter(5,  smpvC , smpvC_error);
  fMinuit->GetParameter(6,  smpvD , smpvD_error);
  fMinuit->GetParameter(7,  smpvR , smpvR_error);
  fMinuit->GetParameter(8,  sgwA  , sgwA_error);
  fMinuit->GetParameter(9,  sgwB  , sgwB_error);
  fMinuit->GetParameter(10, sgwC  , sgwC_error);
  fMinuit->GetParameter(11, blwA  , blwA_error);
  fMinuit->GetParameter(12, blwB  , blwB_error);
  fMinuit->GetParameter(13, blwC  , blwC_error);
  fMinuit->GetParameter(14, bshift, bshift_error);
  fMinuit->GetParameter(15, bgwA  , bgwA_error);
  fMinuit->GetParameter(16, bgwB  , bgwB_error);
  fMinuit->GetParameter(17, bgwC  , bgwC_error);
  std::vector<std::pair<double,double>> sn;
  std::vector<std::pair<double,double>> bn;
  /*for(int ipar = 0; ipar < NBINS; ipar++){
    double snorm, snorm_error;//, bnorm, bnorm_error;
    fMinuit->GetParameter(ipar+CPAR, snorm, snorm_error);
    sn.push_back(std::make_pair(snorm, snorm_error));
    //fMinuit->GetParameter(ipar+CPAR+NBINS, bnorm, bnorm_error);
    bn.push_back(std::make_pair(1-snorm, snorm_error));
    }*/
  double norm, norm_error;//, bnorm, bnorm_error;
  fMinuit->GetParameter(18, norm, norm_error);
  sn.push_back(std::make_pair(norm, norm_error));
  bn.push_back(std::make_pair(1-norm, norm_error));

  GetSignal()->SetClwA(std::make_pair(slwA,slwA_error));
  GetSignal()->SetClwB(std::make_pair(slwB,slwB_error));
  GetSignal()->SetClwC(std::make_pair(slwC,slwC_error));
  GetSignal()->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  GetSignal()->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  GetSignal()->SetCmpvC(std::make_pair(smpvC,smpvR_error));
  GetSignal()->SetCmpvD(std::make_pair(smpvD,smpvC_error));
  GetSignal()->SetCmpvR(std::make_pair(smpvR,smpvD_error));
  GetSignal()->SetCgwA(std::make_pair(sgwA,sgwA_error));
  GetSignal()->SetCgwB(std::make_pair(sgwB,sgwB_error));
  GetSignal()->SetCgwC(std::make_pair(sgwC,sgwC_error));
  GetSignal()->SetCnorm(sn);

  GetBackground()->SetClwA(std::make_pair(blwA,blwA_error));
  GetBackground()->SetClwB(std::make_pair(blwB,blwB_error));
  GetBackground()->SetClwC(std::make_pair(blwC,blwC_error));
  GetBackground()->SetCmpvA(GetSignal()->GetCmpvA());
  GetBackground()->SetCmpvB(GetSignal()->GetCmpvB());
  GetBackground()->SetCmpvC(GetSignal()->GetCmpvC());
  GetBackground()->SetCmpvD(GetSignal()->GetCmpvD());
  GetBackground()->SetCmpvR(GetSignal()->GetCmpvR());
  GetBackground()->SetCshift(std::make_pair(bshift,bshift_error));
  GetBackground()->SetCgwA(std::make_pair(bgwA,bgwA_error));
  GetBackground()->SetCgwB(std::make_pair(bgwB,bgwB_error));
  GetBackground()->SetCgwC(std::make_pair(bgwC,bgwC_error));
  GetBackground()->SetCnorm(bn);
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundQuadraticWidths(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  const int CPAR = 18+1;
  const int NBINS = CSMinuit->fh.size();
  
  //set each parameter
  double slwA  = par[0];
  double slwB  = par[1];
  double slwC  = par[2];
  double smpvA = par[3];
  double smpvB = par[4];
  double smpvC = par[5];
  double smpvD = par[6];
  double smpvR = par[7];
  double sgwA  = par[8];
  double sgwB  = par[9];
  double sgwC  = par[10];
  double blwA  = par[11];
  double blwB  = par[12];
  double blwC  = par[13];
  double shift = par[14];
  double bgwA  = par[15];
  double bgwB  = par[16];
  double bgwC  = par[17];
  double norm  = par[18];
  
  //define function to fit to each histogram
  double slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw;
  //std::vector<double> vnorm;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < NBINS; ibin++){
    //get RR depending parameters
    double slw_par[]  = {slwA,slwB,slwC,0};
    double smpv_par[] = {smpvA,smpvB,smpvC,smpvD,smpvR,0};
    double sgw_par[]  = {sgwA,sgwB,sgwC,0};
    double blw_par[]  = {slwA,slwB,slwC,blwA,blwB,blwC};
    double bmpv_par[] = {smpvA,smpvB,smpvC,smpvD,smpvR,shift};
    double bgw_par[]  = {sgwA,sgwB,sgwC,bgwA,bgwB,bgwC};
    double RR_par[]   = {CSMinuit->fRR[ibin].first};
    
    slw   = CoherentFitUtils::ABCParametrization(RR_par,slw_par);
    smpv  = CoherentFitUtils::ABCDRParametrization(RR_par,smpv_par);
    sgw   = CoherentFitUtils::ABCParametrization(RR_par,sgw_par);
    snorm = norm;//par[CPAR+ibin];
    blw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,blw_par);
    bmpv  = CoherentFitUtils::ABCDRParametrization(RR_par,bmpv_par);
    bgw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,bgw_par);
    bnorm = 1-norm;//par[CPAR+ibin];

    //vnorm.push_back(snorm);
    
    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }

  //double regularization = CoherentFitUtils::NormRegularization(vnorm,CSMinuit->GetIntegralVector());
  //regularization
  //Likelihood = Likelihood - CoherentFitUtils::NormRegularization(vnorm);
  //std::cout << Likelihood << " " << regularization << std::endl;
  //Likelihood = Likelihood - regularization;
  f = -Likelihood;
}

//shift on data and lw
/*
//********************************************************************
void CoherentSample::CoherentFitBackground(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 4;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackground);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {0, fCbgwA.first, fCbgwB.first, fCbgwC.first};
  Double_t step[NPAR]   = {1,    0.1,     0.1,   0.1};
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fIbnorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
  }
  fMinuit->mnparm(0,  "shift",  vstart[0],  step[0],  0, 0, ierflg);
  fMinuit->mnparm(1,  "gwA"  ,  vstart[1],  step[1],  0, 0, ierflg);
  fMinuit->mnparm(2,  "gwB"  ,  vstart[2],  step[2],  0, 0, ierflg);
  fMinuit->mnparm(3,  "gwC"  ,  vstart[3],  step[3],  0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fCbshift.first,  fCbshift.second);
  fMinuit->GetParameter(1,  fCbgwA.first  ,  fCbgwA.second);
  fMinuit->GetParameter(2,  fCbgwB.first  ,  fCbgwB.second);
  fMinuit->GetParameter(3,  fCbgwC.first  ,  fCbgwC.second);
  for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCbnorm.push_back(std::make_pair(norm, norm_error));
  }

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 4;
  
  //set each parameter
  double shift = par[0];
  double lwA   = CSMinuit->GetAssociatedSignal()->GetCslwA().first;
  double lwB   = CSMinuit->GetAssociatedSignal()->GetCslwB().first;
  double lwC   = CSMinuit->GetAssociatedSignal()->GetCslwC().first;
  double mpvA  = CSMinuit->GetAssociatedSignal()->GetCsmpvA().first;
  double mpvB  = CSMinuit->GetAssociatedSignal()->GetCsmpvB().first;
  double mpvC  = CSMinuit->GetAssociatedSignal()->GetCsmpvC().first;
  double mpvD  = CSMinuit->GetAssociatedSignal()->GetCsmpvD().first;
  double mpvR  = CSMinuit->GetAssociatedSignal()->GetCsmpvR().first;
  double gwA   = par[1];
  double gwB   = par[2];
  double gwC   = par[3];
  
  //define function to fit to each histogram
  double lw,mpv,norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    lw   = CoherentFitUtils::ABCParametrization(CSMinuit->fRR[ibin].first+shift,lwA,lwB,lwC);
    mpv  = CoherentFitUtils::ABCDRParametrization(CSMinuit->fRR[ibin].first+shift,mpvA,mpvB,mpvC,mpvD,mpvR);
    gw   = CoherentFitUtils::ABCParametrization(CSMinuit->fRR[ibin].first,gwA,gwB,gwC);
    norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}
 */

//shift everywhere
/*//********************************************************************
void CoherentSample::CoherentFitBackground(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 1;
  const int NBINS = fh.size();
  const int NPAR = CPAR + NBINS;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackground);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {0};
  Double_t step[NPAR]   = {10};
  for(int ipar = 0; ipar < NBINS; ipar++){
    vstart[ipar+CPAR] = fIbnorm[ipar].first;
    step[ipar+CPAR]   = 0.1;
  }
  fMinuit->mnparm(0,  "shift",  vstart[0],  step[0],  0, 0, ierflg);
  for(int ipar = 0; ipar < NBINS; ipar++){
    std::stringstream ssi;
    ssi << ipar;
    fMinuit->mnparm(ipar+CPAR, ("Ar "+ssi.str()+"").c_str(),
		    vstart[ipar+CPAR], step[ipar+CPAR], 0, 1, ierflg);
  }
  
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fCbshift.first,  fCbshift.second);
  for(int ipar = 0; ipar < NBINS; ipar++){
    double norm, norm_error;
    fMinuit->GetParameter(ipar+CPAR, norm, norm_error);
    fCbnorm.push_back(std::make_pair(norm, norm_error));
  }

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  const int CPAR = 1;
  
  //set each parameter
  double lwA   = CSMinuit->GetAssociatedSignal()->GetCslwA().first;
  double lwB   = CSMinuit->GetAssociatedSignal()->GetCslwB().first;
  double lwC   = CSMinuit->GetAssociatedSignal()->GetCslwC().first;
  double mpvA  = CSMinuit->GetAssociatedSignal()->GetCsmpvA().first;
  double mpvB  = CSMinuit->GetAssociatedSignal()->GetCsmpvB().first;
  double mpvC  = CSMinuit->GetAssociatedSignal()->GetCsmpvC().first;
  double mpvD  = CSMinuit->GetAssociatedSignal()->GetCsmpvD().first;
  double mpvR  = CSMinuit->GetAssociatedSignal()->GetCsmpvR().first;
  double gwA   = CSMinuit->GetAssociatedSignal()->GetCsgwA().first;
  double gwB   = CSMinuit->GetAssociatedSignal()->GetCsgwB().first;
  double gwC   = CSMinuit->GetAssociatedSignal()->GetCsgwC().first;
  double shift = par[0];
  
  //define function to fit to each histogram
  double lw,mpv,norm,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    lw   = CoherentFitUtils::ABCParametrization(CSMinuit->fRR[ibin].first+shift,lwA,lwB,lwC);
    mpv  = CoherentFitUtils::ABCDRParametrization(CSMinuit->fRR[ibin].first+shift,mpvA,mpvB,mpvC,mpvD,mpvR);
    gw   = CoherentFitUtils::ABCParametrization(CSMinuit->fRR[ibin].first+shift,gwA,gwB,gwC);
    norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}*/

//SIGNAL PLUS BACKGROUND MODELS
//lw,gw free, mpv shift.

