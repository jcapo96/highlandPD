#include "CoherentSample.hxx"
#include "CoherentFitUtils.hxx"

#include "TVector.h"
#include "TParameter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "Math/Math.h"
#include "Math/ChebyshevPol.h"
#include "Math/SpecFuncMathMore.h"

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

  fCnorm = std::make_pair(0,0);
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

  if(fType != SampleTypeEnum::kSignalPlusBackground){
    if(fClwFit)fClwFit->Write();
    if(fCmpvFit)fCmpvFit->Write();
    if(fCgwFit)fCgwFit->Write();
  }
  
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

  rfile->GetObject("Cnorm",p);
  if(p)fCnorm = *p;
  
  rfile->Close();
}

//********************************************************************
void CoherentSample::ResetAllButHistograms(){
//********************************************************************

  for(int i = 0; i < (int)fIFit.size(); i++)
    delete fIFit[i];
  for(int i = 0; i < (int)fCFit.size(); i++)
    delete fCFit[i];
  
  fIFit.clear();
  fCFit.clear();

  delete fClwFit;  
  delete fCmpvFit; 
  delete fCgwFit;  
  fClwFit  = NULL;
  fCmpvFit = NULL;
  fCgwFit  = NULL;
  
  fIlw.clear();
  fImpv.clear();
  fInorm.clear();
  fIgw.clear();
  fIfw.clear();
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
void CoherentSample::ChangeHistogramsToVariableBinning(const int nmin){
//********************************************************************
  
  if(fh.empty()){
    std::cout << "no histogram stored yet" << std::endl;
    return;
  }

  for(int ihist = 0; ihist < (int)fh.size(); ihist++)
    fh[ihist] = CoherentFitUtils::ChangeHistogramToVariableBinning(fh[ihist],fIntegral[ihist],nmin);
}

//********************************************************************
void CoherentSample::StoreCoherentFits(){
//********************************************************************

  if(fType != SampleTypeEnum::kSignalPlusBackground){

    if(fType == SampleTypeEnum::kSignal || fType == SampleTypeEnum::kTrueSignal){
      fClwFit = new TF1("lw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fClwFit->SetParameters(fClwA.first, fClwB.first, fClwC.first,0,0,0);
      
      fCmpvFit = new TF1("mpv_par",CoherentFitUtils::ABCDRParametrization,0,60,6);
      fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,0);
      
      fCgwFit = new TF1("gw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fCgwFit->SetParameters(fCgwA.first, fCgwB.first, fCgwC.first, 0,0,0);
    }

    else{
      CoherentSample* satb = NULL;
      if(fType == SampleTypeEnum::kTrueBackground)satb = GetTrueSignal();
      if(fType == SampleTypeEnum::kBackground)    satb = GetSignal();
      
      fClwFit = new TF1("lw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fClwFit->SetParameters(satb->GetClwA().first,satb->GetClwB().first,satb->GetClwC().first, fClwA.first, fClwB.first, fClwC.first);
      
      fCmpvFit = new TF1("mpv_par",CoherentFitUtils::ABCDRParametrization,0,60,6);
      fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,fCshift.first);
      
      fCgwFit = new TF1("gw_par" ,CoherentFitUtils::QuadraticABCParametrization,0,60,6);
      fCgwFit->SetParameters(satb->GetCgwA().first,satb->GetCgwB().first,satb->GetCgwC().first, fCgwA.first, fCgwB.first, fCgwC.first);
    }

    // if(fType == SampleTypeEnum::kSignal || fType == SampleTypeEnum::kTrueSignal){
    //   fClwFit = new TF1("lw_par","ROOT::Math::Chebyshev2((x-1)/(x+1),[0],[1],[2])",0,30);
    //   fClwFit->SetParameters(fClwA.first,fClwB.first,fClwC.first);
    //   fCmpvFit = new TF1("mpv_par","ROOT::Math::Chebyshev4((x-1)/(x+1),[0],[1],[2],[3],[4])",0,30);
    //   fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first);
    //   fCgwFit = new TF1("gw_par","ROOT::Math::Chebyshev2((x-1)/(x+1),[0],[1],[2])",0,30);
    //   fCgwFit->SetParameters(fCgwA.first,fCgwB.first,fCgwC.first);
    // }
    // else if(fType == SampleTypeEnum::kBackground || fType == SampleTypeEnum::kTrueBackground){
    //   fClwFit = new TF1("lw_par","ROOT::Math::Chebyshev2((x+[3]-1)/(x+[3]+1),[0],[1],[2])",0,30);
    //   fClwFit->SetParameters(satb->GetClwA().first,satb->GetClwB().first,satb->GetClwC().first,fClwA.first);
    //   //fClwFit->SetParameters(satb->GetClwA().first,satb->GetClwB().first,satb->GetClwC().first,-fCshift.first);
    //   fCmpvFit = new TF1("mpv_par","ROOT::Math::Chebyshev4((x+[5]-1)/(x+[5]+1),[0],[1],[2],[3],[4])",0,30);
    //   fCmpvFit->SetParameters(satb->GetCmpvA().first,satb->GetCmpvB().first,satb->GetCmpvC().first,satb->GetCmpvD().first,satb->GetCmpvR().first,fCshift.first);
    //   fCgwFit = new TF1("gw_par","ROOT::Math::Chebyshev2((x+[3]-1)/(x+[3]+1),[0],[1],[2])",0,30);
    //   fCgwFit->SetParameters(satb->GetCgwA().first,satb->GetCgwB().first,satb->GetCgwC().first,fClwA.first);
    //   //fCgwFit->SetParameters(satb->GetCgwA().first,satb->GetCgwB().first,satb->GetCgwC().first,-fCshift.first);
    //   //fCgwFit->SetParameters(fCgwA.first,fCgwB.first,fCgwC.first,0);

    //   std::cout << "MIGUE " << fClwA.first << " " << fCshift.first << std::endl;
    // }
    
    //store invidual histogram fits
    for(int i = 0; i < (int)fRR.size(); i++){
      std::stringstream ssi;
      ssi << i;
      TF1* f = new TF1(("CF_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);

      double lw,mpv,norm,gw;
      
      lw   = fClwFit->Eval(fRR[i].first);
      mpv  = fCmpvFit->Eval(fRR[i].first);
      norm = fCnorm.first;
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
      snorm = GetSignal()->GetCnorm().first;
      sgw   = GetSignal()->GetCgwFit()->Eval(fRR[i].first);
      blw   = GetBackground()->GetClwFit()->Eval(fRR[i].first);
      bmpv  = GetBackground()->GetCmpvFit()->Eval(fRR[i].first);
      bnorm = GetBackground()->GetCnorm().first;
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
  this->fCnorm.first    = sample->GetCnorm().first;
  
  //apply desired variations
  if(apply_all_var || apply_only_lw_var){
    //this->fClwA.first     = r->Gaus(sample->GetClwA().first,sample->GetClwA().first*sigma);
    this->fClwB.first     = r->Gaus(sample->GetClwB().first,sample->GetClwB().first*sigma);
    //this->fClwC.first     = r->Gaus(sample->GetClwC().first,sample->GetClwC().first*sigma);
  }
  if(apply_all_var || apply_only_mpv_var){
    //this->fCmpvA.first    = r->Gaus(sample->GetCmpvA().first,sample->GetCmpvA().first*sigma);
    //this->fCmpvB.first    = r->Gaus(sample->GetCmpvB().first,sample->GetCmpvB().first*sigma);
    //this->fCmpvC.first    = r->Gaus(sample->GetCmpvC().first,sample->GetCmpvC().first*sigma);
    //this->fCmpvD.first    = r->Gaus(sample->GetCmpvD().first,sample->GetCmpvD().first*sigma);
    this->fCmpvR.first    = r->Gaus(sample->GetCmpvR().first,sample->GetCmpvR().first*sigma);
  }
  if(apply_all_var || apply_only_gw_var){
    //this->fCgwA.first     = r->Gaus(sample->GetCgwA().first,sample->GetCgwA().first*sigma);
    this->fCgwB.first     = r->Gaus(sample->GetCgwB().first,sample->GetCgwB().first*sigma);
    //this->fCgwC.first     = r->Gaus(sample->GetCgwC().first,sample->GetCgwC().first*sigma);
  }
  if(apply_all_var || apply_only_shift_var)
    this->fCshift.first   = r->Gaus(sample->GetCshift().first,sample->GetCshift().first*sigma);
  if(apply_all_var || apply_only_norm_var)
    this->fCnorm.first = r->Gaus(sample->GetCnorm().first,sample->GetCnorm().first*sigma);

  std::cout << "MIGUEEEEE " << this->fCnorm.first << std::endl;
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
      fIFit.push_back(CoherentFitUtils::LangausFit(fh[i],GetSampleType()));
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
  std::cout << fType << std::endl;
  
  if(fType == SampleTypeEnum::kTrueSignal){
    //CoherentFitSignal();
    //CoherentFitSignalCheb();
    //CoherentFitSignalLag();
    CoherentFitSignalLagMPV();
  }
  else if(fType == SampleTypeEnum::kTrueBackground){
    //CoherentFitBackgroundQuadraticWidths();
    //CoherentFitBackgroundCheb();
    //CoherentFitBackgroundLag();
    CoherentFitBackgroundLagMPV();
  }
  else if(fType == SampleTypeEnum::kSignalPlusBackground){
    //CoherentFitSignalPlusBackgroundQuadraticWidths();
    //CoherentFitSignalPlusBackgroundCheb();
    //CoherentFitSignalPlusBackgroundLag();
    CoherentFitSignalPlusBackgroundLagMPV();
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
  for(int irr = 0; irr < (int)fRR.size(); irr++){
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
    mpv_error.push_back(sqrt(error));
    error = 0;
  }

  return new TGraphErrors(mpv.size(),&rr[0],&mpv[0],&rr_error[0],&mpv_error[0]);
}

//********************************************************************
void CoherentSample::CopyHistogramsBinning(std::vector<TH1F*> vh){
//********************************************************************

  if(fh.size() != vh.size()){
    std::cout << "histogram vector does not match size" << std::endl;
    std::exit(1);
  }

  for(UInt_t i = 0; i < fh.size(); i++)
    CoherentFitUtils::CopyHistogramBinning(fh[i],vh[i]);
}

//********************************************************************
TGraphErrors* CoherentSample::GetImpvGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fImpv);
}

//********************************************************************
TGraphErrors* CoherentSample::GetIlwGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fIlw);
}

//********************************************************************
TGraphErrors* CoherentSample::GetIgwGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fIgw);
}

//********************************************************************
TGraphErrors* CoherentSample::GetInormGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fInorm);
}

//********************************************************************
TGraphErrors* CoherentSample::GetIfwGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fIfw);
}

//********************************************************************
std::vector<double> CoherentSample::GetRRVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fRR.size(); i++)
    v.push_back(fRR[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetRRVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fRR.size(); i++)
    v.push_back(fRR[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIlwVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIlw.size(); i++)
    v.push_back(fIlw[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIlwVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIlw.size(); i++)
    v.push_back(fIlw[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetImpvVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fImpv.size(); i++)
    v.push_back(fImpv[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetImpvVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fImpv.size(); i++)
    v.push_back(fImpv[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetInormVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fInorm.size(); i++)
    v.push_back(fInorm[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetInormVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fInorm.size(); i++)
    v.push_back(fInorm[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIgwVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIgw.size(); i++)
    v.push_back(fIgw[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIgwVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIgw.size(); i++)
    v.push_back(fIgw[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIfwVectorValue() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIfw.size(); i++)
    v.push_back(fIfw[i].first);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIfwVectorError() const {
//********************************************************************

  std::vector<double> v;
  v.clear();
  for(int i = 0; i < (int)fIfw.size(); i++)
    v.push_back(fIfw[i].second);

  return v;
}

//********************************************************************
std::vector<double> CoherentSample::GetIalphaVector() const {
//********************************************************************

  std::vector<double> alpha;
  alpha.clear();
  for(int i = 0; i < (int)fIfw.size(); i++)
    alpha.push_back(sqrt(1-pow(fIlw[i].first/fIfw[i].first,2)));

  return alpha;
  
}

//DEFAULT MODEL

//********************************************************************
void CoherentSample::CoherentFitSignal(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 11+1;
  const int NPAR = CPAR;
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
  fMinuit->mnparm(11, "norm",  mean_norm, 0.01, 0, 1.2, ierflg);
  
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
  fMinuit->GetParameter(11, fCnorm.first, fCnorm.second);
  
  StoreCoherentFits(); 
}

//********************************************************************
void CoherentSample::fcnSignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

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
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundQuadraticWidths);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[NPAR] = {2.6,  4.60, -1.98, 0.063, 36.9, -2.01, 0.09};
  Double_t step[NPAR]   = {0.1,   0.1,   0.1,   0.1,  0.1,  0.01, 0.01};
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = mean_norm / fInorm.size();
  fMinuit->mnparm(0,  "shift", vstart[0], step[0],   0, 10, ierflg);
  fMinuit->mnparm(1,  "lwA"  , vstart[1], step[1],   0, 10, ierflg);
  fMinuit->mnparm(2,  "lwB"  , vstart[2], step[2], -10,  0, ierflg);
  fMinuit->mnparm(3,  "lwC"  , vstart[3], step[3],   0,  1, ierflg);
  fMinuit->mnparm(4,  "gwA"  , vstart[4], step[4],   0, 90, ierflg);
  fMinuit->mnparm(5,  "gwB"  , vstart[5], step[5], -10,  0, ierflg);
  fMinuit->mnparm(6,  "gwC"  , vstart[6], step[6],   0,  1, ierflg);
  fMinuit->mnparm(7,  "norm" , mean_norm, 0.01,      0,  1, ierflg);
  
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
  fMinuit->GetParameter(7, fCnorm.first,  fCnorm.second);

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
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundQuadraticWidths);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,smpvR,sgwA,sgwB,sgwC,blwA,blwB,blwC,bshift,bgwA,bgwB,bgwC,norm,slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,smpvR_error,sgwA_error,sgwB_error,sgwC_error,blwA_error,blwB_error,blwC_error,bshift_error,bgwA_error,bgwB_error,bgwC_error,norm_error;

  slwA   = GetSignal()->GetClwA().first;         
  slwB   = GetSignal()->GetClwB().first;	 
  slwC   = GetSignal()->GetClwC().first;	 
  smpvA  = GetSignal()->GetCmpvA().first;	 
  smpvB  = GetSignal()->GetCmpvB().first;	 
  smpvC  = GetSignal()->GetCmpvC().first;	 
  smpvD  = GetSignal()->GetCmpvD().first;	 
  smpvR  = GetSignal()->GetCmpvR().first;	 
  sgwA   = GetSignal()->GetCgwA().first;	 
  sgwB   = GetSignal()->GetCgwB().first;	 
  sgwC   = GetSignal()->GetCgwC().first;	 
  blwA   = GetBackground()->GetClwA().first;	 
  blwB   = GetBackground()->GetClwB().first;	 
  blwC   = GetBackground()->GetClwC().first;	 
  bshift = GetBackground()->GetCshift().first;   
  bgwA   = GetBackground()->GetCgwA().first;	 
  bgwB   = GetBackground()->GetCgwB().first;	 
  bgwC   = GetBackground()->GetCgwC().first;
  norm   = (GetSignal()->GetCnorm().first + (1-GetBackground()->GetCnorm().first))/2;
  
  Double_t vstart[NPAR] = {slwA, slwB, slwC, smpvA, smpvB, smpvC, smpvD, smpvR, sgwA, sgwB, sgwC,
			   blwA, blwB, blwC, bshift,                            bgwA, bgwB, bgwC,
                           norm};
  Double_t step[NPAR]   = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  
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
  fMinuit->mnparm(14, "b shift", vstart[14], step[14],   0,  0, ierflg);
  fMinuit->mnparm(15, "b gw A" , vstart[15], step[15],   0, 90, ierflg);
  fMinuit->mnparm(16, "b gw B" , vstart[16], step[16], -10,  0, ierflg);
  fMinuit->mnparm(17, "b gw C" , vstart[17], step[17],   0, 10, ierflg);
  fMinuit->mnparm(18, "norm"   , vstart[18], step[18], 0.95*vstart[18], 1.05*vstart[18], ierflg);
  
  fMinuit->FixParameter(13);
  fMinuit->FixParameter(17);
  
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
  fMinuit->GetParameter(18, norm  , norm_error);

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
  GetSignal()->SetCnorm(std::make_pair(norm,norm_error));

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
  GetBackground()->SetCnorm(std::make_pair(1-norm,norm_error));
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundQuadraticWidths(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

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
    snorm = norm;
    blw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,blw_par);
    bmpv  = CoherentFitUtils::ABCDRParametrization(RR_par,bmpv_par);
    bgw   = CoherentFitUtils::QuadraticABCParametrization(RR_par,bgw_par);
    bnorm = 1-norm;

    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }

  f = -Likelihood;
}

//Chebyshev polynomials models
//********************************************************************
void CoherentSample::CoherentFitSignalCheb(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 11+1;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalCheb);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
  //double mean_norm = 0;
  //for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  //mean_norm = mean_norm / fInorm.size();
 
  fMinuit->mnparm(0,  "lw A" , 3.29  , 100, 0, 0, ierflg);//0.9*3.29  , 1.1*3.29  , ierflg);//
  fMinuit->mnparm(1,  "lw B" , -4.57 , 100, 0, 0, ierflg);//0.9*-4.57 , 1.1*-4.57 , ierflg);//
  fMinuit->mnparm(2,  "lw C" , 1.42  , 100, 0, 0, ierflg);//0.9*1.42  , 1.1*1.42  , ierflg);//
  fMinuit->mnparm(3,  "mpv A", 0.65  , 100, 0, 0, ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(4,  "mpv B", 9.48  , 100, 0, 0, ierflg);//0.9*9.48  , 1.1*9.48  , ierflg);//
  fMinuit->mnparm(5,  "mpv C", -11.50, 100, 0, 0, ierflg);//0.9*-11.50, 1.1*-11.50, ierflg);//
  fMinuit->mnparm(6,  "mpv D", 5.20  , 100, 0, 0, ierflg);//0.9*5.20  , 1.1*5.20  , ierflg);//
  fMinuit->mnparm(7,  "mpv R", -2.10 , 100, 0, 0, ierflg);//0.9*-2.10 , 1.1*-2.10 , ierflg);//
  fMinuit->mnparm(8,  "gw A" , 2.34  , 100, 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(9,  "gw B" , -2.33 , 100, 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(10, "gw C" , 0.0   , 100, 0, 0, ierflg);//0.9*0.05  , 1.1*-0.05 , ierflg);//
  fMinuit->mnparm(11, "norm" , 0.88  , 100, 0, 0, ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//
  fMinuit->FixParameter(10);
  
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
  fMinuit->GetParameter(11, fCnorm.first, fCnorm.second);
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

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
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x = (CSMinuit->fRR[ibin].first-1)/(CSMinuit->fRR[ibin].first+1);
    lw   = ROOT::Math::Chebyshev2(x,lwA,lwB,lwC);
    mpv  = ROOT::Math::Chebyshev4(x,mpvA,mpvB,mpvC,mpvD,mpvR);
    gw   = ROOT::Math::Chebyshev2(x,gwA,gwB,gwC);
    flg->SetParameters(lw,mpv,norm,gw);
    /*if(ibin==0){
      CSMinuit->fh[ibin]->Draw();
      flg->Draw("same");
      gPad->Update();
      }*/
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitBackgroundCheb(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 2+1;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundCheb);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  fMinuit->mnparm(0, "shift_l" , -1.11, 10, 0, 0, ierflg);//-1.5
  fMinuit->mnparm(1, "shift_m" , 3.78 , 10,  0, 0, ierflg);
  //fMinuit->mnparm(2, "shift_g" , -1.5 , 10, -3, 0, ierflg);
  //fMinuit->mnparm(2, "gwA"     , -8.8 , 10,  0, 0, ierflg);
  //fMinuit->mnparm(3, "gwB"     , 14   , 10,  0, 0, ierflg);
  //fMinuit->mnparm(4, "gwV"     , -5.5 , 10,  0, 0, ierflg);
  //fMinuit->mnparm(0, "shift_l" , 2, 10, -5, 5, ierflg);//-1.5
  fMinuit->mnparm(2, "norm"    , 0.13 , 10,  0, 1, ierflg);
  //fMinuit->FixParameter(0);
  //fMinuit->FixParameter(2);
  //fMinuit->FixParameter(3);
  //fMinuit->FixParameter(4);
  
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
  fMinuit->GetParameter(0, fClwA.first,   fClwA.second);
  fMinuit->GetParameter(1, fCshift.first, fCshift.second);
  //fMinuit->GetParameter(2, fCgwA.first,   fCgwA.second);
  //fMinuit->GetParameter(3, fCgwB.first,   fCgwB.second);
  //fMinuit->GetParameter(4, fCgwC.first,   fCgwC.second);
  fMinuit->GetParameter(2, fCnorm.first,  fCnorm.second);

  fCmpvA = fTrueSignal->GetCmpvA();
  fCmpvB = fTrueSignal->GetCmpvB();
  fCmpvC = fTrueSignal->GetCmpvC();
  fCmpvD = fTrueSignal->GetCmpvD();
  fCmpvR = fTrueSignal->GetCmpvR();
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  //set each parameter
  double shift_l = par[0];
  double shift_m = par[1];
  //double shift_g = par[2];
  //double gwA     = par[2];
  //double gwB     = par[3];
  //double gwC     = par[4];
  double norm  = par[2];
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
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    //double x   = (CSMinuit->fRR[ibin].first-1        )/(CSMinuit->fRR[ibin].first+1        );
    double x_l = (CSMinuit->fRR[ibin].first-1+shift_l)/(CSMinuit->fRR[ibin].first+1+shift_l);
    double x_m = (CSMinuit->fRR[ibin].first-1+shift_m)/(CSMinuit->fRR[ibin].first+1+shift_m);
    //double x_g = (CSMinuit->fRR[ibin].first-1+shift_g)/(CSMinuit->fRR[ibin].first+1+shift_g);
    lw   = ROOT::Math::Chebyshev2(x_l,slwA,slwB,slwC);
    mpv  = ROOT::Math::Chebyshev4(x_m,smpvA,smpvB,smpvC,smpvD,smpvR);
    gw   = ROOT::Math::Chebyshev2(x_l,sgwA,sgwB,sgwC);
    //gw   = ROOT::Math::Chebyshev2(x,gwA,gwB,gwC);
    //norm = par[CPAR+ibin];
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundCheb(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 13+1;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundCheb);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //initial par
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,smpvR,sgwA,sgwB,sgwC,bshift_l,bshift_m,bshift_g,norm,slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,smpvR_error,sgwA_error,sgwB_error,sgwC_error,bshift_l_error,bshift_m_error,bshift_g_error,norm_error;

  slwA     = GetSignal()->GetClwA().first;         
  slwB     = GetSignal()->GetClwB().first;	 
  slwC     = GetSignal()->GetClwC().first;	 
  smpvA    = GetSignal()->GetCmpvA().first;	 
  smpvB    = GetSignal()->GetCmpvB().first;	 
  smpvC    = GetSignal()->GetCmpvC().first;	 
  smpvD    = GetSignal()->GetCmpvD().first;	 
  smpvR    = GetSignal()->GetCmpvR().first;	 
  sgwA     = GetSignal()->GetCgwA().first;	 
  sgwB     = GetSignal()->GetCgwB().first;	 
  sgwC     = GetSignal()->GetCgwC().first;	 
  bshift_l = GetBackground()->GetClwA().first;	 
  bshift_m = GetBackground()->GetCshift().first;   
  //bshift_g = GetBackground()->GetCgwA().first;	 
  norm     = (GetSignal()->GetCnorm().first + (1-GetBackground()->GetCnorm().first))/2;
  
  fMinuit->mnparm(0,  "lw A"    , slwA    , 100, 0, 0, ierflg);//0.9*3.29  , 1.1*3.29  , ierflg);//
  fMinuit->mnparm(1,  "lw B"    , slwB    , 100, 0, 0, ierflg);//0.9*-4.57 , 1.1*-4.57 , ierflg);//
  fMinuit->mnparm(2,  "lw C"    , slwC    , 100, 0, 0, ierflg);//0.9*1.42  , 1.1*1.42  , ierflg);//
  fMinuit->mnparm(3,  "mpv A"   , smpvA   , 100, 0, 0, ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(4,  "mpv B"   , smpvB   , 100, 0, 0, ierflg);//0.9*9.48  , 1.1*9.48  , ierflg);//
  fMinuit->mnparm(5,  "mpv C"   , smpvC   , 100, 0, 0, ierflg);//0.9*-11.50, 1.1*-11.50, ierflg);//
  fMinuit->mnparm(6,  "mpv D"   , smpvD   , 100, 0, 0, ierflg);//0.9*5.20  , 1.1*5.20  , ierflg);//
  fMinuit->mnparm(7,  "mpv R"   , smpvR   , 100, 0, 0, ierflg);//0.9*-2.10 , 1.1*-2.10 , ierflg);//
  fMinuit->mnparm(8,  "gw A"    , sgwA    , 100, 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(9,  "gw B"    , sgwB    , 100, 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(10, "gw C"    , sgwC    , 100, 0, 0, ierflg);//0.9*0.05  , 1.1*-0.05 , ierflg);//
  fMinuit->mnparm(11, "bshift_l", bshift_l, 100, 0, 0, ierflg);//-1.5
  fMinuit->mnparm(12, "bshift_m", bshift_m, 100, 0, 0, ierflg);
  //fMinuit->mnparm(13, "bshift_g", bshift_g, 100, -5, 0, ierflg);
  fMinuit->mnparm(13, "norm"    , norm    , 100, 0, 1, ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//
  fMinuit->FixParameter(10);
  
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
  fMinuit->GetParameter(0,  slwA     , slwA_error);
  fMinuit->GetParameter(1,  slwB     , slwB_error);
  fMinuit->GetParameter(2,  slwC     , slwC_error);
  fMinuit->GetParameter(3,  smpvA    , smpvA_error);
  fMinuit->GetParameter(4,  smpvB    , smpvB_error);
  fMinuit->GetParameter(5,  smpvC    , smpvC_error);
  fMinuit->GetParameter(6,  smpvD    , smpvD_error);
  fMinuit->GetParameter(7,  smpvR    , smpvR_error);
  fMinuit->GetParameter(8,  sgwA     , sgwA_error);
  fMinuit->GetParameter(9,  sgwB     , sgwB_error);
  fMinuit->GetParameter(10, sgwC     , sgwC_error);
  fMinuit->GetParameter(11, bshift_l , bshift_l_error);
  fMinuit->GetParameter(12, bshift_m , bshift_m_error);
  //fMinuit->GetParameter(13, bshift_g , bshift_g_error);
  fMinuit->GetParameter(13, norm     , norm_error);

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
  GetSignal()->SetCnorm(std::make_pair(norm,norm_error));

  GetBackground()->SetClwA(std::make_pair(bshift_l,bshift_l_error));
  GetBackground()->SetCmpvA(GetSignal()->GetCmpvA());
  GetBackground()->SetCmpvB(GetSignal()->GetCmpvB());
  GetBackground()->SetCmpvC(GetSignal()->GetCmpvC());
  GetBackground()->SetCmpvD(GetSignal()->GetCmpvD());
  GetBackground()->SetCmpvR(GetSignal()->GetCmpvR());
  GetBackground()->SetCshift(std::make_pair(bshift_m,bshift_m_error));
  GetBackground()->SetCgwA(std::make_pair(bshift_l,bshift_l_error));
  GetBackground()->SetCnorm(std::make_pair(1-norm,norm_error));
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundCheb(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  //set each parameter
  double lwA     = par[0];
  double lwB     = par[1];
  double lwC     = par[2];
  double mpvA    = par[3];
  double mpvB    = par[4];
  double mpvC    = par[5];
  double mpvD    = par[6];
  double mpvR    = par[7];
  double gwA     = par[8];
  double gwB     = par[9];
  double gwC     = par[10];
  double shift_l = par[11];
  double shift_m = par[12];
  //double shift_g = par[13];
  double norm    = par[13];

  //define function to fit to each histogram
  double slw,smpv,sgw,snorm,blw,bmpv,bgw,bnorm;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = (CSMinuit->fRR[ibin].first-1        )/(CSMinuit->fRR[ibin].first+1        );
    double x_l = (CSMinuit->fRR[ibin].first-1+shift_l)/(CSMinuit->fRR[ibin].first+1+shift_l);
    double x_m = (CSMinuit->fRR[ibin].first-1+shift_m)/(CSMinuit->fRR[ibin].first+1+shift_m);
    //double x_g = (CSMinuit->fRR[ibin].first-1+shift_g)/(CSMinuit->fRR[ibin].first+1+shift_g);
    slw   = ROOT::Math::Chebyshev2(x  ,lwA ,lwB ,lwC);
    smpv  = ROOT::Math::Chebyshev4(x  ,mpvA,mpvB,mpvC,mpvD,mpvR);
    sgw   = ROOT::Math::Chebyshev2(x  ,gwA ,gwB ,gwC);
    snorm = norm;
    blw   = ROOT::Math::Chebyshev2(x_l,lwA ,lwB ,lwC);
    bmpv  = ROOT::Math::Chebyshev4(x_m,mpvA,mpvB,mpvC,mpvD,mpvR);
    bgw   = ROOT::Math::Chebyshev2(x_l,gwA ,gwB ,gwC);
    bnorm = 1-norm;

    //    std::cout << slw << " " << smpv << " " << sgw << " " << blw << " " << bmpv << " " << bgw << std::endl;
    
    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//laguerre polynomials models
//********************************************************************
void CoherentSample::CoherentFitSignalLag(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 10;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalLag);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = 0.88;//mean_norm / fInorm.size();
 
  fMinuit->mnparm(0,  "lw A" , 2.9  , 1, 0, 0, ierflg);//0.9*3.29  , 1.1*3.29  , ierflg);//
  fMinuit->mnparm(1,  "lw B" , -6.6 , 1, 0, 0, ierflg);//0.9*-4.57 , 1.1*-4.57 , ierflg);//
  fMinuit->mnparm(2,  "lw C" , 5.5  , 1, 0, 0, ierflg);//0.9*1.42  , 1.1*1.42  , ierflg);//
  fMinuit->mnparm(3,  "mpv A", -67  , 1, 0, 0, ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(4,  "mpv B", 254  , 1, 0, 0, ierflg);//0.9*9.48  , 1.1*9.48  , ierflg);//
  fMinuit->mnparm(5,  "mpv C", -290 , 1, 0, 0, ierflg);//0.9*-11.50, 1.1*-11.50, ierflg);//
  fMinuit->mnparm(6,  "mpv D", 115  , 1, 0, 0, ierflg);//0.9*5.20  , 1.1*5.20  , ierflg);//
  fMinuit->mnparm(7,  "gw A" , 0.07 , 1, 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(8,  "gw B" , 2.17 , 1, 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(9,  "norm" , mean_norm , 1, 0, 0, ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//
    
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
  fMinuit->GetParameter(7,  fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(8,  fCgwB.first,  fCgwB.second);
  fMinuit->GetParameter(9, fCnorm.first, fCnorm.second);
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lwA  = par[0];
  double lwB  = par[1];
  double lwC  = par[2];
  double mpvA = par[3];
  double mpvB = par[4];
  double mpvC = par[5];
  double mpvD = par[6];
  double gwA  = par[7];
  double gwB  = par[8];
  double norm = par[9];

  //define function to fit to each histogram
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x = (CSMinuit->fRR[ibin].first-1)/(CSMinuit->fRR[ibin].first+1);
    lw   = lwA*ROOT::Math::laguerre(0,x) +lwB*ROOT::Math::laguerre(1,x) +lwC*ROOT::Math::laguerre(2,x);
    mpv  = mpvA*ROOT::Math::laguerre(0,x)+mpvB*ROOT::Math::laguerre(1,x)+mpvC*ROOT::Math::laguerre(2,x)+mpvD*ROOT::Math::laguerre(3,x);
    gw   = gwA*ROOT::Math::laguerre(0,x) +gwB*ROOT::Math::laguerre(1,x);
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitBackgroundLag(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 8;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundLag);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  fMinuit->mnparm(0,  "lw A" , 3.9 , 1, 0, 0, ierflg);//0.9*3.29  , 1.1*3.29  , ierflg);//      3.9 
  fMinuit->mnparm(1,  "lw B" , -8.2 , 1, 0, 0, ierflg);//0.9*-4.57 , 1.1*-4.57 , ierflg);//	     -8.2
  fMinuit->mnparm(2,  "lw C" ,  7.6, 1, 0, 0, ierflg);//0.9*1.42  , 1.1*1.42  , ierflg);//	     7.6 
  fMinuit->mnparm(3,  "mpv S", 4.3  , 1, 0, 0, ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(4,  "gw A" , -8.5 , 1, 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//      0.14
  fMinuit->mnparm(5,  "gw B" , 23.4 , 1, 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//	     2.35
  fMinuit->mnparm(6,  "gw C" , -16.6 , 1, 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//	     1.  
  fMinuit->mnparm(7,  "norm" , 0.1  , 1, 0, 0, ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//

  
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
  fMinuit->GetParameter(0, fClwA.first,   fClwA.second);
  fMinuit->GetParameter(1, fClwB.first,   fClwB.second);
  fMinuit->GetParameter(2, fClwC.first,   fClwC.second);
  fMinuit->GetParameter(3, fCshift.first, fCshift.second);
  fMinuit->GetParameter(4, fCgwA.first,   fCgwA.second);
  fMinuit->GetParameter(5, fCgwB.first,   fCgwB.second);
  fMinuit->GetParameter(6, fCgwC.first,   fCgwC.second);
  fMinuit->GetParameter(7, fCnorm.first,  fCnorm.second);

  fCmpvA = fTrueSignal->GetCmpvA();
  fCmpvB = fTrueSignal->GetCmpvB();
  fCmpvC = fTrueSignal->GetCmpvC();
  fCmpvD = fTrueSignal->GetCmpvD();
  fCmpvR = fTrueSignal->GetCmpvR();

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  //set each parameter
  double lwA   = par[0];
  double lwB   = par[1];
  double lwC   = par[2];
  double shift = par[3];
  double gwA   = par[4];
  double gwB   = par[5];
  double gwC   = par[6];
  double norm  = par[7];

  // double slwA = CSMinuit->GetTrueSignal()->GetClwA().first;
  // double slwB = CSMinuit->GetTrueSignal()->GetClwB().first;
  // double slwC = CSMinuit->GetTrueSignal()->GetClwC().first;
  
  double smpvA = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  double smpvB = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  double smpvC = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  double smpvD = CSMinuit->GetTrueSignal()->GetCmpvD().first;
  
  //define function to fit to each histogram
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = (CSMinuit->fRR[ibin].first-1)/(CSMinuit->fRR[ibin].first+1);
    //double x_lw = (CSMinuit->fRR[ibin].first-1+lwA)/(CSMinuit->fRR[ibin].first+1+lwA);
    double x_mpv = (CSMinuit->fRR[ibin].first-1+shift)/(CSMinuit->fRR[ibin].first+1+shift);
    //lw   = sqrt(pow(slwA*ROOT::Math::laguerre(0,x_lw) +slwB*ROOT::Math::laguerre(1,x_lw) +slwC*ROOT::Math::laguerre(2,x_lw),2)+pow(lwB,2));
    lw   = lwA*ROOT::Math::laguerre(0,x) +lwB*ROOT::Math::laguerre(1,x) +lwC*ROOT::Math::laguerre(2,x);
    mpv  = smpvA*ROOT::Math::laguerre(0,x_mpv)+smpvB*ROOT::Math::laguerre(1,x_mpv)+smpvC*ROOT::Math::laguerre(2,x_mpv)+smpvD*ROOT::Math::laguerre(3,x_mpv);
    gw   = gwA*ROOT::Math::laguerre(0,x) +gwB*ROOT::Math::laguerre(1,x) +gwC*ROOT::Math::laguerre(2,x);
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundLag(){
//********************************************************************

  CSMinuit = this;
  
  //create Minuit
  const int CPAR = 17;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundLag);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //initial par
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,sgwA,sgwB,norm;
  double slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,sgwA_error,sgwB_error, norm_error;
  double blwA,blwB,blwC,bmpvS,bgwA,bgwB,bgwC;
  double blwA_error,blwB_error,blwC_error,bmpvS_error,bgwA_error,bgwB_error,bgwC_error;

  slwA = CSMinuit->GetTrueSignal()->GetClwA().first;
  slwB = CSMinuit->GetTrueSignal()->GetClwB().first;
  slwC = CSMinuit->GetTrueSignal()->GetClwC().first;

  smpvA = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  smpvB = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  smpvC = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  smpvD = CSMinuit->GetTrueSignal()->GetCmpvD().first;

  sgwA = CSMinuit->GetTrueSignal()->GetCgwA().first;
  sgwB = CSMinuit->GetTrueSignal()->GetCgwB().first;

  blwA = CSMinuit->GetTrueBackground()->GetClwA().first;
  blwB = CSMinuit->GetTrueBackground()->GetClwB().first;
  blwC = CSMinuit->GetTrueBackground()->GetClwC().first;

  bmpvS = CSMinuit->GetTrueBackground()->GetCshift().first;

  bgwA = CSMinuit->GetTrueBackground()->GetCgwA().first;
  bgwB = CSMinuit->GetTrueBackground()->GetCgwB().first;
  bgwC = CSMinuit->GetTrueBackground()->GetCgwC().first;

  norm = 1-CSMinuit->GetTrueBackground()->GetCnorm().first;

  slwA_error = CSMinuit->GetTrueSignal()->GetClwA().second;
  slwB_error = CSMinuit->GetTrueSignal()->GetClwB().second;
  slwC_error = CSMinuit->GetTrueSignal()->GetClwC().second;
      
  smpvA_error = CSMinuit->GetTrueSignal()->GetCmpvA().second;
  smpvB_error = CSMinuit->GetTrueSignal()->GetCmpvB().second;
  smpvC_error = CSMinuit->GetTrueSignal()->GetCmpvC().second;
  smpvD_error = CSMinuit->GetTrueSignal()->GetCmpvD().second;
      
  sgwA_error = CSMinuit->GetTrueSignal()->GetCgwA().second;
  sgwB_error = CSMinuit->GetTrueSignal()->GetCgwB().second;
      
  blwA_error = CSMinuit->GetTrueBackground()->GetClwA().second;
  blwB_error = CSMinuit->GetTrueBackground()->GetClwB().second;
  blwC_error = CSMinuit->GetTrueBackground()->GetClwC().second;
      
  bmpvS_error = CSMinuit->GetTrueBackground()->GetCshift().second;
      
  bgwA_error = CSMinuit->GetTrueBackground()->GetCgwA().second;
  bgwB_error = CSMinuit->GetTrueBackground()->GetCgwB().second;
  bgwC_error = CSMinuit->GetTrueBackground()->GetCgwC().second;

  norm_error = 1-CSMinuit->GetTrueBackground()->GetCnorm().second;
  
  fMinuit->mnparm(0,  "lw A"      , slwA  , 10*slwA_error  , 0, 0, ierflg);  // 2.73938e+00 
  fMinuit->mnparm(1,  "lw B"      , slwB  , 10*slwB_error  , 0, 0, ierflg);  //-6.30265e+00
  fMinuit->mnparm(2,  "lw C"      , slwC  , 10*slwC_error  , 0, 0, ierflg);  // 5.22968e+00
  fMinuit->mnparm(3,  "mpv A"     , smpvA , 10*smpvA_error , 0, 0, ierflg);  //-7.80112e+01
  fMinuit->mnparm(4,  "mpv B"     , smpvB , 10*smpvB_error , 0, 0, ierflg);  // 2.95159e+02
  fMinuit->mnparm(5,  "mpv C"     , smpvC , 10*smpvC_error , 0, 0, ierflg);  //-3.41434e+02
  fMinuit->mnparm(6,  "mpv D"     , smpvD , 10*smpvD_error , 0, 0, ierflg);  // 1.36370e+02
  fMinuit->mnparm(7,  "gw A"      , sgwA  , 10*sgwA_error  , 0, 0, ierflg);  // 7.61715e-02
  fMinuit->mnparm(8,  "gw B"      , sgwB  , 10*sgwB_error  , 0, 0, ierflg);  // 2.02543e+00
  fMinuit->mnparm(9,  "blw A"     , blwA  , 10*blwA_error  , 0, 0, ierflg);  // 5.96330e-01
  fMinuit->mnparm(10, "blw B"     , blwB  , 10*blwB_error  , 0, 0, ierflg);  // 4.02872e-02
  fMinuit->mnparm(11, "blw C"     , blwC  , 10*blwC_error  , 0, 0, ierflg);  // 4.02872e-02
  fMinuit->mnparm(12, "mpv shift" , bmpvS , 10*bmpvS_error , 0.9*bmpvS, 1.1*bmpvS, ierflg);  // 3.14637e+00
  fMinuit->mnparm(13,  "bgw A"    , bgwA  , 10*bgwA_error  , 0, 0, ierflg);  // 4.21610e+00
  fMinuit->mnparm(14,  "bgw B"    , bgwB  , 10*bgwB_error  , 0, 0, ierflg);  //-5.62652e+00
  fMinuit->mnparm(15,  "bgw C"    , bgwC  , 10*bgwC_error  , 0, 0, ierflg);  // 8.24227e+00
  fMinuit->mnparm(16, "norm"      , norm  , 10*norm_error  , 0.9*norm, 1.1*norm, ierflg);// 8.74262e-01

  // fMinuit->FixParameter(0);
  // fMinuit->FixParameter(1);
  // fMinuit->FixParameter(2);
  // fMinuit->FixParameter(3);
  // fMinuit->FixParameter(4);
  // fMinuit->FixParameter(5);
  // fMinuit->FixParameter(6);
  // fMinuit->FixParameter(7);
  // fMinuit->FixParameter(8);
  // // fMinuit->FixParameter(9);
  // // fMinuit->FixParameter(10);
  // // fMinuit->FixParameter(11);
  // // fMinuit->FixParameter(12);
  // // fMinuit->FixParameter(13);
  // // fMinuit->FixParameter(14);    
  // // fMinuit->FixParameter(15);
  // fMinuit->FixParameter(16); 

  // Set strategy to 2, we want TMinuit to try its best
  arglist[0] = 2;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  
  // // Now ready for minimization step
  // arglist[0] = 50000;
  // arglist[1] = 1.;
  // fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  //first fit exclusively norm and shift
  // fMinuit->FixParameter(0); 
  // fMinuit->FixParameter(1); 
  // fMinuit->FixParameter(2); 
  // fMinuit->FixParameter(3); 
  // fMinuit->FixParameter(4); 
  // fMinuit->FixParameter(5); 
  // fMinuit->FixParameter(6); 
  // fMinuit->FixParameter(7); 
  // fMinuit->FixParameter(8); 
  // fMinuit->FixParameter(9); 
  // fMinuit->FixParameter(10);
  // fMinuit->FixParameter(11);
  // fMinuit->FixParameter(12);
  // fMinuit->FixParameter(13);
  // fMinuit->FixParameter(14);
  // fMinuit->FixParameter(15);
  // fMinuit->FixParameter(16);
      
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
  fMinuit->GetParameter(7,  sgwA  , sgwA_error);
  fMinuit->GetParameter(8,  sgwB  , sgwB_error);
  fMinuit->GetParameter(9,  blwA  , blwA_error);
  fMinuit->GetParameter(10, blwB  , blwB_error);
  fMinuit->GetParameter(11, blwC  , blwC_error);
  fMinuit->GetParameter(12, bmpvS , bmpvS_error);
  fMinuit->GetParameter(13, bgwA  , bgwA_error);
  fMinuit->GetParameter(14, bgwB  , bgwB_error);
  fMinuit->GetParameter(15, bgwC  , bgwC_error);
  fMinuit->GetParameter(16, norm  , norm_error);

  GetSignal()->SetClwA(std::make_pair(slwA,slwA_error));
  GetSignal()->SetClwB(std::make_pair(slwB,slwB_error));
  GetSignal()->SetClwC(std::make_pair(slwC,slwC_error));
  GetSignal()->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  GetSignal()->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  GetSignal()->SetCmpvC(std::make_pair(smpvC,smpvC_error));
  GetSignal()->SetCmpvD(std::make_pair(smpvD,smpvD_error));
  GetSignal()->SetCgwA(std::make_pair(sgwA,sgwA_error));
  GetSignal()->SetCgwB(std::make_pair(sgwB,sgwB_error));
  GetSignal()->SetCnorm(std::make_pair(norm,norm_error));

  GetBackground()->SetClwA(std::make_pair(blwA,blwA_error));
  GetBackground()->SetClwB(std::make_pair(blwB,blwB_error));
  GetBackground()->SetClwC(std::make_pair(blwC,blwC_error));
  GetBackground()->SetCshift(std::make_pair(bmpvS,bmpvS_error));
  GetBackground()->SetCgwA(std::make_pair(bgwA,bgwA_error));
  GetBackground()->SetCgwB(std::make_pair(bgwB,bgwB_error));
  GetBackground()->SetCgwC(std::make_pair(bgwC,bgwC_error));
  GetBackground()->SetCnorm(std::make_pair(1-norm,norm_error));
  
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundLag(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************
  
  //set each parameter
  double slwA  = par[0];
  double slwB  = par[1];
  double slwC  = par[2];
  double smpvA = par[3];
  double smpvB = par[4];
  double smpvC = par[5];
  double smpvD = par[6];
  double sgwA  = par[7];
  double sgwB  = par[8];
  double blwA  = par[9];
  double blwB  = par[10];
  double blwC  = par[11];
  double bmpvS = par[12];
  double bgwA  = par[13];
  double bgwB  = par[14];
  double bgwC  = par[15];
  double norm  = par[16];

  //define function to fit to each histogram
  double slw,smpv,sgw,snorm,blw,bmpv,bgw,bnorm;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = (CSMinuit->fRR[ibin].first-1        )/(CSMinuit->fRR[ibin].first+1     );
    //double x_l = (CSMinuit->fRR[ibin].first-1+blwS) /(CSMinuit->fRR[ibin].first+1+blwS );
    double x_m = (CSMinuit->fRR[ibin].first-1+bmpvS)/(CSMinuit->fRR[ibin].first+1+bmpvS);

    slw   =  slwA*ROOT::Math::laguerre(0,x) +  slwB*ROOT::Math::laguerre(1,x) +  slwC*ROOT::Math::laguerre(2,x);
    smpv  = smpvA*ROOT::Math::laguerre(0,x) + smpvB*ROOT::Math::laguerre(1,x) + smpvC*ROOT::Math::laguerre(2,x) + smpvD*ROOT::Math::laguerre(3,x);
    sgw   =  sgwA*ROOT::Math::laguerre(0,x) +  sgwB*ROOT::Math::laguerre(1,x);
    snorm = norm;

    blw   = blwA*ROOT::Math::laguerre(0,x) + blwB*ROOT::Math::laguerre(1,x) + blwC*ROOT::Math::laguerre(2,x);
    bmpv  = smpvA*ROOT::Math::laguerre(0,x_m) + smpvB*ROOT::Math::laguerre(1,x_m) + smpvC*ROOT::Math::laguerre(2,x_m) + smpvD*ROOT::Math::laguerre(3,x_m);
    bgw   = bgwA*ROOT::Math::laguerre(0,x) + bgwB*ROOT::Math::laguerre(1,x) + bgwC*ROOT::Math::laguerre(2,x);
    bnorm = 1-norm;

    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//laguerre polynomials models mpv dependence
//********************************************************************
void CoherentSample::CoherentFitSignalLagMPV(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 11;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalLagMPV);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = 0.88;//mean_norm / fInorm.size();

  fMinuit->mnparm(0,  "lw A" , 1.95      , 0.1  , 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(1,  "lw B" , -1.94     , 0.1  , 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(2,  "lw C" , 0.49      , 0.1  , 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(3,  "mpv A", -78       , 10   , 0, 0, ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(4,  "mpv B",  295      , 10   , 0, 0, ierflg);//0.9*9.48  , 1.1*9.48  , ierflg);//
  fMinuit->mnparm(5,  "mpv C", -341      , 10   , 0, 0, ierflg);//0.9*-11.50, 1.1*-11.50, ierflg);//
  fMinuit->mnparm(6,  "mpv D",  136      , 10   , 0, 0, ierflg);//0.9*5.20  , 1.1*5.20  , ierflg);//
  fMinuit->mnparm(7,  "gw A" , 3.56      , 0.1  , 0, 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(8,  "gw B" , -3.85     , 0.1  , 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(9,  "gw C" , 1.11      , 0.1  , 0, 0, ierflg);//0.9*-2.18 , 1.1*-2.18 , ierflg);//
  fMinuit->mnparm(10, "norm" , mean_norm , 0.1  , 0, 1, ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//

  // Now ready for minimization step
  arglist[0] = 500000;
  arglist[1] = 1.;
  //fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
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
  fMinuit->GetParameter(7,  fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(8,  fCgwB.first,  fCgwB.second);
  fMinuit->GetParameter(9,  fCgwC.first,  fCgwC.second);
  fMinuit->GetParameter(10, fCnorm.first, fCnorm.second);
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lwA  = par[0];
  double lwB  = par[1];
  double lwC  = par[2];
  double mpvA = par[3];
  double mpvB = par[4];
  double mpvC = par[5];
  double mpvD = par[6];
  double gwA  = par[7];
  double gwB  = par[8];
  double gwC  = par[9];
  double norm = par[10];

  //define function to fit to each histogram
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = CSMinuit->fRR[ibin].first;
    double x_t = (x-1)/(x+1);
    mpv  = mpvA*ROOT::Math::laguerre(0,x_t)+mpvB*ROOT::Math::laguerre(1,x_t)+mpvC*ROOT::Math::laguerre(2,x_t)+mpvD*ROOT::Math::laguerre(3,x_t);

    lw   = lwA + lwB*sqrt(mpv)+lwC*mpv;
    gw   = gwA + gwB*sqrt(mpv)+gwC*mpv;
    //std::cout << lw << " " << mpv << " " << norm << " " << gw << std::endl;
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
    // CSMinuit->fh[ibin]->Draw();
    // flg->Draw("same");
    // gPad->Update();gPad->WaitPrimitive();
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitBackgroundLagMPV(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 4;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackgroundLagMPV);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
  double mean_norm = 0;
  for(int i = 0; i < (int)fInorm.size(); i++)mean_norm += fInorm[i].first;
  mean_norm = 0.1;//mean_norm / fInorm.size();

  fMinuit->mnparm(0,  "lw Q" , 0.1       , 10 , 0 , 0, ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(1,  "mpv S", 3.4       , 10 , 0 , 0 , ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(2,  "gw Q" , 0.1       , 10 , 0 , 0 , ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(3,  "norm" , mean_norm , 1  , 0 , 1 , ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//

  //fMinuit->FixParameter(0);
  //fMinuit->FixParameter(3);
  
  std::vector<TH1F*> h_rebin, h_save;
  for(int i = 0; i < (int)fh.size(); i++){
    h_rebin.push_back((TH1F*)fh[i]->Clone());
    h_rebin.back()->Rebin(2);
  }
  h_save = fh;
  fh = h_rebin;
    
  // Now ready for minimization step
  arglist[0] = 500000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
 
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fClwA.first,  fClwA.second);
  fMinuit->GetParameter(1,  fCshift.first,  fCshift.second);
  fMinuit->GetParameter(2,  fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(3,  fCnorm.first, fCnorm.second);
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackgroundLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lwQ   = par[0];
  double shift = par[1];
  double gwQ   = par[2];
  double norm  = par[3];

  double smpvA = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  double smpvB = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  double smpvC = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  double smpvD = CSMinuit->GetTrueSignal()->GetCmpvD().first;

  double slwA = CSMinuit->GetTrueSignal()->GetClwA().first;
  double slwB = CSMinuit->GetTrueSignal()->GetClwB().first;
  double slwC = CSMinuit->GetTrueSignal()->GetClwC().first;

  double sgwA = CSMinuit->GetTrueSignal()->GetCgwA().first;
  double sgwB = CSMinuit->GetTrueSignal()->GetCgwB().first;
  double sgwC = CSMinuit->GetTrueSignal()->GetCgwC().first;
      
  //define function to fit to each histogram
  double lw,mpv,gw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = CSMinuit->fRR[ibin].first+shift;
    double x_t = (x-1)/(x+1);
    mpv  = smpvA*ROOT::Math::laguerre(0,x_t)+smpvB*ROOT::Math::laguerre(1,x_t)+smpvC*ROOT::Math::laguerre(2,x_t)+smpvD*ROOT::Math::laguerre(3,x_t);
    lw   = sqrt(pow(slwA + slwB*sqrt(mpv)+slwC*mpv,2) + pow(lwQ*mpv,2));
    gw   = sqrt(pow(sgwA + sgwB*sqrt(mpv)+sgwC*mpv,2) + pow(gwQ*mpv,2));
    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundLagMPV(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 14;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackgroundLagMPV);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //initial par
  double slwA,slwB,slwC,smpvA,smpvB,smpvC,smpvD,sgwA,sgwB,sgwC,norm;
  double slwA_error,slwB_error,slwC_error,smpvA_error,smpvB_error,smpvC_error,smpvD_error,sgwA_error,sgwB_error,sgwC_error, norm_error;
  double blwQ,bmpvS,bgwQ;
  double blwQ_error,bmpvS_error,bgwQ_error;

  slwA = CSMinuit->GetTrueSignal()->GetClwA().first;
  slwB = CSMinuit->GetTrueSignal()->GetClwB().first;
  slwC = CSMinuit->GetTrueSignal()->GetClwC().first;

  smpvA = CSMinuit->GetTrueSignal()->GetCmpvA().first;
  smpvB = CSMinuit->GetTrueSignal()->GetCmpvB().first;
  smpvC = CSMinuit->GetTrueSignal()->GetCmpvC().first;
  smpvD = CSMinuit->GetTrueSignal()->GetCmpvD().first;

  sgwA = CSMinuit->GetTrueSignal()->GetCgwA().first;
  sgwB = CSMinuit->GetTrueSignal()->GetCgwB().first;
  sgwC = CSMinuit->GetTrueSignal()->GetCgwC().first;

  norm = CSMinuit->GetTrueSignal()->GetCnorm().first;
  
  blwQ = CSMinuit->GetTrueBackground()->GetClwA().first;

  bmpvS = CSMinuit->GetTrueBackground()->GetCshift().first;

  bgwQ = CSMinuit->GetTrueBackground()->GetCgwA().first;

  slwA_error = CSMinuit->GetTrueSignal()->GetClwA().second;
  slwB_error = CSMinuit->GetTrueSignal()->GetClwB().second;
  slwC_error = CSMinuit->GetTrueSignal()->GetClwC().second;
   
  smpvA_error = CSMinuit->GetTrueSignal()->GetCmpvA().second;
  smpvB_error = CSMinuit->GetTrueSignal()->GetCmpvB().second;
  smpvC_error = CSMinuit->GetTrueSignal()->GetCmpvC().second;
  smpvD_error = CSMinuit->GetTrueSignal()->GetCmpvD().second;
   
  sgwA_error = CSMinuit->GetTrueSignal()->GetCgwA().second;
  sgwB_error = CSMinuit->GetTrueSignal()->GetCgwB().second;
   
  blwQ_error = CSMinuit->GetTrueBackground()->GetClwA().second;
   
  bmpvS_error = CSMinuit->GetTrueBackground()->GetCshift().second;
   
  bgwQ_error = CSMinuit->GetTrueBackground()->GetCgwA().second;

  norm_error = CSMinuit->GetTrueBackground()->GetCnorm().second;

  fMinuit->mnparm(0,  "lw A"      , slwA  , 0.1 ,0.95*slwA   ,1.05*slwA   , ierflg);  // 2.73938e+00 
  fMinuit->mnparm(1,  "lw B"      , slwB  , 0.1 ,0.95*slwB   ,1.05*slwB   , ierflg);  //-6.30265e+00
  fMinuit->mnparm(2,  "lw C"      , slwC  , 0.1 ,0.95*slwC   ,1.05*slwC   , ierflg);  // 5.22968e+00
  fMinuit->mnparm(3,  "mpv A"     , smpvA , 0.1 ,0.95*smpvA  ,1.05*smpvA  , ierflg);  //-7.80112e+01
  fMinuit->mnparm(4,  "mpv B"     , smpvB , 0.1 ,0.95*smpvB  ,1.05*smpvB  , ierflg);  // 2.95159e+02
  fMinuit->mnparm(5,  "mpv C"     , smpvC , 0.1 ,0.95*smpvC  ,1.05*smpvC  , ierflg);  //-3.41434e+02
  fMinuit->mnparm(6,  "mpv D"     , smpvD , 0.1 ,0.95*smpvD  ,1.05*smpvD  , ierflg);  // 1.36370e+02
  fMinuit->mnparm(7,  "gw A"      , sgwA  , 0.1 ,0.95*sgwA   ,1.05*sgwA   , ierflg);  // 7.61715e-02
  fMinuit->mnparm(8,  "gw B"      , sgwB  , 0.1 ,0.95*sgwB   ,1.05*sgwB   , ierflg);  // 2.02543e+00
  fMinuit->mnparm(9,  "gw C"      , sgwC  , 0.1 ,0.95*sgwC   ,1.05*sgwC   , ierflg);  // 2.02543e+00a
  fMinuit->mnparm(10, "blw Q"     , blwQ  , 0.1 ,0,0  , ierflg);  // 5.96330e-01
  fMinuit->mnparm(11, "mpv shift" , bmpvS , 0.1 ,0,0  , ierflg);
  fMinuit->mnparm(12, "bgw Q"     , bgwQ  , 0.1 ,0,0  , ierflg);  // 4.21610e+00
  fMinuit->mnparm(13, "norm"      , norm  , 0.1 ,0.95*norm   ,1   , ierflg);// 8.74262e-01

  // fMinuit->FixParameter(0);
  // fMinuit->FixParameter(1);
  // fMinuit->FixParameter(2);
  // fMinuit->FixParameter(3);
  // fMinuit->FixParameter(4);
  // fMinuit->FixParameter(5);
  // fMinuit->FixParameter(6);
  // fMinuit->FixParameter(7);
  // fMinuit->FixParameter(8);
  // fMinuit->FixParameter(9);
  fMinuit->FixParameter(10);
  fMinuit->FixParameter(11);
  fMinuit->FixParameter(12);
  // fMinuit->FixParameter(13);

  // // Calcule hesse
  // arglist[0] = 50000;
  // fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  
  // // Now ready for minimization step
  // // Set strategy to 2, we want TMinuit to try its best
  // arglist[0] = 2;
  // fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 50000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // arglist[1] = 3;
  // arglist[2] = 4;
  // arglist[3] = 5;
  // arglist[4] = 6;
  // fMinuit->mnexcm("MINOS", arglist, 5, ierflg);

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
  fMinuit->GetParameter(7,  sgwA  , sgwA_error);
  fMinuit->GetParameter(8,  sgwB  , sgwB_error);
  fMinuit->GetParameter(9,  sgwC  , sgwC_error);
  fMinuit->GetParameter(10,  blwQ  , blwQ_error);
  fMinuit->GetParameter(11, bmpvS , bmpvS_error);
  fMinuit->GetParameter(12, bgwQ  , bgwQ_error);
  fMinuit->GetParameter(13, norm  , norm_error);

  GetSignal()->SetClwA(std::make_pair(slwA,slwA_error));
  GetSignal()->SetClwB(std::make_pair(slwB,slwB_error));
  GetSignal()->SetClwC(std::make_pair(slwC,slwC_error));
  GetSignal()->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  GetSignal()->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  GetSignal()->SetCmpvC(std::make_pair(smpvC,smpvC_error));
  GetSignal()->SetCmpvD(std::make_pair(smpvD,smpvD_error));
  GetSignal()->SetCgwA(std::make_pair(sgwA,sgwA_error));
  GetSignal()->SetCgwB(std::make_pair(sgwB,sgwB_error));
  GetSignal()->SetCgwC(std::make_pair(sgwC,sgwC_error));
  GetSignal()->SetCnorm(std::make_pair(norm,norm_error));

  GetBackground()->SetClwA(std::make_pair(blwQ,blwQ_error));
  // GetBackground()->SetClwB(std::make_pair(blwB,blwB_error));
  // GetBackground()->SetClwC(std::make_pair(blwC,blwC_error));
  GetBackground()->SetCshift(std::make_pair(bmpvS,bmpvS_error));
  GetBackground()->SetCgwA(std::make_pair(bgwQ,bgwQ_error));
  // GetBackground()->SetCgwB(std::make_pair(bgwB,bgwB_error));
  // GetBackground()->SetCgwC(std::make_pair(bgwC,bgwC_error));
  GetBackground()->SetCnorm(std::make_pair(1-norm,norm_error));

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackgroundLagMPV(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double slwA  = par[0];
  double slwB  = par[1];
  double slwC  = par[2];
  double smpvA = par[3];
  double smpvB = par[4];
  double smpvC = par[5];
  double smpvD = par[6];
  double sgwA  = par[7];
  double sgwB  = par[8];
  double sgwC  = par[9];
  double blwQ  = par[10];
  double bmpvS = par[11];
  double bgwQ  = par[12];
  double norm  = par[13];

  //define function to fit to each histogram
  double slw,smpv,sgw,snorm,blw,bmpv,bgw,bnorm;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = (CSMinuit->fRR[ibin].first-1      )/(CSMinuit->fRR[ibin].first+1     );
    double x_m = (CSMinuit->fRR[ibin].first-1+bmpvS)/(CSMinuit->fRR[ibin].first+1+bmpvS);

    smpv  = smpvA*ROOT::Math::laguerre(0,x) + smpvB*ROOT::Math::laguerre(1,x) + smpvC*ROOT::Math::laguerre(2,x) + smpvD*ROOT::Math::laguerre(3,x);
    slw   = slwA + slwB*sqrt(smpv) + slwC*smpv;
    sgw   = sgwA + sgwB*sqrt(smpv) + sgwC*smpv;
    snorm = norm;

    bmpv  = smpvA*ROOT::Math::laguerre(0,x_m) + smpvB*ROOT::Math::laguerre(1,x_m) + smpvC*ROOT::Math::laguerre(2,x_m) + smpvD*ROOT::Math::laguerre(3,x_m);
    blw   = sqrt(pow(slwA + slwB*sqrt(bmpv) + slwC*bmpv,2)+pow(blwQ*bmpv,2));
    bgw   = sqrt(pow(sgwA + sgwB*sqrt(bmpv) + sgwC*bmpv,2)+pow(bgwQ*bmpv,2));
    bnorm = 1-norm;

    flg->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}
