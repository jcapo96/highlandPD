#include "CoherentSample.hxx"
#include "CoherentFitUtils.hxx"
#include <vector>
#include "TVector.h"
#include "TParameter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "Math/Math.h"
#include "Math/ChebyshevPol.h"
#include "Math/SpecFuncMathMore.h"

#include <iostream>

static CoherentSample* CSMinuit;

bool debug = false;

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

  fLikelihood = -999.;
  fPartialLikelihood.clear();

  for(int i = 0; i < 18; i++)
    for(int j = 0; j < 18; j++)
      fcorr_matrix[i][j] = -999.;
}

//********************************************************************
CoherentSample::~CoherentSample(){
//********************************************************************

  //if(!(fMinuit==NULL))delete fMinuit;

  for(int i = 0; i < (int)fh.size(); i++)
    delete fh[i];
  fh.clear();

  for(int i = 0; i < (int)fIFit.size(); i++)
    delete fIFit[i];
  fIFit.clear();

  for(int i = 0; i < (int)fCFit.size(); i++)
    delete fCFit[i];
  fCFit.clear();
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
  fSemiSignal = c.GetSemiSignal();
  fBackground = c.GetBackground();
  fSemiBackground = c.GetSemiBackground();
  fTrueSignal = c.GetTrueSignal();
  fTrueSemiSignal = c.GetTrueSemiSignal();
  fTrueBackground = c.GetTrueBackground();
  fTrueSemiBackground = c.GetTrueSemiBackground();

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
<<<<<<< HEAD

=======
  fCconst = c.GetCconst();
>>>>>>> develop
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

<<<<<<< HEAD
  if(fClwFit)fClwFit->Write();
  if(fCmpvFit)fCmpvFit->Write();
  if(fCgwFit)fCgwFit->Write();

=======
  if(fType != SampleTypeEnum::kSignalPlusBackground){
    if(fClwFit)fClwFit->Write();
    if(fCmpvFit)fCmpvFit->Write();
    if(fCgwFit)fCgwFit->Write();
  }
  
>>>>>>> develop
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
<<<<<<< HEAD

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
=======
    if(fType == SampleTypeEnum::kTrueSignal || fType == SampleTypeEnum::kSignal)
      StoreCoherentFitsSignal();
    if(fType == SampleTypeEnum::kTrueBackground || fType == SampleTypeEnum::kBackground)
      StoreCoherentFitsBackground();
>>>>>>> develop
  }
  else{
    //store parametric functions
    GetSignal()->StoreCoherentFits();
    GetBackground()->StoreCoherentFits();
    
    //store individual histogram fits
    for(int i = 0; i < (int)fh.size(); i++){
      std::stringstream ssi;
      ssi << i;
      TF1* f = new TF1(("CF_sb"+ssi.str()+"").c_str(),CoherentFitUtils::TripleLangausPlusConstant,1,30,13);

<<<<<<< HEAD
      double slw, smpv, snorm, sgw, blw, bmpv, bnorm, bgw;

=======
      double slw, smpv, snorm, sgw, blw, bmpv, bnorm, bgw, sslw, ssmpv, ssgw, ssnorm, bC;
    
>>>>>>> develop
      slw   = GetSignal()->GetClwFit()->Eval(fRR[i].first);
      smpv  = GetSignal()->GetCmpvFit()->Eval(fRR[i].first);
      snorm = GetSignal()->GetCnorm().first;
      sgw   = GetSignal()->GetCgwFit()->Eval(fRR[i].first);
      blw   = GetBackground()->GetClwA().first;
      bmpv  = GetBackground()->GetCmpvFit()->Eval(fRR[i].first);
      bnorm = GetBackground()->GetCnorm().first;
<<<<<<< HEAD
      bgw   = GetBackground()->GetCgwFit()->Eval(fRR[i].first);

      f->SetParameters(slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw);
=======
      bgw   = GetBackground()->GetCgwA().first;
      bC    = GetBackground()->GetCconst().first;
      sslw   = GetBackground()->GetSemiBackground()->GetClwFit()->Eval(fRR[i].first);
      ssmpv  = GetBackground()->GetSemiBackground()->GetCmpvA().first;
      ssnorm = GetBackground()->GetSemiBackground()->GetCnorm().first;
      ssgw   = GetBackground()->GetSemiBackground()->GetCgwA().first;
      double params[] = {slw,smpv,snorm,sgw,blw,bmpv,bnorm,bgw,sslw,ssmpv,ssnorm,ssgw,bC};
      f->SetParameters(params);
>>>>>>> develop
      fCFit.push_back(f);
      // fCIntegral.push_back(CoherentFitUtils::GetFunctionNormalizationInsideHistogramBoundaries(fh[i],f));
    }
  }
}

//********************************************************************
void CoherentSample::StoreCoherentFitsSignal(){
//********************************************************************

  fClwFit = new TF1("s_lw_par" ,"([0]/x-1)/x+[1]",0,61);
  fClwFit->SetParameters(fClwA.first, fClwB.first);
  
  //fCmpvFit = new TF1("s_mpv_par","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  fCmpvFit = new TF1("s_mpv_par","[0]*(([1]*x-1)/([1]*x+1)+[2])",1,61);
  fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first,fCmpvC.first);
  
  fCgwFit = new TF1("s_gw_par" ,"[0]/(x+1)+[1]",0,61);
  fCgwFit->SetParameters(fCgwA.first, fCgwB.first);

  //store invidual histogram fits
  for(int i = 0; i < (int)fRR.size(); i++){
    std::stringstream ssi;
    ssi << i;
    TF1* f = new TF1(("CF_s_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);
    
    double lw,mpv,norm,gw;
    
    lw   = fClwFit->Eval(fRR[i].first);
    mpv  = fCmpvFit->Eval(fRR[i].first);
    norm = fCnorm.first;
    gw   = fCgwFit->Eval(fRR[i].first);
    
    f->SetParameters(lw,mpv,norm,gw);		       
    fCFit.push_back(f);
  }

  if(fType == SampleTypeEnum::kTrueSignal)StoreCoherentFitsTrueSemiSignal();
}

//********************************************************************
void CoherentSample::StoreCoherentFitsTrueSemiSignal(){
//********************************************************************

  fTrueSemiSignal->fClwFit = new TF1("ss_lw_par" ,"[0]+[1]*x",0,61);
  fTrueSemiSignal->fClwFit->SetParameters(fTrueSemiSignal->fClwA.first, fTrueSemiSignal->fClwB.first);

  fTrueSemiSignal->fCmpvFit = new TF1("ss_mpv_par" ,"[0]+[1]*TMath::DiLog(x+[2])",0,61);
  fTrueSemiSignal->fCmpvFit->SetParameters(fTrueSemiSignal->fCmpvA.first, fTrueSemiSignal->fCmpvB.first, fTrueSemiSignal->fCmpvC.first);

  //store invidual histogram fits
  for(int i = 0; i < (int)fRR.size(); i++){
    std::stringstream ssi;
    ssi << i;
    TF1* f = new TF1(("CF_ss_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);
    
    double lw,mpv,norm,gw;
    
    lw   = fTrueSemiSignal->fClwFit->Eval(fRR[i].first);
    mpv  = fTrueSemiSignal->fCmpvFit->Eval(fRR[i].first);
    norm = fTrueSemiSignal->fCnorm.first;
    gw   = fTrueSemiSignal->fCgwA.first;
    
    f->SetParameters(lw,mpv,norm,gw);		       
    fTrueSemiSignal->fCFit.push_back(f);
  }
}

//********************************************************************
void CoherentSample::StoreCoherentFitsBackground(){
//********************************************************************

  fCmpvFit = new TF1("b_mpw_par","[0]+[1]*x",1,61);
  fCmpvFit->SetParameters(fCmpvA.first,fCmpvB.first);
  
  //store invidual histogram fits
  for(int i = 0; i < (int)fRR.size(); i++){
    std::stringstream ssi;
    ssi << i;
    TF1* f = new TF1(("CF_b_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);
    
    double lw,mpv,norm,gw;
    
    lw   = fClwA.first;
    mpv  = fCmpvFit->Eval(fRR[i].first);
    norm = fCnorm.first;
    gw   = fCgwA.first;
    
    f->SetParameters(lw,mpv,norm,gw);		       
    fCFit.push_back(f);
  }

  if(fType == SampleTypeEnum::kBackground)StoreCoherentFitsSemiBackground();
}

//********************************************************************
void CoherentSample::StoreCoherentFitsSemiBackground(){
//********************************************************************

  fSemiBackground->fClwFit = new TF1("sb_lw_par" ,"[0]/(x+1)+[1]",1,61);
  fSemiBackground->fClwFit->SetParameters(fSemiBackground->fClwA.first, fSemiBackground->fClwB.first);

  //store invidual histogram fits
  for(int i = 0; i < (int)fRR.size(); i++){
    std::stringstream ssi;
    ssi << i;
    TF1* f = new TF1(("CF_sb_"+ssi.str()+"").c_str(),CoherentFitUtils::Langaus,1,30,4);
    
    double lw,mpv,norm,gw;
    
    lw   = fSemiBackground->fClwFit->Eval(fRR[i].first);
    mpv  = fSemiBackground->fCmpvA.first;
    norm = fSemiBackground->fCnorm.first;
    gw   = fSemiBackground->fCgwA.first;
    
    f->SetParameters(lw,mpv,norm,gw);		       
    fSemiBackground->fCFit.push_back(f);
  }
}

//********************************************************************
void CoherentSample::SetCFitParameters(const CoherentSample* sample){
//********************************************************************

<<<<<<< HEAD
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
=======
  if(sample->GetSampleType() == SampleTypeEnum::kTrueSignal){
    GetSignal()->SetClwA(sample->GetClwA());
    GetSignal()->SetClwB(sample->GetClwB());
    GetSignal()->SetClwC(sample->GetClwC());
    GetSignal()->SetCmpvA(sample->GetCmpvA());
    GetSignal()->SetCmpvB(sample->GetCmpvB());
    GetSignal()->SetCmpvC(sample->GetCmpvC());
    GetSignal()->SetCmpvD(sample->GetCmpvD());
    GetSignal()->SetCgwA(sample->GetCgwA());
    GetSignal()->SetCgwB(sample->GetCgwB());
    GetSignal()->SetCgwC(sample->GetCgwC());
    GetSignal()->SetCnorm(sample->GetCnorm());
    GetSignal()->GetSemiSignal()->SetClwA( sample->GetTrueSemiSignal()->GetClwA());
    GetSignal()->GetSemiSignal()->SetClwB( sample->GetTrueSemiSignal()->GetClwB());
    GetSignal()->GetSemiSignal()->SetCmpvA(sample->GetTrueSemiSignal()->GetCmpvA());
    GetSignal()->GetSemiSignal()->SetCmpvB(sample->GetTrueSemiSignal()->GetCmpvB());
    GetSignal()->GetSemiSignal()->SetCmpvC(sample->GetTrueSemiSignal()->GetCmpvC());
    GetSignal()->GetSemiSignal()->SetCgwA( sample->GetTrueSemiSignal()->GetCgwA());
    GetSignal()->GetSemiSignal()->SetCnorm(sample->GetTrueSemiSignal()->GetCnorm());
  }
  else if(sample->GetSampleType() == SampleTypeEnum::kSignal){
    this->SetClwA(sample->GetClwA());
    this->SetClwB(sample->GetClwB());
    this->SetClwC(sample->GetClwC());
    this->SetCmpvA(sample->GetCmpvA());
    this->SetCmpvB(sample->GetCmpvB());
    this->SetCmpvC(sample->GetCmpvC());
    this->SetCmpvD(sample->GetCmpvD());
    this->SetCgwA(sample->GetCgwA());
    this->SetCgwB(sample->GetCgwB());
    this->SetCgwC(sample->GetCgwC());
    this->SetCnorm(sample->GetCnorm());
    this->GetSemiSignal()->SetClwA( sample->GetSemiSignal()->GetClwA());
    this->GetSemiSignal()->SetClwB( sample->GetSemiSignal()->GetClwB());
    this->GetSemiSignal()->SetCmpvA(sample->GetSemiSignal()->GetCmpvA());
    this->GetSemiSignal()->SetCmpvB(sample->GetSemiSignal()->GetCmpvB());
    this->GetSemiSignal()->SetCgwA( sample->GetSemiSignal()->GetCgwA());
    this->GetSemiSignal()->SetCnorm(sample->GetSemiSignal()->GetCnorm());
  }
  else if(sample->GetSampleType() == SampleTypeEnum::kTrueBackground){
    GetBackground()->SetClwA(sample->GetClwA());
    GetBackground()->SetCmpvA(sample->GetCmpvA());
    GetBackground()->SetCmpvB(sample->GetCmpvB());
    GetBackground()->SetCgwA(sample->GetCgwA());
    GetBackground()->SetCnorm(sample->GetCnorm());
    GetBackground()->SetCconst(sample->GetCconst());
    GetBackground()->GetSemiBackground()->SetClwA(sample->GetTrueSemiBackground()->GetClwA());
    GetBackground()->GetSemiBackground()->SetClwB(sample->GetTrueSemiBackground()->GetClwB());
    GetBackground()->GetSemiBackground()->SetCmpvA(sample->GetTrueSemiBackground()->GetCmpvA());
    GetBackground()->GetSemiBackground()->SetCmpvB(sample->GetTrueSemiBackground()->GetCmpvB());
    GetBackground()->GetSemiBackground()->SetCgwA(sample->GetTrueSemiBackground()->GetCgwA());
    GetBackground()->GetSemiBackground()->SetCgwB(sample->GetTrueSemiBackground()->GetCgwB());
    GetBackground()->GetSemiBackground()->SetCnorm(sample->GetTrueSemiBackground()->GetCnorm());
    GetBackground()->GetSemiBackground()->SetCconst(sample->GetTrueSemiBackground()->GetCconst());
  }
  else if(sample->GetSampleType() == SampleTypeEnum::kBackground){
    this->SetClwA(sample->GetClwA());
    this->SetCmpvA(sample->GetCmpvA());
    this->SetCmpvB(sample->GetCmpvB());
    this->SetCgwA(sample->GetCgwA());
    this->SetCnorm(sample->GetCnorm());
    this->SetCconst(sample->GetCconst());
    this->GetSemiBackground()->SetClwA(sample->GetSemiBackground()->GetClwA());
    this->GetSemiBackground()->SetClwB(sample->GetSemiBackground()->GetClwB());
    this->GetSemiBackground()->SetCmpvA(sample->GetSemiBackground()->GetCmpvA());
    this->GetSemiBackground()->SetCmpvB(sample->GetSemiBackground()->GetCmpvB());
    this->GetSemiBackground()->SetCgwA(sample->GetSemiBackground()->GetCgwA());
    this->GetSemiBackground()->SetCgwB(sample->GetSemiBackground()->GetCgwB());
    this->GetSemiBackground()->SetCnorm(sample->GetSemiBackground()->GetCnorm());
    this->GetSemiBackground()->SetCconst(sample->GetSemiBackground()->GetCconst());
  }
}

//********************************************************************
void CoherentSample::SetCFitParametersWithVariations(const CoherentSample* sample, TRandom3* r){
//********************************************************************

  //set parameters with variations
  if(sample->GetSampleType() == SampleTypeEnum::kSignal){
    this->SetClwA(ApplyVariation(sample->GetClwA(),r));
    this->SetClwB(ApplyVariation(sample->GetClwB(),r));
    this->SetCmpvA(ApplyVariation(sample->GetCmpvA(),r));
    this->SetCmpvB(ApplyVariation(sample->GetCmpvB(),r));
    std::cout << sample->GetCmpvC().first << " " << sample->GetCmpvC().second << std::endl;
    this->SetCmpvC(ApplyVariation(sample->GetCmpvC(),r));
    std::cout << this->GetCmpvC().first << " " << this->GetCmpvC().second << std::endl;
    this->SetCgwA(ApplyVariation(sample->GetCgwA(),r));
    this->SetCgwB(ApplyVariation(sample->GetCgwB(),r));
    this->SetCnorm(ApplyVariation(sample->GetCnorm(),r));
    this->GetSemiSignal()->SetClwA( ApplyVariation(sample->GetSemiSignal()->GetClwA(),r));
    this->GetSemiSignal()->SetClwB( ApplyVariation(sample->GetSemiSignal()->GetClwB(),r));
    this->GetSemiSignal()->SetCmpvA(ApplyVariation(sample->GetSemiSignal()->GetCmpvA(),r));
    this->GetSemiSignal()->SetCmpvB(ApplyVariation(sample->GetSemiSignal()->GetCmpvB(),r));
    this->GetSemiSignal()->SetCgwA( ApplyVariation(sample->GetSemiSignal()->GetCgwA(),r));
    this->GetSemiSignal()->SetCnorm(ApplyVariation(sample->GetSemiSignal()->GetCnorm(),r));
  }
  else if(sample->GetSampleType() == SampleTypeEnum::kBackground){
    this->SetClwA(ApplyVariation(sample->GetClwA(),r));
    this->SetCmpvA(ApplyVariation(sample->GetCmpvA(),r));
    this->SetCmpvB(ApplyVariation(sample->GetCmpvB(),r));
    this->SetCgwA(ApplyVariation(sample->GetCgwA(),r));
    this->SetCnorm(ApplyVariation(sample->GetCnorm(),r));
    this->SetCconst(ApplyVariation(sample->GetCconst(),r));
    this->GetSemiBackground()->SetClwA(ApplyVariation(sample->GetSemiBackground()->GetClwA(),r));
    this->GetSemiBackground()->SetClwB(ApplyVariation(sample->GetSemiBackground()->GetClwB(),r));
    this->GetSemiBackground()->SetCmpvA(ApplyVariation(sample->GetSemiBackground()->GetCmpvA(),r));
    this->GetSemiBackground()->SetCmpvB(ApplyVariation(sample->GetSemiBackground()->GetCmpvB(),r));
    this->GetSemiBackground()->SetCgwA(ApplyVariation(sample->GetSemiBackground()->GetCgwA(),r));
    this->GetSemiBackground()->SetCgwB(ApplyVariation(sample->GetSemiBackground()->GetCgwB(),r));
    this->GetSemiBackground()->SetCnorm(ApplyVariation(sample->GetSemiBackground()->GetCnorm(),r));
    this->GetSemiBackground()->SetCconst(ApplyVariation(sample->GetSemiBackground()->GetCconst(),r));
>>>>>>> develop
  }
}

//********************************************************************
std::pair<double,double> CoherentSample::ApplyVariation(const std::pair<double,double> p, TRandom3* r){
//********************************************************************

  return std::make_pair(r->Gaus(p.first,p.second),0);
}

//********************************************************************
void CoherentSample::SequentialCoherentFit(bool minos){
//********************************************************************

  if(fh.empty()){
    std::cout << "No histograms to fit" << std::endl;
    std::exit(1);
  }

  std::cout << "---------------------------------" << std::endl;
  std::cout << "Starting sequential coherent fit" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  IncoherentFit();
  GetInitialParValuesForCoherentFit();
  CoherentFit(minos);
}

//********************************************************************
void CoherentSample::IncoherentFit(){
//********************************************************************

<<<<<<< HEAD
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
=======
  std::cout << "---------------------------------" << std::endl;
  std::cout << "Incoherent fit" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  if(fType == SampleTypeEnum::kTrueSignal)         IncoherentFitSignal();
  else if(fType == SampleTypeEnum::kTrueBackground)IncoherentFitBackground();
  else                                             return;
}

//********************************************************************
void CoherentSample::IncoherentFitSignal(){
//********************************************************************

  //loop over signal histograms
  for(int ihist = 0; ihist < (int)fh.size(); ihist++){
    //double langaus fit
    fIFit.push_back(CoherentFitUtils::DoubleLangausFit(fh[ihist]));

    //save parameters and errors
    double lw1   = fIFit.back()->GetParameter(0);   
    double mpv1  = fIFit.back()->GetParameter(1);
    double norm1 = fIFit.back()->GetParameter(2);
    double gw1   = fIFit.back()->GetParameter(3);
    double fw1   = sqrt(pow(lw1,2)+pow(gw1,2));
    double lw2   = fIFit.back()->GetParameter(4);
    double mpv2  = fIFit.back()->GetParameter(5);
    double norm2 = fIFit.back()->GetParameter(6);
    double gw2   = fIFit.back()->GetParameter(7);
    double fw2   = sqrt(pow(lw2,2)+pow(gw2,2));

    double error_lw1   = fIFit.back()->GetParError(0);
    double error_mpv1  = fIFit.back()->GetParError(1);
    double error_norm1 = fIFit.back()->GetParError(2);
    double error_gw1   = fIFit.back()->GetParError(3);
    double error_fw1   = sqrt(pow(lw1*error_lw1,2)+pow(gw1*error_gw1,2))/fw1;
    double error_lw2   = fIFit.back()->GetParError(4);
    double error_mpv2  = fIFit.back()->GetParError(5);
    double error_norm2 = fIFit.back()->GetParError(6);
    double error_gw2   = fIFit.back()->GetParError(7);
    double error_fw2   = sqrt(pow(lw2*error_lw2,2)+pow(gw2*error_gw2,2))/fw2;
    
    fIlw.push_back(  std::make_pair(lw1   , error_lw1));
    fImpv.push_back( std::make_pair(mpv1  , error_mpv1));
    fInorm.push_back(std::make_pair(norm1 , error_norm1));
    fIgw.push_back(  std::make_pair(gw1   , error_gw1));
    fIfw.push_back(  std::make_pair(fw1   , error_fw1));
    fIalpha.push_back(std::make_pair(gw1/fw1,gw1/fw1*sqrt(pow(error_gw1/gw1,2)+pow(error_fw1/fw1,2))));

    fTrueSemiSignal->AddToIlwVector(  std::make_pair(lw2   , error_lw2));    
    fTrueSemiSignal->AddToImpvVector( std::make_pair(mpv2  , error_mpv2));  
    fTrueSemiSignal->AddToInormVector(std::make_pair(norm2 , error_norm2));
    fTrueSemiSignal->AddToIgwVector(  std::make_pair(gw2   , error_gw2));    
    fTrueSemiSignal->AddToIfwVector(  std::make_pair(fw2   , error_fw2));    
>>>>>>> develop
  }
}

//********************************************************************
void CoherentSample::IncoherentFitBackground(){
//********************************************************************

<<<<<<< HEAD
  std::cout << "Initial par values" << std::endl;

  if(fType == SampleTypeEnum::kTrueSignal || fType == SampleTypeEnum::kTrueBackground){
    CoherentFitUtils::GetABCParametrization(fClwA.first,fClwB.first,fClwC.first,
					    fRR,fIlw,true,draw_fits);
    CoherentFitUtils::GetABCParametrization(fCgwA.first,fCgwB.first,fCgwC.first,
					    fRR,fIgw,equal_weights,draw_fits);
    CoherentFitUtils::GetABCDRParametrization(fCmpvA.first,fCmpvB.first,fCmpvC.first,fCmpvD.first,fCmpvR.first,
					      fRR,fImpv,equal_weights,draw_fits);
=======
  //loop over signal histograms
  for(int ihist = 0; ihist < (int)fh.size(); ihist++){
    //double langaus fit
    fIFit.push_back(CoherentFitUtils::LangausPlusConstantFit(fh[ihist]));
    
    //save parameters and errors
    double lw   = fIFit.back()->GetParameter(0);   
    double mpv  = fIFit.back()->GetParameter(1);
    double norm = fIFit.back()->GetParameter(2);
    double gw   = fIFit.back()->GetParameter(3);
    double fw   = sqrt(pow(lw,2)+pow(gw,2));
    double c    = fIFit.back()->GetParameter(4);

    double error_lw   = fIFit.back()->GetParError(0);
    double error_mpv  = fIFit.back()->GetParError(1);
    double error_norm = fIFit.back()->GetParError(2);
    double error_gw   = fIFit.back()->GetParError(3);
    double error_fw   = sqrt(pow(lw*error_lw,2)+pow(gw*error_gw,2))/fw;
    double error_c    = fIFit.back()->GetParError(4);
    
    fIlw.push_back(  std::make_pair(lw   , error_lw));
    fImpv.push_back( std::make_pair(mpv  , error_mpv));
    fInorm.push_back(std::make_pair(norm , error_norm));
    fIgw.push_back(  std::make_pair(gw   , error_gw));
    fIfw.push_back(  std::make_pair(fw   , error_fw));
    fIconst.push_back(  std::make_pair(c , error_c));    
>>>>>>> develop
  }
}

//********************************************************************
void CoherentSample::GetInitialParValuesForCoherentFit(){
//********************************************************************

  std::cout << "---------------------------------" << std::endl;
  std::cout << "Estimating initial coherent pars" << std::endl;
  std::cout << "---------------------------------" << std::endl;
  
  if(fType == SampleTypeEnum::kTrueSignal)
    GetInitialParValuesForSignal();
  else if(fType == SampleTypeEnum::kTrueBackground)
    GetInitialParValuesForBackground();
  else if(fType == SampleTypeEnum::kSignalPlusBackground){
    SetCFitParameters(GetTrueSignal());
    SetCFitParameters(GetTrueBackground());
    if(fTrueSignal && fTrueBackground)EstimateWhatIsLeft();
  }
}

//********************************************************************
void CoherentSample::GetInitialParValuesForSignal(){
//********************************************************************

<<<<<<< HEAD
  std::cout << "coherent fit" << std::endl;

  if(fType == SampleTypeEnum::kTrueSignal){
    CoherentFitSignal();
    //CoherentFitSignalCheb();
=======
  //signal lw
  double lwA1, lwB1;
  double error_lwA1, error_lwB1;
  CoherentFitUtils::GetSignalLWEstimation(lwA1,lwB1,
					  error_lwA1,error_lwB1,
					  fRR, fIlw);
  fClwA = std::make_pair(lwA1,error_lwA1);
  fClwB = std::make_pair(lwB1,error_lwB1);

  //signal mpv
  double mpvA1, mpvB1, mpvC1;
  double error_mpvA1, error_mpvB1, error_mpvC1;
  CoherentFitUtils::GetSignalMPVEstimation(mpvA1,mpvB1,mpvC1,
					   error_mpvA1,error_mpvB1,error_mpvC1,
					   fRR, fImpv);
  fCmpvA = std::make_pair(mpvA1,error_mpvA1);
  fCmpvB = std::make_pair(mpvB1,error_mpvB1);
  fCmpvC = std::make_pair(mpvC1,error_mpvC1);

  //signal gw
  double gwA1, gwB1;
  double error_gwA1, error_gwB1;
  CoherentFitUtils::GetSignalGWEstimation(gwA1,gwB1,
					  error_gwA1,error_gwB1,
					  fRR, fIgw);
  fCgwA = std::make_pair(gwA1,error_gwA1);
  fCgwB = std::make_pair(gwB1,error_gwB1);

  //signal norm
  double mean = 0, error_mean = 0;
  int counter = 0;
  for(int i = 0; i < (int)fInorm.size(); i++){
    mean += fInorm[i].first;
    error_mean += fInorm[i].second;
    counter++;
  }
  mean /= counter;
  error_mean /= counter;
  fCnorm = std::make_pair(mean,error_mean);

  //semisignal lw
  double lwA2, lwB2;
  double error_lwA2, error_lwB2;
  CoherentFitUtils::GetLinearParametrization(lwA2,lwB2,
  					     error_lwA2,error_lwB2,
  					     fTrueSemiSignal->fImpv, fTrueSemiSignal->fIlw);
  fTrueSemiSignal->fClwA = std::make_pair(lwA2,error_lwA2);
  fTrueSemiSignal->fClwB = std::make_pair(lwB2,error_lwB2);
  
  //semisignal mpv
  double mpvA2, mpvB2, mpvC2;
  double error_mpvA2, error_mpvB2, error_mpvC2;
  CoherentFitUtils::GetDiLogParametrization(mpvA2,mpvB2,mpvC2,
  					     error_mpvA2,error_mpvB2,error_mpvC2,
  					     fTrueSemiSignal->fRR, fTrueSemiSignal->fImpv);
  fTrueSemiSignal->fCmpvA = std::make_pair(mpvA2,error_mpvA2);
  fTrueSemiSignal->fCmpvB = std::make_pair(mpvB2,error_mpvB2);
  fTrueSemiSignal->fCmpvC = std::make_pair(mpvC2,error_mpvC2);

  //semisignal gw
  mean = 0;
  error_mean = 0;
  counter = 0;
  for(int i = 0; i < (int)fTrueSemiSignal->fInorm.size(); i++){
    mean += fTrueSemiSignal->fIgw[i].first;
    error_mean += fTrueSemiSignal->fIgw[i].second;
    counter++;
>>>>>>> develop
  }
  mean /= counter;
  error_mean /= counter;
  fTrueSemiSignal->fCgwA = std::make_pair(mean,error_mean);

  //semisignal norm
  mean = 0;
  error_mean = 0;
  counter = 0;
  for(int i = 0; i < (int)fTrueSemiSignal->fInorm.size(); i++){
    mean += fTrueSemiSignal->fInorm[i].first;
    error_mean += fTrueSemiSignal->fInorm[i].second;
    counter++;
  }
  mean /= counter;
  error_mean /= counter;
  fTrueSemiSignal->fCnorm = std::make_pair(mean,error_mean);
}

//********************************************************************
void CoherentSample::GetInitialParValuesForBackground(){
//********************************************************************

  //signal lw
  double lwA1, lwB1;
  double error_lwA1, error_lwB1;
  CoherentFitUtils::GetLinearParametrization(lwA1,lwB1,
					     error_lwA1,error_lwB1,
					     fRR, fIlw);
  fClwA = std::make_pair(lwA1,error_lwA1);
  fClwB = std::make_pair(lwB1,error_lwB1);

  //signal mpv
  double mpvA1, mpvB1;
  double error_mpvA1, error_mpvB1;
  CoherentFitUtils::GetLinearParametrization(mpvA1,mpvB1,
					     error_mpvA1,error_mpvB1,
					     fRR, fImpv);
  fCmpvA = std::make_pair(mpvA1,error_mpvA1);
  fCmpvB = std::make_pair(mpvB1,error_mpvB1);

  //signal gw
  double mean = 0, error_mean = 0;
  int counter = 0;
  for(int i = 0; i < (int)fIgw.size(); i++){
    mean += fIgw[i].first;
    error_mean += fIgw[i].second;
    counter++;
  }
  mean /= counter;
  error_mean /= counter;
  fCgwA = std::make_pair(mean,error_mean);

  //signal norm
  mean = 0;
  error_mean = 0;
  counter = 0;
  for(int i = 0; i < (int)fInorm.size(); i++){
    mean += fInorm[i].first;
    error_mean += fInorm[i].second;
    counter++;
  }
  mean /= counter;
  error_mean /= counter;
  fCnorm = std::make_pair(mean,error_mean);

  //const
  mean = 0;
  error_mean = 0;
  counter = 0;
  for(int i = 0; i < (int)fIconst.size(); i++){
    mean += fIconst[i].first;
    error_mean += fIconst[i].second;
    counter++;
  }
  mean /= counter;
  error_mean /= counter;
  fCconst = std::make_pair(mean,error_mean);
}

//********************************************************************
void CoherentSample::CoherentFit(bool minos){
//********************************************************************

  std::cout << "---------------------------------" << std::endl;
  std::cout << "Coherent fit" << std::endl;
  std::cout << "---------------------------------" << std::endl;
  
  if(fType == SampleTypeEnum::kTrueSignal)
    CoherentFitSignal();
  else if(fType == SampleTypeEnum::kTrueBackground)
    CoherentFitBackground();
  else if(fType == SampleTypeEnum::kSignalPlusBackground)
    CoherentFitSignalPlusBackground(minos);
}

//********************************************************************
TGraphErrors* CoherentSample::GetMPVErrorBandOld(){
//********************************************************************

  //get covariance matrix
  // const int ndim = fMinuit->GetNumPars();
  // double matrix[ndim][ndim] = {0};

  int ndim = fMinuit->GetNumPars();  // No need for 'const'
  std::vector<std::vector<double>> matrix(ndim, std::vector<double>(ndim, 0.0));

  fMinuit->mnemat(&matrix[0][0],ndim);

  std::cout << matrix[2][2] << " " << matrix[3][3] << " " << matrix[4][4] << std::endl; 
  std::cout << fSignal->fCmpvA.first*0.036 << " " << fSignal->fCmpvB.first*0.054 << " " << fSignal->fCmpvC.first*0.036 << std::endl;

  matrix[2][2] = matrix[2][2]+pow(fSignal->fCmpvA.first*0.036,2);
  matrix[3][3] = matrix[3][3]+pow(fSignal->fCmpvB.first*0.054,2);
  matrix[4][4] = matrix[4][4]+pow(fSignal->fCmpvC.first*0.036,2);

  std::cout << matrix[2][2] << " " << matrix[3][3] << " " << matrix[4][4] << std::endl; 
  
  //get hessian
<<<<<<< HEAD
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

=======
  TF1* J[3];
  J[0] = new TF1("fa","([1]*x-1)/([1]*x+1)+[0]*[2]*0",1,61);
  J[1] = new TF1("fb","[0]*2*x/pow([1]*x+1,2)+[2]*0",1,61);
  J[2] = new TF1("fc","1+[0]*[1]*[2]*0",1,61);
  for(int i = 0; i < 3; i++)J[i]->SetParameters(fSignal->fCmpvA.first,fSignal->fCmpvB.first,fSignal->fCmpvC.first);
  TF1* mpvFit = fSignal->GetCmpvFit();
  
>>>>>>> develop
  //compute error at each RR
  std::vector<double> rr,rr_error,mpv,mpv_error;
  rr.clear();
  rr_error.clear();
  mpv.clear();
  mpv_error.clear();
  double error    = 0;
  double aux      = 0;
  double v_aux[3] = {0};
  for(int irr = 0; irr < (int)fRR.size(); irr++){
    rr.push_back(fRR[irr].first);
    rr_error.push_back(fRR[irr].second);
    mpv.push_back(mpvFit->Eval(fRR[irr].first));
    for(int icolumn = 2; icolumn < 5; icolumn++){
      for(int irow = 2; irow < 5; irow++)
	aux += J[irow-2]->Eval(fRR[irr].first)*matrix[irow][icolumn];
      v_aux[icolumn-2] = aux;
      aux = 0;
    }
    for(int i = 0; i < 3; i++)error += v_aux[i]*J[i]->Eval(fRR[irr].first);
    mpv_error.push_back(sqrt(error));
    error = 0;
  }

  return new TGraphErrors(mpv.size(),&rr[0],&mpv[0],&rr_error[0],&mpv_error[0]);
}

//********************************************************************
TGraphAsymmErrors* CoherentSample::GetMPVErrorBand(const double min, const double max, const double step){
//********************************************************************

  std::vector<double> rr, val, low_error, high_error;
  double A   = fSignal->fCmpvA.first;
  double A_l = fCmpvA_error_neg;
  double A_h = fCmpvA_error_pos;
  double B   = fSignal->fCmpvB.first;
  double B_l = fCmpvB_error_neg;
  double B_h = fCmpvB_error_pos;
  double C   = fSignal->fCmpvC.first;
  double C_l = fCmpvC_error_neg;
  double C_h = fCmpvC_error_pos;

  // std::cout << A+A_l << " " << A << " " << A+A_h << std::endl;
  // std::cout << B+B_l << " " << B << " " << B+B_h << std::endl;
  // std::cout << C+C_l << " " << C << " " << C+C_h << std::endl;
  
  // TF1* fc = new TF1("fc","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  // TF1* fl = new TF1("fl","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  // TF1* fh = new TF1("fh","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  
  // fc->SetParameters(A,B,C);
  // fl->SetParameters(A+A_l,B+B_l,C+C_l);
  // fh->SetParameters(A+A_h,B+B_h,C+C_h);

  // int n = (max-min)/step;

  // for(int i = 0; i < n; i++){
  //   rr.push_back(min+i*step);
  //   val.push_back(fc->Eval(min+i*step));
  //   low_error.push_back(fc->Eval(min+i*step)-fl->Eval(min+i*step));
  //   high_error.push_back(fh->Eval(min+i*step)-fc->Eval(min+i*step));
  // }

  TF1* fc = new TF1("fc","[0]*([1]*x-1)/([1]*x+1)+[2]",1,61);
  fc->SetParameters(A,B,C);

  double A_error = (abs(A_l)+abs(A_h))/2;
  double B_error = (abs(B_l)+abs(B_h))/2;
  double C_error = (abs(C_l)+abs(C_h))/2;

  int n = (max-min)/step;
  for(int i = 0; i < n; i++){
    double x = min+i*step;
    rr.push_back(x);
    val.push_back(fc->Eval(x));
    double error = sqrt(pow(A_error*(B*x-1)/(B*x+1),2)+pow(2*x*A*B_error/pow(B*x+1,2),2)+pow(C_error,2));
    low_error.push_back(error);
    high_error.push_back(error);
  }

  return new TGraphAsymmErrors(rr.size(),&rr[0],&val[0],0,0,&low_error[0],&high_error[0]);
}

//********************************************************************
TRatioPlot* CoherentSample::GetFitRatio(const int ibin){
//********************************************************************

  TH1F* hdummy = (TH1F*)fh[ibin]->Clone("hdummy");
  hdummy->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  for(int i = 0; i < 13; i++)fCFit[ibin]->FixParameter(i,fCFit[ibin]->GetParameter(i));

  std::string name = fCFit[ibin]->GetTitle();
  hdummy->Fit(name.c_str(),"");

  return new TRatioPlot(hdummy,"errasym");
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
TGraphErrors* CoherentSample::GetIconstGraph() const {
//********************************************************************

  return CoherentFitUtils::GetGraph(fRR,fIconst);
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

//********************************************************************
void CoherentSample::CoherentFitSignal(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 15;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignal);

<<<<<<< HEAD
  Double_t arglist[10];
=======
  Double_t arglist[11];
>>>>>>> develop
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
<<<<<<< HEAD
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
  fMinuit->mnparm(11, "norm",  mean_norm, 0.01, 0, 1, ierflg);

  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
=======
  double lwA,lwB,lwC,mpvA,mpvB,mpvC,mpvD,gwA,gwB,gwC,norm;
  double lwA_error,lwB_error,lwC_error,mpvA_error,mpvB_error,mpvC_error,mpvD_error,gwA_error,gwB_error,gwC_error,norm_error;
  double sslwA,sslwB,ssmpvA,ssmpvB,ssmpvC,ssgw,ssnorm;
  double sslwA_error,sslwB_error,ssmpvA_error,ssmpvB_error,ssmpvC_error,ssgw_error,ssnorm_error;
    
  mpvA = fCmpvA.first;
  mpvB = fCmpvB.first;
  mpvC = fCmpvC.first;

  gwA = fCgwA.first;
  gwB = fCgwB.first;
  
  lwA = fClwA.first;
  lwB = fClwB.first;

  norm = fCnorm.first;

  ssmpvA = fTrueSemiSignal->fCmpvA.first;
  ssmpvB = fTrueSemiSignal->fCmpvB.first;
  ssmpvC = fTrueSemiSignal->fCmpvC.first;

  sslwA = fTrueSemiSignal->fClwA.first;
  sslwB = fTrueSemiSignal->fClwB.first;

  ssgw = fTrueSemiSignal->fCgwA.first;

  ssnorm = fTrueSemiSignal->fCnorm.first;

  fMinuit->mnparm(0,  "smpv A", mpvA   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(1,  "smpv B", mpvB   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(2,  "smpv C", mpvC   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(3,  "lwA"   , lwA    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(4,  "lwB"   , lwB    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(5,  "gwA"   , gwA    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(6,  "gwB"   , gwB    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(7,  "snorm" , norm   , 0.001, 0, 1, ierflg);
  fMinuit->mnparm(8,  "bmpv A", ssmpvA , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(9,  "bmpv B", ssmpvB , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(10, "bmpv C", ssmpvC , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(11, "blwA"  , sslwA  , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(12, "blwB"  , sslwB  , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(13, "bgw"   , ssgw   , 0.1, 0, 0, ierflg);
  fMinuit->mnparm(14, "bnorm" , ssnorm , 0.01, 0, 1, ierflg);

  
  arglist[0] = 500000;
  arglist[1] = 1;
  fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
>>>>>>> develop
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  //retrieve parameters for next stage
<<<<<<< HEAD
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
=======
  fMinuit->GetParameter(0,  fCmpvA.first, fCmpvA.second);
  fMinuit->GetParameter(1,  fCmpvB.first, fCmpvB.second);
  fMinuit->GetParameter(2,  fCmpvC.first, fCmpvC.second);
  fMinuit->GetParameter(3, fClwA.first,  fClwA.second);
  fMinuit->GetParameter(4, fClwB.first,  fClwB.second);
  fMinuit->GetParameter(5, fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(6, fCgwB.first,  fCgwB.second);
  fMinuit->GetParameter(7, fCnorm.first, fCnorm.second);

  fMinuit->GetParameter(8,  fTrueSemiSignal->fCmpvA.first, fTrueSemiSignal->fCmpvA.second);
  fMinuit->GetParameter(9,  fTrueSemiSignal->fCmpvB.first, fTrueSemiSignal->fCmpvB.second);
  fMinuit->GetParameter(10, fTrueSemiSignal->fCmpvC.first, fTrueSemiSignal->fCmpvC.second);
  fMinuit->GetParameter(11, fTrueSemiSignal->fClwA.first,  fTrueSemiSignal->fClwA.second);
  fMinuit->GetParameter(12, fTrueSemiSignal->fClwB.first,  fTrueSemiSignal->fClwB.second);
  fMinuit->GetParameter(13, fTrueSemiSignal->fCgwA.first,  fTrueSemiSignal->fCgwA.second);
  fMinuit->GetParameter(14, fTrueSemiSignal->fCnorm.first, fTrueSemiSignal->fCnorm.second);
>>>>>>> develop

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnSignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
<<<<<<< HEAD
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
=======
  double mpvA  = par[0];
  double mpvB  = par[1];
  double mpvC  = par[2];
  double lwA   = par[3];
  double lwB   = par[4];
  double gwA   = par[5];
  double gwB   = par[6];
  double norm  = par[7];
  double bmpvA = par[8];
  double bmpvB = par[9];
  double bmpvC = par[10];
  double blwA  = par[11];
  double blwB  = par[12];
  double bgw   = par[13];
  double bnorm = par[14];
>>>>>>> develop

  //define function to fit to each histogram
  double mpv,fw,gw,lw,bmpv,blw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangaus,1,30,8);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
<<<<<<< HEAD
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
=======
    double x   = CSMinuit->fRR[ibin].first;
    //mpv  = mpvA*(mpvB*x-1)/(mpvB*x+1)+mpvC;
    mpv  = mpvA*((mpvB*x-1)/(mpvB*x+1)+mpvC);
    gw   = gwA/(x+1)+gwB;
    lw   = (lwA/x-1)/x+lwB;
>>>>>>> develop

    bmpv = bmpvA+bmpvB*TMath::DiLog(x+bmpvC);
    blw  = blwA+blwB*bmpv;

    flg->SetParameters(lw,mpv,norm,gw,blw,bmpv,bnorm,bgw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
    // CSMinuit->fh[ibin]->Draw();
    // flg->Draw("same");
    // gPad->Update();gPad->WaitPrimitive();
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitBackground(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 6;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnBackground);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
<<<<<<< HEAD
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
=======
  double lw, mpvA, mpvB, gw, norm, C;
  double lw_error, mpvA_error, mpvB_error, norm_error, gw_error, C_error;

  lw   = fClwA.first;
  mpvA = fCmpvA.first;
  mpvB = fCmpvB.first;
  gw   = fCgwA.first;
  norm = fCnorm.first;
  C    = fCconst.first;

  fMinuit->mnparm(0,  "lw"    , lw   , 0.001  , 0 , 0   , ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(1,  "mpvA"  , mpvA , 0.001   , 0 , 0   , ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(2,  "mpvB"  , mpvB , 0.001   , 0 , 0   , ierflg);//0.9*0.65  , 1.1*0.65  , ierflg);//
  fMinuit->mnparm(3,  "gw"    , gw   , 0.001  , 0 , 0   , ierflg);//0.9*2.25  , 1.1*2.25  , ierflg);//
  fMinuit->mnparm(4,  "norm"  , norm , 0.001  , 0 , 1   , ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//
  fMinuit->mnparm(5,  "const" , C    , 0.001  , 0 , 1   , ierflg);//0.9*0.88  , 1.1*0.88  , ierflg);//
>>>>>>> develop

  // Now ready for minimization step
  arglist[0] = 2;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 500000;
  arglist[1] = 0.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fClwA.first   , fClwA.second);
  fMinuit->GetParameter(1,  fCmpvA.first  , fCmpvA.second);
  fMinuit->GetParameter(2,  fCmpvB.first  , fCmpvB.second);
  fMinuit->GetParameter(3,  fCgwA.first   , fCgwA.second);
  fMinuit->GetParameter(4,  fCnorm.first  , fCnorm.second);
  fMinuit->GetParameter(5,  fCconst.first , fCconst.second);
  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lw    = par[0];
  double mpvA  = par[1];
  double mpvB  = par[2];
  double gw    = par[3];
  double norm  = par[4];
  double C     = par[5];
  
  //define function to fit to each histogram
  double Likelihood = 0;
  double mpv;
  TF1* flg = new TF1("flg",CoherentFitUtils::LangausPlusConstant,1,30,5);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    double x = CSMinuit->fRR[ibin].first;
    mpv = mpvA + mpvB*x;
    flg->SetParameters(lw,mpv,norm,gw,C);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackground(bool minos, double step){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 19;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackground);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  double slwA,slwB,smpvA,smpvB,smpvC,sgwA,sgwB,snorm;
  double slwA_error,slwB_error,smpvA_error,smpvB_error,smpvC_error,sgwA_error,sgwB_error,snorm_error;
  double blwA, blwB, bmpvA, bmpvB, bgw, bnorm, bC;
  double blwA_error, blwB_error, bmpvA_error, bmpvB_error, bnorm_error, bgw_error, bC_error;
  double sslwA,sslwB,ssmpvA,ssmpvB,ssmpvC,ssgw,ssnorm;
  double sslwA_error,sslwB_error,ssmpvA_error,ssmpvB_error,ssmpvC_error,ssgw_error,ssnorm_error;

  slwA  = fSignal->fClwA.first;  
  slwB  = fSignal->fClwB.first;	 
  smpvA = fSignal->fCmpvA.first;
  smpvB = fSignal->fCmpvB.first; 
  smpvC = fSignal->fCmpvC.first; 
  sgwA  = fSignal->fCgwA.first;	 
  sgwB  = fSignal->fCgwB.first;	 
  snorm = fSignal->fCnorm.first;
  
  blwA  = fBackground->fClwA.first;
  bmpvA = fBackground->fCmpvA.first;
  bmpvB = fBackground->fCmpvB.first;
  bgw   = fBackground->fCgwA.first;
  bnorm = fBackground->fCnorm.first;
  bC    = 0;//fBackground->fCconst.first;
  	  
  sslwA  = fBackground->fSemiBackground->fClwA.first;  
  sslwB  = fBackground->fSemiBackground->fClwB.first;	 
  ssmpvA = fBackground->fSemiBackground->fCmpvA.first;
  ssgw   = fBackground->fSemiBackground->fCgwA.first;	 
  ssnorm = fBackground->fSemiBackground->fCnorm.first;

  fMinuit->mnparm(0,  "slw A"   , slwA     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(1,  "slw A"   , slwB     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(2,  "smpv A"  , smpvA    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(3,  "smpv B"  , smpvB    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(4,  "smpv C"  , smpvC    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(5,  "sgw A"   , sgwA     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(6,  "sgw B"   , sgwB     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(7,  "snorm"   , snorm    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(8,  "blwA"    , blwA     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(9,  "bmpvA"   , bmpvA    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(10, "bmpvB"   , bmpvB    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(11, "bgw"     , bgw      , step, 0   , 0   , ierflg);
  fMinuit->mnparm(12, "bnorm"   , bnorm    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(13, "bconst"  , bC       , step, 0   , 1   , ierflg);
  fMinuit->mnparm(14, "sslw A"  , sslwA    , step, 0   , 0   , ierflg);
  fMinuit->mnparm(15, "sslw B"  , sslwB    , step, 0   , 2   , ierflg);
  fMinuit->mnparm(16, "ssmpv A" , ssmpvA   , step, 0   , 0   , ierflg);
  fMinuit->mnparm(17, "ssgw"    , ssgw     , step, 0   , 0   , ierflg);
  fMinuit->mnparm(18, "ssnorm"  , ssnorm   , step, 0   , 0   , ierflg);
  fMinuit->FixParameter(13);
  fMinuit->FixParameter(18);
  
  // Now ready for minimization step
  arglist[0] = 1;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 500000;
  arglist[1] = 5.;
  fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  if(minos){
    // // fMinuit->mnfree(0);
    // arglist[0] = 2;
    // fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
    // arglist[0] = 500000;
    // arglist[1] = 1.;
    fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
    // fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
    arglist[0] = 500000;
    arglist[1] = 3;
    arglist[2] = 4;
    arglist[3] = 5;
    fMinuit->mnexcm("MINOS", arglist, 4, ierflg);
  }
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  //retrieve parameters for next stage

  //SIGNAL
  fMinuit->GetParameter(0 , slwA , slwA_error );
  fMinuit->GetParameter(1 , slwB , slwB_error );
  fMinuit->GetParameter(2 , smpvA, smpvA_error);
  fMinuit->GetParameter(3 , smpvB, smpvB_error);
  fMinuit->GetParameter(4 , smpvC, smpvC_error);
  fMinuit->GetParameter(5 , sgwA , sgwA_error );
  fMinuit->GetParameter(6 , sgwB , sgwB_error );
  fMinuit->GetParameter(7, snorm, snorm_error);

  if(minos){
    double dummy;
    fMinuit->mnerrs(2,fCmpvA_error_pos,fCmpvA_error_neg,smpvA_error,dummy);
    fMinuit->mnerrs(3,fCmpvB_error_pos,fCmpvB_error_neg,smpvB_error,dummy);
    fMinuit->mnerrs(4,fCmpvC_error_pos,fCmpvC_error_neg,smpvC_error,dummy);
    smpvA_error = (fCmpvA_error_pos-fCmpvA_error_neg)/2;
    smpvB_error = (fCmpvB_error_pos-fCmpvB_error_neg)/2;
    smpvC_error = (fCmpvC_error_pos-fCmpvC_error_neg)/2;
  }

  fSignal->SetClwA(std::make_pair(slwA  ,slwA_error ));
  fSignal->SetClwB(std::make_pair(slwB  ,slwB_error ));
  fSignal->SetCmpvA(std::make_pair(smpvA,smpvA_error));
  fSignal->SetCmpvB(std::make_pair(smpvB,smpvB_error));
  fSignal->SetCmpvC(std::make_pair(smpvC,smpvC_error));
  fSignal->SetCgwA(std::make_pair(sgwA  ,sgwA_error ));
  fSignal->SetCgwB(std::make_pair(sgwB  ,sgwB_error ));
  fSignal->SetCnorm(std::make_pair(snorm,snorm_error));

  //BACKGROUND
  fMinuit->GetParameter(8 , blwA  , blwA_error);
  fMinuit->GetParameter(9 , bmpvA,bmpvA_error);
  fMinuit->GetParameter(10 , bmpvB,bmpvB_error);
  fMinuit->GetParameter(11, bgw  , bgw_error);
  fMinuit->GetParameter(12, bnorm, bnorm_error);
  fMinuit->GetParameter(13, bC   , bC_error);

  fBackground->SetClwA(std::make_pair(blwA   ,blwA_error));
  fBackground->SetCmpvA(std::make_pair(bmpvA ,bmpvA_error));
  fBackground->SetCmpvB(std::make_pair(bmpvB ,bmpvB_error));
  fBackground->SetCgwA(std::make_pair(bgw   ,bgw_error));
  fBackground->SetCnorm(std::make_pair(bnorm,bnorm_error));
  fBackground->SetCconst(std::make_pair(bC  ,bC_error));
  
  fMinuit->GetParameter(14, sslwA  ,sslwA_error);
  fMinuit->GetParameter(15, sslwB  ,sslwB_error);
  fMinuit->GetParameter(16, ssmpvA  ,ssmpvA_error);
  fMinuit->GetParameter(17, ssgw   ,ssgw_error);
  fMinuit->GetParameter(18, ssnorm ,ssnorm_error);

  fBackground->fSemiBackground->SetClwA(std::make_pair(sslwA  ,sslwA_error));
  fBackground->fSemiBackground->SetClwB(std::make_pair(sslwB  ,sslwB_error));
  fBackground->fSemiBackground->SetCmpvA(std::make_pair(ssmpvA ,ssmpvA_error));
  fBackground->fSemiBackground->SetCgwA(std::make_pair(ssgw,ssgw_error));
  //fBackground->fSemiBackground->SetCnorm(std::make_pair(ssnorm ,ssnorm_error));
  fBackground->fSemiBackground->SetCnorm(std::make_pair(1-snorm-bnorm ,ssnorm_error));

  //Print covariance matrix
  fMinuit->mnmatu(1);

  //get correlation matrix
  double cov[18][18] = {0};
  fMinuit->mnemat(&cov[0][0],18);
  for(int i = 0; i < 18; i++)
    for(int j = 0; j < 18; j++)
      fcorr_matrix[i][j] = cov[i][j]/sqrt(cov[i][i]*cov[j][j]);
  
  StoreCoherentFits();  
}

//********************************************************************
void CoherentSample::fcnSignalPlusBackground(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lwA   = par[0];
  double lwB   = par[1];
  double mpvA  = par[2];
  double mpvB  = par[3];
  double mpvC  = par[4];
  double gwA   = par[5];
  double gwB   = par[6];
  double norm  = par[7];
  double blwA  = par[8];
  double bmpvA = par[9];
  double bmpvB = par[10];
  double bgw   = par[11];
  double bnorm = par[12];
  double bC    = par[13];
  double sslwA  = par[14];
  double sslwB  = par[15];
  double ssmpv  = par[16];
  double ssgw   = par[17];
  double ssnorm = par[18];

  //set likelihood vector size
  if(CSMinuit->fPartialLikelihood.empty())
    CSMinuit->fPartialLikelihood.resize(CSMinuit->GetSize());
  
  //define function to fit to each histogram
  double lw,mpv,gw,bmpv,blw,sslw;
  double Likelihood = 0;
  double PartialLikelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::TripleLangausPlusConstant,1,30,13);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = CSMinuit->fRR[ibin].first;
    //mpv  = mpvA*(mpvB*x-1)/(mpvB*x+1)+mpvC;
    mpv  = mpvA*((mpvB*x-1)/(mpvB*x+1)+mpvC);
    gw   = gwA/(x+1)+gwB;
    lw   = (lwA/x-1)/x+lwB;

    blw  = blwA;
    bmpv = bmpvA+bmpvB*x;

    sslw = sslwA/(x+1)+sslwB;

    ssnorm = 1-norm-bnorm;

    double params[] = {lw,mpv,norm,gw,blw,bmpv,bnorm,bgw,sslw,ssmpv,ssnorm,ssgw,bC};
    //double params[] = {lw,mpv,norm,gw,blw,bmpv,bnorm,bgw,sslw,ssmpv,1-norm-bnorm,ssgw,bC};
    flg->SetParameters(params);

    PartialLikelihood = CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
    CSMinuit->fPartialLikelihood[ibin] = PartialLikelihood;
    Likelihood = Likelihood + PartialLikelihood;
    if(debug){
      CSMinuit->fh[ibin]->Draw();
      flg->Draw("same");
      gPad->Update();gPad->WaitPrimitive();
    }
  }

  CSMinuit->fLikelihood = Likelihood;
  for(int i = 0; i < CSMinuit->GetSize(); i++)CSMinuit->fPartialLikelihood[i] /= CSMinuit->fLikelihood;

  f = -Likelihood;
}


//********************************************************************
void CoherentSample::EstimateWhatIsLeft(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 19;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
<<<<<<< HEAD
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
=======
  fMinuit->SetFCN(fcnEstimateWhatIsLeft);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  double slwA,slwB,smpvA,smpvB,smpvC,sgwA,sgwB,snorm;
  double blwA, blwB, bmpvA, bmpvB, bgw, bnorm, bC;

  slwA  = fTrueSignal->fClwA.first;  
  slwB  = fTrueSignal->fClwB.first;	 
  smpvA = fTrueSignal->fCmpvA.first;
  smpvB = fTrueSignal->fCmpvB.first; 
  smpvC = fTrueSignal->fCmpvC.first; 
  sgwA  = fTrueSignal->fCgwA.first;	 
  sgwB  = fTrueSignal->fCgwB.first;	 
  snorm = fTrueSignal->fCnorm.first;

  blwA  = fTrueBackground->fClwA.first;
  blwB  = fTrueBackground->fClwA.first;
  bmpvA = fTrueBackground->fCmpvA.first;
  bmpvB = fTrueBackground->fCmpvB.first;
  bgw   = fTrueBackground->fCgwA.first;
  bnorm = fTrueBackground->fCnorm.first;
  bC    = 0;//fTrueBackground->fCconst.first;

  //define parameters and step
  fMinuit->mnparm(0,  "slw A"   , slwA      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(1,  "slw B"   , slwB      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(2,  "smpv A"  , smpvA     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(3,  "smpv B"  , smpvB     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(4,  "smpv C"  , smpvC     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(5,  "sgw A"   , sgwA      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(6,  "sgw B"   , sgwB      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(7,  "snorm"   , snorm     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(8,  "blwA"    , blwA      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(9, "bmpvA"    , bmpvA     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(10, "bmpvB"   , bmpvB     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(11, "bgw"     , bgw       , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(12, "bnorm"   , bnorm     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(13, "bconst"  , bC        , 0.1, 0   , 1   , ierflg);
  fMinuit->mnparm(14, "sslw A"  , 2.2       , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(15, "sslw B"  , 0.09      , 0.1, 0   , 2   , ierflg);
  fMinuit->mnparm(16, "ssmpv A" , 2.4       , 0.1, 0   , 2.2 , ierflg);
  fMinuit->mnparm(17, "ssgw"    , 0.43      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(18, "ssnorm"  , 0.0945    , 0.1, 0   , 1-snorm-bnorm, ierflg);

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
  fMinuit->FixParameter(11);
  fMinuit->FixParameter(12);
  fMinuit->FixParameter(13);
>>>>>>> develop

  // Now ready for minimization step
  arglist[0] = 0;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 500000;
  arglist[1] = 1.;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);
<<<<<<< HEAD

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
=======

  double sslwA,sslwB,ssmpvA,ssmpvB,ssmpvC,ssgw,ssnorm;
  double sslwA_error,sslwB_error,ssmpvA_error,ssmpvB_error,ssmpvC_error,ssgw_error,ssnorm_error;

  fMinuit->GetParameter(14 , sslwA  , sslwA_error  );
  fMinuit->GetParameter(15 , sslwB  , sslwB_error  );
  fMinuit->GetParameter(16 , ssmpvA  , ssmpvA_error );
  fMinuit->GetParameter(17 , ssgw   , ssgw_error  );
  fMinuit->GetParameter(18 , ssnorm , ssnorm_error);

  CSMinuit->GetBackground()->GetSemiBackground()->SetClwA( std::make_pair(sslwA  , sslwA_error ));
  CSMinuit->GetBackground()->GetSemiBackground()->SetClwB( std::make_pair(sslwB  , sslwB_error ));
  CSMinuit->GetBackground()->GetSemiBackground()->SetCmpvA(std::make_pair(ssmpvA  , ssmpvA_error ));
  CSMinuit->GetBackground()->GetSemiBackground()->SetCgwA( std::make_pair(ssgw   , ssgw_error  ));
  CSMinuit->GetBackground()->GetSemiBackground()->SetCnorm(std::make_pair(ssnorm , ssnorm_error));

  // CSMinuit = this;

  // //create Minuit
  // const int CPAR = 47;
  // const int NPAR = CPAR;
  // fMinuit = new TMinuit(NPAR);
  // fMinuit->SetFCN(fcnEstimateWhatIsLeft);

  // Double_t arglist[11];
  // Int_t ierflg = 0;

  // arglist[0] = 0.5;
  // fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // double slwA,slwB,smpvA,smpvB,smpvC,sgwA,sgwB,snorm;
  // double blwA, blwB, bmpvA, bmpvB, bgw, bnorm, bC;

  // slwA  = fTrueSignal->fClwA.first;  
  // slwB  = fTrueSignal->fClwB.first;	 
  // smpvA = fTrueSignal->fCmpvA.first;
  // smpvB = fTrueSignal->fCmpvB.first; 
  // smpvC = fTrueSignal->fCmpvC.first; 
  // sgwA  = fTrueSignal->fCgwA.first;	 
  // sgwB  = fTrueSignal->fCgwB.first;	 
  // snorm = fTrueSignal->fCnorm.first;

  // blwA  = fTrueBackground->fClwA.first;
  // blwB  = fTrueBackground->fClwA.first;
  // bmpvA = fTrueBackground->fCmpvA.first;
  // bmpvB = fTrueBackground->fCmpvB.first;
  // bgw   = fTrueBackground->fCgwA.first;
  // bnorm = fTrueBackground->fCnorm.first;
  // bC    = 0;//fTrueBackground->fCconst.first;

  // //define parameters and step
  // fMinuit->mnparm(0,  "slw A"   , slwA      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(1,  "slw B"   , slwB      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(2,  "smpv A"  , smpvA     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(3,  "smpv B"  , smpvB     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(4,  "smpv C"  , smpvC     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(5,  "sgw A"   , sgwA      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(6,  "sgw B"   , sgwB      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(7,  "snorm"   , snorm     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(8,  "blwA"    , blwA      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(9, "bmpvA"    , bmpvA     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(10, "bmpvB"   , bmpvB     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(11, "bgw"     , bgw       , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(12, "bnorm"   , bnorm     , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(13, "bconst"  , bC        , 0.1, 0   , 1   , ierflg);
  // fMinuit->mnparm(14, "ssmpv A" , 2.       , 0.1, 0   , 5 , ierflg);
  // fMinuit->mnparm(15, "ssgw"    , 0.1      , 0.1, 0   , 0   , ierflg);
  // fMinuit->mnparm(16, "ssnorm"  , 0.1    , 0.1, 0   , 1, ierflg);
  // for(int i = 0; i < 30; i++)
  //   fMinuit->mnparm(17+i, "e"  , 1    , 100, 0   , 0, ierflg);
    
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

  // // Now ready for minimization step
  // arglist[0] = 1;
  // fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  // arglist[0] = 500000;
  // arglist[1] = 1.;
  // fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // // Print results
  // Double_t amin,edm,errdef;
  // Int_t nvpar,nparx,icstat;
  // fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  // fMinuit->mnprin(3,amin);

  // double dummy1, dummy2;
  // std::vector<double> x,x_error,y,y_error;

  // for(int i = 0; i < 30; i++){
  //   x.push_back(fRR[i].first);
  //   x_error.push_back(fRR[i].second);
  //   fMinuit->GetParameter(17+i,dummy1,dummy2);
  //   y.push_back(dummy1);
  //   y_error.push_back(dummy2);
  // }

  // TGraphErrors* tg = new TGraphErrors(x.size(),&x[0],&y[0],&x_error[0],&y_error[0]);
  // tg->SetMarkerStyle(20);
  // tg->Draw("ap");
  // gPad->Update();gPad->WaitPrimitive();
}

//********************************************************************
void CoherentSample::fcnEstimateWhatIsLeft(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double lwA    = par[0];
  double lwB    = par[1];
  double mpvA   = par[2];
  double mpvB   = par[3];
  double mpvC   = par[4];
  double gwA    = par[5];
  double gwB    = par[6];
  double norm   = par[7];
  double blwA   = par[8];
  double bmpvA  = par[9];
  double bmpvB  = par[10];
  double bgw    = par[11];
  double bnorm  = par[12];
  double bC     = par[13];
  double sslwA  = par[14];
  double sslwB  = par[15];
  double ssmpv  = par[16];
  double ssgw   = par[17];
  double ssnorm = par[18];
  
  //define function to fit to each histogram
  double lw,mpv,gw,blw,bmpv,sslw;
  double Likelihood = 0;
  // TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangausPlusGausPlusConstant,1,30,12);
  TF1* flg = new TF1("flg",CoherentFitUtils::TripleLangausPlusConstant,1,30,13);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = CSMinuit->fRR[ibin].first;
    //mpv  = mpvA*(mpvB*x-1)/(mpvB*x+1)+mpvC;
    mpv  = mpvA*((mpvB*x-1)/(mpvB*x+1)+mpvC);
    gw   = gwA/(x+1)+gwB;
    lw   = (lwA/x-1)/x+lwB;

    bmpv = bmpvA+bmpvB*x;
    blw  = blwA;//  + blwB*x;

    sslw = sslwA/(x+1)+sslwB;
  
    double params[] = {lw,mpv,norm,gw,blw,bmpv,bnorm,bgw,sslw,ssmpv,ssnorm,ssgw,bC};
    flg->SetParameters(params);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  }
  f = -Likelihood;

  // //set each parameter
  // double lwA    = par[0];
  // double lwB    = par[1];
  // double mpvA   = par[2];
  // double mpvB   = par[3];
  // double mpvC   = par[4];
  // double gwA    = par[5];
  // double gwB    = par[6];
  // double norm   = par[7];
  // double blwA   = par[8];
  // double bmpvA  = par[9];
  // double bmpvB  = par[10];
  // double bgw    = par[11];
  // double bnorm  = par[12];
  // double bC     = par[13];
  // double ssmpv  = par[14];
  // double ssgw   = par[15];
  // double ssnorm = par[16];
  // double sslwv[30];
  // for(int i = 0; i < 30; i++)sslwv[i]=par[17+1];
   
  // //define function to fit to each histogram
  // double lw,mpv,gw,blw,bmpv,sslw;
  // double Likelihood = 0;
  // // TF1* flg = new TF1("flg",CoherentFitUtils::DoubleLangausPlusGausPlusConstant,1,30,12);
  // TF1* flg = new TF1("flg",CoherentFitUtils::TripleLangausPlusConstant,1,30,13);
  // for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
  //   //get RR depending parameters
  //   double x   = CSMinuit->fRR[ibin].first;
  //   //mpv  = mpvA*(mpvB*x-1)/(mpvB*x+1)+mpvC;
  //   mpv  = mpvA*((mpvB*x-1)/(mpvB*x+1)+mpvC);
  //   gw   = gwA/(x+1)+gwB;
  //   lw   = (lwA/x-1)/x+lwB;

  //   bmpv = bmpvA+bmpvB*x;
  //   blw  = blwA;//  + blwB*x;

  //   sslw = sslwv[ibin];//sslwA/(x+1)+sslwB;
    
  //   double params[] = {lw,mpv,norm,gw,blw,bmpv,bnorm,bgw,sslw,ssmpv,ssnorm,ssgw,bC};
  //   flg->SetParameters(params);
  //   Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);
  // }
  // f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitOnlySignal(){
//********************************************************************

  CSMinuit = this;

  //create Minuit
  const int CPAR = 8;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnOnlySignal);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set starting values and step sizes for parameters
  double lwA,lwB,lwC,mpvA,mpvB,mpvC,mpvD,gwA,gwB,gwC,norm;
  double lwA_error,lwB_error,lwC_error,mpvA_error,mpvB_error,mpvC_error,mpvD_error,gwA_error,gwB_error,gwC_error,norm_error;
    
  mpvA = fCmpvA.first;
  mpvB = fCmpvB.first;
  mpvC = fCmpvC.first;

  gwA = fCgwA.first;
  gwB = fCgwB.first;
  
  lwA = fClwA.first;
  lwB = fClwB.first;

  norm = fCnorm.first;

  fMinuit->mnparm(0,  "smpv A", mpvA   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(1,  "smpv B", mpvB   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(2,  "smpv C", mpvC   , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(3,  "lwA"   , lwA    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(4,  "lwB"   , lwB    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(5,  "gwA"   , gwA    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(6,  "gwB"   , gwB    , 0.001, 0, 0, ierflg);
  fMinuit->mnparm(7,  "snorm" , norm   , 0.001, 0, 1, ierflg);

  arglist[0] = 500000;
  arglist[1] = 1;
  fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  //retrieve parameters for next stage
  fMinuit->GetParameter(0,  fCmpvA.first, fCmpvA.second);
  fMinuit->GetParameter(1,  fCmpvB.first, fCmpvB.second);
  fMinuit->GetParameter(2,  fCmpvC.first, fCmpvC.second);
  fMinuit->GetParameter(3, fClwA.first,  fClwA.second);
  fMinuit->GetParameter(4, fClwB.first,  fClwB.second);
  fMinuit->GetParameter(5, fCgwA.first,  fCgwA.second);
  fMinuit->GetParameter(6, fCgwB.first,  fCgwB.second);
  fMinuit->GetParameter(7, fCnorm.first, fCnorm.second);
>>>>>>> develop

  StoreCoherentFits();
}

//********************************************************************
void CoherentSample::fcnOnlySignal(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//********************************************************************

  //set each parameter
  double mpvA  = par[0];
  double mpvB  = par[1];
  double mpvC  = par[2];
  double lwA   = par[3];
  double lwB   = par[4];
  double gwA   = par[5];
  double gwB   = par[6];
  double norm  = par[7];

  //define function to fit to each histogram
  double mpv,fw,gw,lw;
  double Likelihood = 0;
  TF1* flg = new TF1("flg",CoherentFitUtils::Langaus,1,30,4);
  for(int ibin = 0; ibin < (int)CSMinuit->fh.size(); ibin++){
    //get RR depending parameters
    double x   = CSMinuit->fRR[ibin].first;
    mpv  = mpvA*(mpvB*x-1)/(mpvB*x+1)+mpvC;
    gw   = gwA/(x+1)+gwB;
    lw   = (lwA/x-1)/x+lwB;

    flg->SetParameters(lw,mpv,norm,gw);
    Likelihood = Likelihood + CoherentFitUtils::ComputeLikelihood(CSMinuit->fh[ibin],flg,CSMinuit->fIntegral[ibin]);

    // CSMinuit->fh[ibin]->Draw();
    // flg->Draw("same");
    // gPad->Update();
    // gPad->WaitPrimitive();
      
  }
  f = -Likelihood;
}

//********************************************************************
void CoherentSample::CoherentFitSignalPlusBackgroundToy(bool &result){
//********************************************************************

<<<<<<< HEAD
=======
  // debug = true;
  
>>>>>>> develop
  CSMinuit = this;

  //create Minuit
  const int CPAR = 19;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
<<<<<<< HEAD
  fMinuit->SetFCN(fcnSignalPlusBackgroundCheb);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
=======
  fMinuit->SetFCN(fcnSignalPlusBackground);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
>>>>>>> develop
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  double slwA,slwB,smpvA,smpvB,smpvC,sgwA,sgwB,snorm;
  double slwA_error,slwB_error,smpvA_error,smpvB_error,smpvC_error,sgwA_error,sgwB_error,snorm_error;
  double blwA, blwB, bmpvA, bmpvB, bgw, bnorm, bC;
  double blwA_error, blwB_error, bmpvA_error, bmpvB_error, bnorm_error, bgw_error, bC_error;
  double sslwA,sslwB,ssmpvA,ssmpvB,ssmpvC,ssgw,ssnorm;
  double sslwA_error,sslwB_error,ssmpvA_error,ssmpvB_error,ssmpvC_error,ssgw_error,ssnorm_error;

<<<<<<< HEAD
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

=======
  slwA  = fSignal->fClwA.first;  
  slwB  = fSignal->fClwB.first;	 
  smpvA = fSignal->fCmpvA.first;
  smpvB = fSignal->fCmpvB.first; 
  smpvC = fSignal->fCmpvC.first; 
  sgwA  = fSignal->fCgwA.first;	 
  sgwB  = fSignal->fCgwB.first;	 
  snorm = fSignal->fCnorm.first;
  
  blwA  = fBackground->fClwA.first;
  bmpvA = fBackground->fCmpvA.first;
  bmpvB = fBackground->fCmpvB.first;
  bgw   = fBackground->fCgwA.first;
  bnorm = fBackground->fCnorm.first;
  bC    = 0;//fBackground->fCconst.first;
  	  
  sslwA  = fBackground->fSemiBackground->fClwA.first;  
  sslwB  = fBackground->fSemiBackground->fClwB.first;	 
  ssmpvA = fBackground->fSemiBackground->fCmpvA.first;
  ssgw   = fBackground->fSemiBackground->fCgwA.first;	 
  ssnorm = fBackground->fSemiBackground->fCnorm.first;

  double conv    = 0.001;
  fMinuit->mncler();
  fMinuit->mnparm(0,  "slw A"   , slwA     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(1,  "slw A"   , slwB     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(2,  "smpv A"  , smpvA    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(3,  "smpv B"  , smpvB    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(4,  "smpv C"  , smpvC    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(5,  "sgw A"   , sgwA     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(6,  "sgw B"   , sgwB     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(7,  "snorm"   , snorm    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(8,  "blwA"    , blwA     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(9,  "bmpvA"   , bmpvA    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(10, "bmpvB"   , bmpvB    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(11, "bgw"     , bgw      , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(12, "bnorm"   , bnorm    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(13, "bconst"  , bC       , 0.1, 0   , 1   , ierflg);
  fMinuit->mnparm(14, "sslw A"  , sslwA    , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(15, "sslw B"  , sslwB    , 0.1, 0   , 2   , ierflg);
  fMinuit->mnparm(16, "ssmpv A" , ssmpvA   , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(17, "ssgw"    , ssgw     , 0.1, 0   , 0   , ierflg);
  fMinuit->mnparm(18, "ssnorm"  , ssnorm   , 0.1, 0   , 0   , ierflg);

  fMinuit->FixParameter(13);
  fMinuit->FixParameter(18);
  
>>>>>>> develop
  // Now ready for minimization step
  arglist[0] = 1;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 500000;
  arglist[1] = 1.;
  fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

<<<<<<< HEAD
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
=======
  if(TMath::IsNaN(edm))
    result = false;
  if(edm < conv)
    result = true;
  
  //retrieve parameters for next stage
  if(result){
    fMinuit->GetParameter(2 , smpvA, smpvA_error);
    fMinuit->GetParameter(3 , smpvB, smpvB_error);
    fMinuit->GetParameter(4 , smpvC, smpvC_error);
    
    fSignal->SetCmpvA(std::make_pair(smpvA,smpvA_error));
    fSignal->SetCmpvB(std::make_pair(smpvB,smpvB_error));
    fSignal->SetCmpvC(std::make_pair(smpvC,smpvC_error));

    // fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
    // // fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
    // arglist[0] = 500000;
    // arglist[1] = 3;
    // arglist[2] = 4;
    // arglist[3] = 5;
    // fMinuit->mnexcm("MINOS", arglist, 4, ierflg);
  }
  
  if(!result)TrySoloPeaks(result);
}

//********************************************************************
void CoherentSample::TrySoloPeaks(bool &result){
//********************************************************************

  std::cout << "-------------------------------------" << std::endl;
  std::cout << "STANDARD FIT DIT NOT WORK" << std::endl;
  std::cout << "TRYING SOLO PEAKS" << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  
  CSMinuit = this;

  //create Minuit
  const int CPAR = 19;
  const int NPAR = CPAR;
  fMinuit = new TMinuit(NPAR);
  fMinuit->SetFCN(fcnSignalPlusBackground);

  Double_t arglist[11];
  Int_t ierflg = 0;

  arglist[0] = 0.5;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  double slwA,slwB,smpvA,smpvB,smpvC,sgwA,sgwB,snorm;
  double slwA_error,slwB_error,smpvA_error,smpvB_error,smpvC_error,sgwA_error,sgwB_error,snorm_error;
  double blwA, blwB, bmpvA, bmpvB, bgw, bnorm, bC;
  double blwA_error, blwB_error, bmpvA_error, bmpvB_error, bnorm_error, bgw_error, bC_error;
  double sslwA,sslwB,ssmpvA,ssmpvB,ssmpvC,ssgw,ssnorm;
  double sslwA_error,sslwB_error,ssmpvA_error,ssmpvB_error,ssmpvC_error,ssgw_error,ssnorm_error;

  slwA  = fSignal->fClwA.first;  
  slwB  = fSignal->fClwB.first;	 
  smpvA = fSignal->fCmpvA.first;
  smpvB = fSignal->fCmpvB.first; 
  smpvC = fSignal->fCmpvC.first; 
  sgwA  = fSignal->fCgwA.first;	 
  sgwB  = fSignal->fCgwB.first;	 
  snorm = fSignal->fCnorm.first;
  
  blwA  = fBackground->fClwA.first;
  bmpvA = fBackground->fCmpvA.first;
  bmpvB = fBackground->fCmpvB.first;
  bgw   = fBackground->fCgwA.first;
  bnorm = fBackground->fCnorm.first;
  bC    = 0;//fBackground->fCconst.first;
  	  
  sslwA  = fBackground->fSemiBackground->fClwA.first;  
  sslwB  = fBackground->fSemiBackground->fClwB.first;	 
  ssmpvA = fBackground->fSemiBackground->fCmpvA.first;
  ssgw   = fBackground->fSemiBackground->fCgwA.first;	 
  ssnorm = fBackground->fSemiBackground->fCnorm.first;

  fMinuit->mncler();
  fMinuit->mnparm(0,  "slw A"   , slwA     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(1,  "slw A"   , slwB     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(2,  "smpv A"  , smpvA    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(3,  "smpv B"  , smpvB    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(4,  "smpv C"  , smpvC    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(5,  "sgw A"   , sgwA     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(6,  "sgw B"   , sgwB     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(7,  "snorm"   , snorm    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(8,  "blwA"    , blwA     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(9,  "bmpvA"   , bmpvA    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(10, "bmpvB"   , bmpvB    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(11, "bgw"     , bgw      , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(12, "bnorm"   , bnorm    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(13, "bconst"  , bC       , 1, 0   , 1   , ierflg);
  fMinuit->mnparm(14, "sslw A"  , sslwA    , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(15, "sslw B"  , sslwB    , 1, 0   , 2   , ierflg);
  fMinuit->mnparm(16, "ssmpv A" , ssmpvA   , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(17, "ssgw"    , ssgw     , 1, 0   , 0   , ierflg);
  fMinuit->mnparm(18, "ssnorm"  , ssnorm   , 1, 0   , 0   , ierflg);

  fMinuit->FixParameter(0);
  fMinuit->FixParameter(1);
  fMinuit->FixParameter(5);
  fMinuit->FixParameter(6);
  fMinuit->FixParameter(7);
  fMinuit->FixParameter(8);
  fMinuit->FixParameter(11);
  fMinuit->FixParameter(12);
  fMinuit->FixParameter(13);
  fMinuit->FixParameter(14);
  fMinuit->FixParameter(15);
  fMinuit->FixParameter(17);
  fMinuit->FixParameter(18);
    
  // Now ready for minimization step
  arglist[0] = 0;
  fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  arglist[0] = 500000;
  arglist[1] = 100.;
  double conv = 0.001*100;
  fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  fMinuit->mnprin(3,amin);

  if(TMath::IsNaN(edm))result = false;
  else if(edm < conv)  result = true;

  if(!result){
    std::cout << "SOLO PEAKS DID NOT WORK, SKIPPING TOY" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    return;
  }

  std::cout << "SOLO PEAKS WORKED! GOING FOR FULL FIT" << std::endl;
  std::cout << "-------------------------------------" << std::endl;

  result = false;
  
  //retrieve parameters for next stage
  fMinuit->GetParameter(2 , smpvA, smpvA_error);
  fMinuit->GetParameter(3 , smpvB, smpvB_error);
  fMinuit->GetParameter(4 , smpvC, smpvC_error);
  fMinuit->GetParameter(9 , bmpvA, bmpvA_error);
  fMinuit->GetParameter(10, bmpvB, bmpvB_error);
  fMinuit->GetParameter(16, ssmpvA, ssmpvA_error);

  double precision[4] = {1,10,100,1000};
  for(int i = 0; i < 4; i++){
    fMinuit->mncler();
    fMinuit->mnfree(0);
    fMinuit->mnparm(0,  "slw A"   , slwA     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(1,  "slw A"   , slwB     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(2,  "smpv A"  , smpvA    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(3,  "smpv B"  , smpvB    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(4,  "smpv C"  , smpvC    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(5,  "sgw A"   , sgwA     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(6,  "sgw B"   , sgwB     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(7,  "snorm"   , snorm    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(8,  "blwA"    , blwA     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(9,  "bmpvA"   , bmpvA    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(10, "bmpvB"   , bmpvB    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(11, "bgw"     , bgw      , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(12, "bnorm"   , bnorm    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(13, "bconst"  , bC       , 0.1, 0   , 1   , ierflg);
    fMinuit->mnparm(14, "sslw A"  , sslwA    , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(15, "sslw B"  , sslwB    , 0.1, 0   , 2   , ierflg);
    fMinuit->mnparm(16, "ssmpv A" , ssmpvA   , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(17, "ssgw"    , ssgw     , 0.1, 0   , 0   , ierflg);
    fMinuit->mnparm(18, "ssnorm"  , ssnorm   , 0.1, 0   , 0   , ierflg);
    
    fMinuit->FixParameter(13);
    fMinuit->FixParameter(18);
    
    // Now ready for minimization step
    arglist[0] = 1;
    fMinuit->mnexcm("SET STR", arglist, 1, ierflg);
    arglist[0] = 500000;
    arglist[1] = precision[i];
    conv       = 0.001*precision[i];
    fMinuit->mnexcm("HESSE", arglist, 1, ierflg);
    fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    
    // Print results
    fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    fMinuit->mnprin(3,amin);
    
    if(TMath::IsNaN(edm))result = false;
    else if(edm < conv)result = true;
    if(!result){
      std::cout << "SOLO PEAKS FULL FIT DID NOT WORK, DECREASING PRECISION" << std::endl;
      std::cout << "-------------------------------------" << std::endl;
      continue;
    }
    else{
      std::cout << "SOLO PEAKS FULL FIT DID WORK! RETRIEVING PARAMETERS" << std::endl;
      //retrieve parameters for next stage
      fMinuit->GetParameter(2 , smpvA, smpvA_error);
      fMinuit->GetParameter(3 , smpvB, smpvB_error);
      fMinuit->GetParameter(4 , smpvC, smpvC_error);
      
      fSignal->SetCmpvA(std::make_pair(smpvA,smpvA_error));
      fSignal->SetCmpvB(std::make_pair(smpvB,smpvB_error));
      fSignal->SetCmpvC(std::make_pair(smpvC,smpvC_error));
      break;
    }
>>>>>>> develop
  }
}
