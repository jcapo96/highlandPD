#include "dEdxUtils.hxx"

//********************************************
Double_t dEdxUtils::Langaus(Double_t *x, Double_t *par) {
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
  Double_t np = 100.0;      // number of convolution steps
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
void dEdxUtils::SetFitParameters(TF1* f, Double_t max) {
//********************************************

  f->SetParameters(1,60,5000,5);
  f->SetParLimits(0,1,5);
  f->SetParLimits(3,1,5);
  f->SetParameter(1,max);
  f->SetParLimits(1,max-5,max+10);
  f->SetRange(max-7.5,max+10);
}

//********************************************************************
Double_t dEdxUtils::Recombination(Double_t E){
//********************************************************************

  Double_t LAr_density=1.39;
  Double_t alp=0.93;
  Double_t bet=0.212;
  Double_t dedx=1.9;
  Double_t xsi=bet*dedx/(LAr_density*E);
  Double_t xsi0=bet*dedx/(LAr_density*0.4867);
  Double_t rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}

//********************************************************************
Double_t dEdxUtils::ElectricField(const SpaceCharge* sce, Double_t x, Double_t y, Double_t z){
//********************************************************************
  
  Double_t E0 = 0.4867;
  Double_t ex = E0+E0*sce->GetEFieldXAtPoint(x,y,z);
  Double_t ey = 0.0+E0*sce->GetEFieldYAtPoint(x,y,z);
  Double_t ez = 0.0+E0*sce->GetEFieldZAtPoint(x,y,z);
  return sqrt(ex*ex+ey*ey+ez*ez);
}

//********************************************************************
bool dEdxUtils::IsInterestingHit(AnaHitPD& hit){
//********************************************************************

  bool ItIs = false;
  
  if(hit.Position.X() > XMIN && hit.Position.X() < XMAX &&
     hit.Position.Y() > YMIN && hit.Position.Y() < YMAX &&
     hit.Position.Z() > ZMIN && hit.Position.Z() < ZMAX)
    ItIs = true;
  
  return ItIs;
}
