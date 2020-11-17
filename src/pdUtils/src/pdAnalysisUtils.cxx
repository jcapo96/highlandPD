#include "pdAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include <TH3F.h>
#include <TH2F.h>


bool debug = false;

//data for range-momentum conversion, muons
//http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf divided by LAr density for cm

float Range_grampercm_new[29] = {
  9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396,
  1.063E2/1.396,  1.725E2/1.396, 2.385E2/1.396, 4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396,
  2.297E3/1.396,  4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396, 1.013E4/1.396, 1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396,
  4.326E4/1.396,  5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396};

float KE_MeV_new[29]= {
  10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
  400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
  20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

//data for range-momentum conversion, protons

double Range_gpercm_P_Nist[31]={
  1.887E-1/1.396, 3.823E-1/1.396, 6.335E-1/1.396, 1.296/1.396,   2.159/1.396,   7.375/1.396,   1.092E1/1.396, 2.215E1/1.396,
  3.627E1/1.396,  5.282E1/1.396,  7.144E1/1.396,  9.184E1/1.396, 1.138E2/1.396, 1.370E2/1.396, 1.614E2/1.396, 1.869E2/1.396, 
  2.132E2/1.396,  2.403E2/1.396,  2.681E2/1.396,  2.965E2/1.396, 3.254E2/1.396, 3.548E2/1.396, 3.846E2/1.396, 4.148E2/1.396, 
  4.454E2/1.396,  7.626E2/1.396,  1.090E3/1.396,  1.418E3/1.396, 1.745E3/1.396, 2.391E3/1.396, 3.022E3/1.396};

double KE_MeV_P_Nist[31]={
  10,   15,   20,   30,   40,   80,   100, 150,
  200,  250,  300,  350,  400,  450,  500, 550, 
  600,  650,  700,  750,  800,  850,  900, 950, 
  1000, 1500, 2000, 2500, 3000, 4000, 5000};

TGraph const KEvsR(29, Range_grampercm_new, KE_MeV_new);
TSpline3 const KEvsR_spline3("KEvsRS", &KEvsR);

TGraph const RvsKE(29, KE_MeV_new, Range_grampercm_new);
TSpline3 const RvsKE_spline3("RvsKES", &RvsKE);

TGraph const RvsKE_P(31, KE_MeV_P_Nist, Range_gpercm_P_Nist);
TSpline3 const RvsKE_P_spline3("RvsKE_P_S", &RvsKE_P);


TFile* dEdX_template_file = new TFile( (std::string(getenv("PIONANALYSISROOT"))+"/data/dEdxrestemplates.root").c_str(), "OPEN" );
std::map< int, TProfile* > templates;

TProfile* ProtonTemplate = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
/*
templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
*/


TFile* E_field_file = new TFile( (std::string(getenv("PIONANALYSISROOT"))+"/data/SCE_DataDriven_180kV_v3.root").c_str(), "OPEN" );

TH3F* ex_neg = (TH3F*)E_field_file->Get("Reco_ElecField_X_Neg");
TH3F* ey_neg = (TH3F*)E_field_file->Get("Reco_ElecField_Y_Neg");
TH3F* ez_neg = (TH3F*)E_field_file->Get("Reco_ElecField_Z_Neg");
TH3F* ex_pos = (TH3F*)E_field_file->Get("Reco_ElecField_X_Pos");
TH3F* ey_pos = (TH3F*)E_field_file->Get("Reco_ElecField_Y_Pos");
TH3F* ez_pos = (TH3F*)E_field_file->Get("Reco_ElecField_Z_Pos");


std::string X_correction_name = std::string(getenv("PIONANALYSISROOT"))+"/data/Xcalo_r5387.root";
std::string YZ_correction_name = std::string(getenv("PIONANALYSISROOT"))+"/data/YZcalo_r5387.root";


TFile* X_correction_file  = new TFile( X_correction_name.c_str(), "OPEN" );
TFile* YZ_correction_file = new TFile( YZ_correction_name.c_str(), "OPEN" );

UInt_t planeID=2;
std::string hist_name = "dqdx_X_correction_hist_" + std::to_string(planeID);
TH1F* X_correction_hist = (TH1F*)X_correction_file->Get( hist_name.c_str() );

TH2F* YZ_neg = (TH2F*)YZ_correction_file->Get("correction_dqdx_ZvsY_negativeX_hist_2");
TH2F* YZ_pos = (TH2F*)YZ_correction_file->Get("correction_dqdx_ZvsY_positiveX_hist_2");


//*****************************************************************************
Float_t pdAnaUtils::ComputeRangeMomentum(double trkrange, int pdg){
//*****************************************************************************
  
  /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
     website:
     http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
     
     CSDA table values:
     float Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0,
     6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
     2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
     4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
     4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}; float KE_MeV[30] = {10, 14,
     20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000,
     4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000,
     200000, 300000, 400000};
     
     Functions below are obtained by fitting polynomial fits to KE_MeV vs
     Range (cm) graph. A better fit was obtained by splitting the graph into
     two: Below range<=200cm,a polynomial of power 4 was a good fit; above
     200cm, a polynomial of power 6 was a good fit
     
     Fit errors for future purposes:
     Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2
     err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7; Above 200cm,
     Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3
     err=5.9155E-9; p4 err=9.71709E-13; p5 err=7.22381E-17;p6
     err=1.9709E-21;*/
  
  //*********For muon, the calculations are valid up to 1.91E4 cm range
  //corresponding to a Muon KE of 40 GeV**********//
  
  /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
    website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
    
    CSDA values:
    double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
    1500, 2000, 2500, 3000, 4000, 5000};
    
    double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
    2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
    9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
    2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
    7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};
    
    Functions below are obtained by fitting power and polynomial fits to
    KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
    graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
    polynomial of power 6 was a good fit
    
    Fit errors for future purposes:
    For power function fit: a=0.388873; and b=0.00347075
    Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
    err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
    err=2.96958E-17;*/
  
  //*********For proton, the calculations are valid up to 3.022E3 cm range
  //corresponding to a Muon KE of 5 GeV**********//

  if (trkrange < 0 || std::isnan(trkrange)) {
    //    std::cout << "TrackMomentumCalculator   " 
    //              << "Invalid track range " << trkrange << " return -1" << std::endl;
    return -1.;
  }

    
  double KE, Momentum, M;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;
  
  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = KEvsR_spline3.Eval(trkrange);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    if (trkrange > 0 && trkrange <= 80)
      KE = 29.9317 * std::pow(trkrange, 0.586304);
    else if (trkrange > 80 && trkrange <= 3.022E3)
      KE =
        149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
        (4.34587E-6 * trkrange * trkrange * trkrange) +
        (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
        (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
        (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *
         trkrange);
    else
      KE = -999;
  } else
    KE = -999;
  
  if (KE < 0)
    Momentum = -999;
  else
    Momentum = std::sqrt((KE * KE) + (2 * M * KE));
  
  Momentum = Momentum / 1000;
  
  return Momentum;
}

//*****************************************************************************
Float_t pdAnaUtils::ComputeCSDARange(double beammom, int pdg){
//*****************************************************************************
    
  if (beammom < 0 || std::isnan(beammom)) {
    //    std::cout << "CSDARangeCalculator   " 
    //              << "Invalid beam mom " << beammom << " return -1" << std::endl;
    return -1.;
  }

  double KE, M, CSDARange;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;
  
  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_spline3.Eval(KE);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_P_spline3.Eval(KE);
  }
  else 
    CSDARange=-1;
  
  return CSDARange;
}


//********************************************************************
Float_t pdAnaUtils::ComputePIDA(const AnaParticlePD &track) {
//********************************************************************

  Float_t cut=30;

  Float_t PIDA=0;
  Int_t ncontrib=0;
  for (Int_t i=0;i<3;i++){
    for (UInt_t j=0;j<track.Hits[i].size();j++){
      if (track.Hits[i][j].ResidualRange<cut && track.Hits[i][j].ResidualRange>0){
        ncontrib++;
        PIDA += track.Hits[i][j].dEdx*pow(track.Hits[i][j].ResidualRange,0.42);
      }
    }
  }
  if (ncontrib>0) PIDA /= ncontrib*1.;

  return PIDA;
}

//********************************************************************
Float_t pdAnaUtils::ComputeKineticEnergy(const AnaParticlePD &part) {
//********************************************************************


  Int_t plane=2;

  int nhits=part.Hits[plane].size();
  double kinetic=0;
  double res=0;
  for (int j=0;j<nhits;j++){
    double dedxi = part.Hits[plane][j].dEdx_corr;
//    double dedxi = part.dEdx[plane][i];
    double Residualrangei = part.Hits[plane][j].ResidualRange;
    kinetic = kinetic + dedxi * (Residualrangei - res);
    res = Residualrangei;
  }

  // convert to GeV
  return kinetic/units::GeV;
}

//********************************************************************
Float_t pdAnaUtils::ComputeDeDxFromDqDx(Float_t dqdx_adc, Int_t plane, Float_t x, Float_t y, Float_t z) {
//********************************************************************

  // Formula found for example at https://arxiv.org/abs/1306.1712v1
  
  
  // Paremeters from /cvmfs/dune.opensciencegrid.org/products/dune/protoduneana/v09_01_00/job/ProtoDUNECalibration.fcl
  
  /***************modified box model parameters and function*****************/
  double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
  double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
  double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
  double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm. In MeV/e
  //  double Efield1=0.50;//kV/cm protoDUNE electric filed

  // TODO: hit position not available
  double Efield1 = ComputeTotalEField(0,0,0);//kV/cm protoDUNE electric filed
  //  std::cout << x << " " << y << " " << z << " " << Efield1 << std::endl;
  //  double calib_factor =6.155e-3; //right cali constant for the run 5387. This converts from ADC to e
  //  double calib_factor[3] = {4.81e-3, 4.81e-3, 4.86e-3}; //right cali constant for the run 5387. This converts from ADC to e
  double calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //
  //  double norm_factor[3] = {1.0078, 1.0082, 0.9947};
    double norm_factor[3] = {1.0078, 1.0082, 0.9946};

   
  // dq/dx should be in e/cm
  // dE/dx is returned in MeV/cm

  double dqdx = dqdx_adc/calib_factor[plane]*norm_factor[plane];
  double beta = betap/(Rho*Efield1); // cm/MeV

  if (debug){
    std::cout << "  --> dqdx_adc, dqdx_e (pdAnaUtils)= " << dqdx_adc << ", " << dqdx << std::endl;
    
    std::cout << beta << " " << Wion << " " << alpha << std::endl;
    std::cout << "  --> dedx (pdAnaUtils)= " << (exp(dqdx*beta*Wion)-alpha)/beta << std::endl;
  }
  return (exp(dqdx*beta*Wion)-alpha)/beta;
}

//********************************************************************
Float_t pdAnaUtils::ComputeDqDxFromDeDx(Float_t dedx, Int_t plane) {
//********************************************************************

  // Paremeters from /cvmfs/dune.opensciencegrid.org/products/dune/protoduneana/v09_01_00/job/ProtoDUNECalibration.fcl
  
  /***************modified box model parameters and function*****************/
  double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
  double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
  double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
  double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm. In MeV/e
  //  double Efield1=0.50;//kV/cm protoDUNE electric filed

  // TODO: hit position not available
  double Efield1 = ComputeTotalEField(-50,450,40);//kV/cm protoDUNE electric filed

  if (debug)
    std::cout << Efield1 << std::endl;
  //  double calib_factor =6.155e-3; //right cali constant for the run 5387. This converts from ADC to e
  double calib_factor[3] = {4.81e-3, 4.81e-3, 4.86e-3}; //right cali constant for the run 5387. This converts from ADC to e
  double norm_factor[3] = {1.0078, 1.0082, 0.9947};

   
  // dq/dx should be in e/cm
  // dE/dx is returned in MeV/cm

  double beta = betap/(Rho*Efield1); // cm/MeV
  
  return  log(dedx*beta + alpha)/(beta*Wion)*calib_factor[plane]/norm_factor[plane];
}

//********************************************************************
Float_t pdAnaUtils::ComputeTotalEField( Float_t x, Float_t y, Float_t z ){
//********************************************************************

  if( x >= 0 ){
    Float_t ex = 0.5 + 0.5 * ex_pos->GetBinContent( ex_pos->FindBin( x, y, z ) );
    Float_t ey =       0.5 * ey_pos->GetBinContent( ey_pos->FindBin( x, y, z ) );
    Float_t ez =       0.5 * ez_pos->GetBinContent( ez_pos->FindBin( x, y, z ) );
    return sqrt( (ex*ex) + (ey*ey) + (ez*ez) );
  }
  else if( x < 0 ){
    Float_t ex = 0.5 + 0.5 * ex_neg->GetBinContent( ex_neg->FindBin( x, y, z ) );
    Float_t ey =       0.5 * ey_neg->GetBinContent( ey_neg->FindBin( x, y, z ) );
    Float_t ez =       0.5 * ez_neg->GetBinContent( ez_neg->FindBin( x, y, z ) );
    return sqrt( (ex*ex) + (ey*ey) + (ez*ez) );
  }
  else return 0.5;
}

//********************************************************************
Float_t* pdAnaUtils::ExtrapolateToZ(const AnaParticlePD* part, Float_t z, Float_t* posz) {
//********************************************************************

  // Get track starting point
  double xi,yi,zi;
  xi = part->PositionStart[0];
  yi = part->PositionStart[1];
  zi = part->PositionStart[2];

  // Get track direction
  double ux, uy, uz;
  ux = part->DirectionStart[0];
  uy = part->DirectionStart[1];
  uz = part->DirectionStart[2];

  // zi = z + n * uz; get n
  double n = (zi-z)/uz;

  // Get position at z plane
  posz[0] = xi-n*ux;
  posz[1] = yi-n*uy;
  posz[2] = zi-n*uz;
  
  //std::cout << posz[0] << " " << posz[1] << " " << posz[2] << std::endl;

  return posz;
}



//********************************************************************
void pdAnaUtils::ComputeBinnedDeDx(const AnaParticlePD* part, Float_t max_resrange, Int_t nbins, Float_t** avg_dedx){  
//********************************************************************

  for (Int_t i=0;i<3;i++){
    for (Int_t j=0;j<nbins;j++){
      avg_dedx[i][j]=0;
    }
  }
    
  
  Float_t bin_width = max_resrange/nbins;
  
  for (Int_t i=0;i<3;i++){
    for (Int_t k=0;k<nbins;k++){
      Float_t cut_min = k*bin_width;
      Float_t cut_max = (k+1)*bin_width;      
      Float_t ncontrib=0;
      //      std::cout << "k = " << k << std::endl;
      for (Int_t j=0;j<std::min((Int_t)part->NHitsPerPlane[i],(Int_t)NMAXHITSPERPLANE);j++){
        // a protection against crazy values
        if (part->Hits[i][j].ResidualRange<0.01 || part->Hits[i][j].dEdx<0.01 || part->Hits[i][j].dEdx>100) continue;
        if (part->Hits[i][j].ResidualRange<cut_max && part->Hits[i][j].ResidualRange>cut_min){
          ncontrib++;
          avg_dedx[i][k] +=part->Hits[i][j].dEdx;
        }
      }      
      if (ncontrib>0)
        avg_dedx[i][k] /=ncontrib;

    }
  }
  
}

//********************************************************************
AnaTrueParticle* pdAnaUtils::FindBeamTrueParticle(const AnaSpillB& spill){  
//********************************************************************


  AnaTrueParticle* beampart=NULL;
  
  AnaBeamPD* beam         = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;

  Float_t beammom=0;
  if (beamPart){
    if (beamPart->TrueObject){
      beammom = static_cast<AnaTrueParticleB*>(beamPart->TrueObject)->Momentum;
    }
  }
  
  if (spill.TrueParticles.size() > 0){
    for (UInt_t i =0; i< spill.TrueParticles.size();i++){
      if (beammom == spill.TrueParticles[i]->Momentum){
        beampart = static_cast<AnaTrueParticle*> (spill.TrueParticles[i]);
        break;
      }
    }
  }

  return beampart;

}


//********************************************************************
void pdAnaUtils::AddParticles(AnaParticlePD* part1, AnaParticlePD* part2){  
//********************************************************************

  part1->Length     += part2->Length;
  part1->NHits      += part2->NHits;
  
  anaUtils::CopyArray(part2->DirectionEnd,    part1->DirectionEnd,   3);
  anaUtils::CopyArray(part2->PositionEnd,     part1->PositionEnd,    4);
  
  AnaParticlePD* part1c = part1->Clone();
    
    
  for (Int_t i=0;i<3;i++){
    part1->NHitsPerPlane[i] += part2->NHitsPerPlane[i];

    Int_t last_hit=0;
    for (Int_t j=0;j<std::min((Int_t)NMAXHITSPERPLANE,part2->NHitsPerPlane[i]);j++){
      Int_t offset = 0;
      part1->Hits[i][j+offset].dEdx           = part2->Hits[i][j].dEdx;          
      part1->Hits[i][j+offset].dQdx           = part2->Hits[i][j].dQdx;          
      part1->Hits[i][j+offset].dEdx_corr      = part2->Hits[i][j].dEdx_corr;     
      part1->Hits[i][j+offset].dQdx_corr      = part2->Hits[i][j].dQdx_corr;     
      part1->Hits[i][j+offset].ResidualRange  = part2->Hits[i][j].ResidualRange; 
      part1->Hits[i][j+offset].Position.SetX(   part2->Hits[i][j].Position.X());
      last_hit=j;
    }

    for (Int_t j=0;j<std::min((Int_t)part1->NHitsPerPlane[i],(Int_t)NMAXHITSPERPLANE-part2->NHitsPerPlane[i]);j++){
      Int_t offset = part1->NHitsPerPlane[i];
      part1->Hits[i][j+offset].dEdx           = part1c->Hits[i][j].dEdx;          
      part1->Hits[i][j+offset].dQdx           = part1c->Hits[i][j].dQdx;          
      part1->Hits[i][j+offset].dEdx_corr      = part1c->Hits[i][j].dEdx_corr;     
      part1->Hits[i][j+offset].dQdx_corr      = part1c->Hits[i][j].dQdx_corr;     
      part1->Hits[i][j+offset].Position.SetX(   part1c->Hits[i][j].Position.X());   
      part1->Hits[i][j+offset].ResidualRange  = part1c->Hits[i][j].ResidualRange+part2->Hits[i][last_hit].ResidualRange;
    }    
  }

  delete part1c;
  
  for (int i=0; i<3; i++) {
    part1->PIDA[i]     = part2->PIDA[i];
    part1->ReconPDG[i] = part2->ReconPDG[i];

    for (int j=0; j<10; j++) {
      part1->PID[i][j]   = part2->PID[i][j];
      part1->CALO[i][j] += part2->CALO[i][j];
    }
  }

  part1->RangeMomentum[0] = pdAnaUtils::ComputeRangeMomentum(part1->Length, 13);
  part1->RangeMomentum[1] = pdAnaUtils::ComputeRangeMomentum(part1->Length, 2212);

  // A pointer to the original particle
  //    Original = &part;
  
  
  
  //  part1->TrueEff       = part->TrueEff;
  //  part1->TruePur       = part->TruePur;
  
  
  //  AveragedQdx   = part->AveragedQdx;
  //  AveragedEdx   = part->AveragedEdx;
  //  MomentumError = part->MomentumError;

  // TODO: just to mark it as a broken track
  part1->NDOF       = 8888;
  //  Chi2          = part->Chi2;
  //  FitPDG        = part->FitPDG;
  //  Bunch         = part->Bunch;

  
  // TODO
  //  Daughters.clear();
  //  for (UInt_t i=0;i<part.Daughters.size();i++){
  //    Daughters.push_back(part.Daughters[i]->Clone());
  //  }

}

//********************************************************************
void pdAnaUtils::ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance){
//********************************************************************

  distance.clear();
  if (!part) return;
  if (part->Daughters.empty()) return;

  for(int i = 0; i < (int)part->Daughters.size(); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);

    //compute distance between part and daughter. Assume dau can be reconstructed backwards
    Float_t dis1 = sqrt(pow(part->PositionEnd[0]-dau->PositionStart[0],2)+pow(part->PositionEnd[1]-dau->PositionStart[1],2)+pow(part->PositionEnd[2]-dau->PositionStart[2],2));
    Float_t dis2 = sqrt(pow(part->PositionEnd[0]-dau->PositionEnd[0],2)+pow(part->PositionEnd[1]-dau->PositionEnd[1],2)+pow(part->PositionEnd[2]-dau->PositionEnd[2],2));

    //return lowest value
    if(dis1<dis2)distance.push_back(dis1);
    else distance.push_back(dis2);
  }
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(AnaEventB* event, Int_t ID){
//*****************************************************************************

  // Get all reconstructed tracks in the event
  AnaTrueParticleB** trueParticles = event->TrueParticles;
  Int_t nTrueParts                 = event->nTrueParticles;

  for (Int_t i=0;i<nTrueParts;i++){
    if (trueParticles[i]->ID == ID){
      return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
    }
  }

  return NULL;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID){
 //*****************************************************************************
 
  // Get all reconstructed tracks in the event
  for (UInt_t i=0;i<trueParticles.size();i++){
    if (trueParticles[i]->ID == ID){
      return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
    }
  }

  return NULL;
}


//********************************************************************
void pdAnaUtils::FillBeamDaughterCounters(AnaEventB& event, PDCounters& counters){
//********************************************************************

  AnaBeamPD* beam  = static_cast<AnaBeamPD*>(event.Beam); 
  AnaParticleB* beamPart = beam->BeamParticle;
  if (!beamPart) return;
  AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(beamPart->TrueObject);  
  if (!trueBeamPart) return;
  
  counters.ntrue_beamdaughter_piplus=0;
  counters.ntrue_beamdaughter_piminus=0;
  counters.ntrue_beamdaughter_pi0=0;
  counters.ntrue_beamdaughter_proton=0;
  counters.ntrue_beamdaughter_neutron=0;
  counters.ntrue_beamdaughter_nucleus=0;

  for(UInt_t i = 0; i < trueBeamPart->Daughters.size(); i++){
    AnaTrueParticlePD* trueBeamDaughter = pdAnaUtils::GetTrueParticle(&event, trueBeamPart->Daughters[i]);
    if (!trueBeamDaughter) continue;
    if      (trueBeamDaughter->PDG==211)  counters.ntrue_beamdaughter_piplus++;
    else if (trueBeamDaughter->PDG==-211) counters.ntrue_beamdaughter_piminus++;
    else if (trueBeamDaughter->PDG==111)  counters.ntrue_beamdaughter_pi0++;
    else if (trueBeamDaughter->PDG==2212) counters.ntrue_beamdaughter_proton++;
    else if (trueBeamDaughter->PDG==2112) counters.ntrue_beamdaughter_neutron++;
    else if (trueBeamDaughter->PDG>2212)  counters.ntrue_beamdaughter_nucleus++;
  }
}


//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID(const AnaParticlePD& part, TProfile * profile ){
//*****************************************************************************	

  double pid_chi2 = 0.; 
  int npt = 0;

  Int_t plane=2;
  
  profile = ProtonTemplate;
  
  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);
  
  //Ignore first and last point
  for( Int_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx_corr > 1000. )
      continue;

    int bin = profile->FindBin( part.Hits[plane][i].ResidualRange );

    if( bin >= 1 && bin <= profile->GetNbinsX() ){
      
      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;        
      }
      
      
      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;        
      }

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx_corr * part.Hits[plane][i].dEdx_corr;      
      dedx_res *= part.Hits[plane][i].dEdx_corr; 
      
      
      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )      
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx_corr - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) ); 
            
      ++npt;      
    }	
  }
		
  if( npt == 0 )	
    return std::make_pair(9999., -1);
	  		
  return std::make_pair(pid_chi2, npt); 	
}

//*****************************************************************************
bool pdAnaUtils::isBeamLike(AnaParticlePD* part, AnaBeamPD* beam ){
//*****************************************************************************	

  if (!part) return false;
  
  //From Owen Goodwins studies
  Float_t mccuts[7]  ={-3.,  7., -8.,  7., 27.5, 32.5, 0.93};
  Float_t datacuts[7]={ 0., 10., -5., 10., 30.,  35.0, 0.93};

  //cast the beam particle
  AnaParticlePD* beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if (!beampart) return false;
  
  Float_t beampos[3],beamdir[3], dist[3], dcos=0, cuts[7];      

  // different way of obtaining the beam position and angle for DATA and MC
  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  // For MC
  if (trueBeamPart){
    for (int i=0;i<3;i++){
      beampos[i] = trueBeamPart->Position[i]-trueBeamPart->Position[2]*(trueBeamPart->Direction[i]/trueBeamPart->Direction[2]);
      beamdir[i] = trueBeamPart->Direction[i];
    }
    for (int i=0;i<7;i++) cuts[i] = mccuts[i];
  }
  else{
    // For Data
    if(beam->nMomenta != 1 || beam->nTracks != 1)return false;
    for (int i=0;i<3;i++){
      beampos[i] = beampart->PositionEnd[i];
      beamdir[i] = beampart->DirectionEnd[i];
    }
    for (int i=0;i<7;i++) cuts[i] = datacuts[i];
  }
  
  // compute the difference in position and cos(angle) considering that particle could have been reconstructed backwards
  if(part->PositionEnd[2]<part->PositionStart[2] && part->PositionEnd[2]!=-999){
    for (int i=0;i<3;i++){
      dist[i] = part->PositionEnd[i] - beampos[i];
      dcos   += -1*(part->DirectionEnd[i])*beamdir[i];
    }
  }
  else{
    for (int i=0;i<3;i++){
      dist[i] = part->PositionStart[i] - beampos[i];
      dcos   += part->DirectionStart[i]*beamdir[i];
    }
  }

  if(dist[0] < cuts[0] || dist[0] > cuts[1] || dist[1] < cuts[2] || dist[1] > cuts[3] || dist[2] < cuts[4] || dist[2] > cuts[5] || dcos < cuts[6]) return false;
  else return true;
}


//***************************************************************
std::vector< int > pdAnaUtils::GetPID( const AnaBeamPD& beam, double nominal_momentum ){
//***************************************************************
  const auto& thePIDCands = GetPIDCandidates(beam, nominal_momentum);
  std::vector< int > thePIDs = thePIDCands.getPDGCodes();
  return thePIDs;        
}

//***************************************************************
PossibleParticleCands2 pdAnaUtils::GetPIDCandidates( const AnaBeamPD& beam, double nominal_momentum ){
//***************************************************************
  return GetPIDCandidates_CERNCalib(beam,nominal_momentum);
}
//***************************************************************
PossibleParticleCands2 pdAnaUtils::GetPIDCandidates_CERNCalib( const AnaBeamPD& beam, double nominal_momentum ){
//***************************************************************
  PossibleParticleCands2 candidates;


  bool fUseCERNCalibSelection=true;
  
  //Check if momentum is in valid set
  std::vector< double > valid_momenta = {1., 2., 3., 6., 7.};
  if( std::find(valid_momenta.begin(), valid_momenta.end(), nominal_momentum) == valid_momenta.end() ){
    std::cout << "Reference momentum " << nominal_momentum << " not valid" << std::endl;
    return candidates;
  }
  //Get the high/low pressure Cerenkov info
  //int high_pressure_status, low_pressure_status; 
  int high_pressure_status = beam.CerenkovStatus[0];
  int low_pressure_status  = beam.CerenkovStatus[1];
  //std::cout << "Pressures: " << beam.GetCKov0Pressure() << " " << beam.GetCKov1Pressure() << std::endl;

  //if( beam.GetCKov0Pressure() < beam.GetCKov1Pressure() ){
  //  high_pressure_status = beam.GetCKov1Status();
  //  low_pressure_status = beam.GetCKov0Status();
  //}
  //else{
  //  high_pressure_status = beam.GetCKov0Status();
  //  low_pressure_status = beam.GetCKov1Status();
  //}

  if( nominal_momentum == 1. ){
    if( beam.TOF < 0 ){
      std::cout << "TOF invalid" << std::endl;
      return candidates;
    }
    if( low_pressure_status == -1 ){
      std::cout << "High pressure status invalid" << std::endl;
      return candidates;
    }
    const double & tof = beam.TOF;
    if ( 
        ((fUseCERNCalibSelection && tof < 105.) 
            || (!fUseCERNCalibSelection && tof < 170.))
        && low_pressure_status == 1 
       ) {
      candidates.electron = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof < 110.) 
            || (!fUseCERNCalibSelection && tof < 170.))
        && low_pressure_status == 0 ){
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof > 110. && tof < 160.) 
            || (!fUseCERNCalibSelection && tof > 170.))
        && low_pressure_status == 0 ) {
      candidates.proton = true;
    }
  }
  else if( nominal_momentum == 2. ){
    if( beam.TOF < 0 ){
      std::cout << "TOF invalid" << std::endl;
      return candidates;
    }
    if( low_pressure_status == -1 ){
      std::cout << "High pressure Cerenkov status invalid" << std::endl;
      return candidates;
    }
    const double & tof = beam.TOF;
    if ( 
        ((fUseCERNCalibSelection && tof < 105.) 
            || (!fUseCERNCalibSelection && tof < 160.))
        && low_pressure_status == 1 
       ) {
      candidates.electron = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof < 103.) 
            || (!fUseCERNCalibSelection && tof < 160.))
        && low_pressure_status == 0 ){
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof > 103. && tof < 160.) 
            || (!fUseCERNCalibSelection && tof > 160.))
        && low_pressure_status == 0 ) {
      candidates.proton = true;
    }
  }
  else if( nominal_momentum == 3. ){
    if( high_pressure_status == -1 || low_pressure_status == -1 ){
      std::cout << "At least one Cerenkov status invalid " << std::endl;
      std::cout << "High: " << high_pressure_status << " Low: " << low_pressure_status << std::endl;
      return candidates;
    }
    else if ( low_pressure_status == 1 && high_pressure_status == 1 ) 
      candidates.electron = true;
    else if ( low_pressure_status == 0 && high_pressure_status == 1 ){
      candidates.muon = true;
      candidates.pion = true;
    }
    else{ // low, high = 0, 0
      candidates.proton = true;
      candidates.kaon = true; 
    }
  }
  else if( nominal_momentum == 6. || nominal_momentum == 7. ){
    if( high_pressure_status == -1 || low_pressure_status == -1 ){
      std::cout << "At least one Cerenkov status invalid " << std::endl;
      std::cout << "High: " << high_pressure_status << " Low: " << low_pressure_status << std::endl;
      return candidates;
    }
    else if ( low_pressure_status == 1 && high_pressure_status == 1 ){
      candidates.electron = true;
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( low_pressure_status == 0 && high_pressure_status == 1 ) 
      candidates.kaon = true; 
    else  // low, high = 0, 0
      candidates.proton = true;
  }
  return candidates;

}

//***************************************************************
AnaParticlePD* pdAnaUtils::GetBeamParticle(const AnaEventC& event){
//***************************************************************

  //Get the beam 
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<const AnaEventB*>(&event)->Beam);
  if (!beam->BeamParticle) return NULL;

  // Get the beam particle
  return static_cast<AnaParticlePD*>(beam->BeamParticle);
}

//***************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueBeamParticle(const AnaEventC& event){
//***************************************************************

  // Get the beam particle
  AnaParticlePD* beampart = GetBeamParticle(event);
  if (!beampart) return NULL;

  // Get the true beam particle
  return static_cast<AnaTrueParticlePD*>(beampart->TrueObject);
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromHitPosition(const AnaParticlePD* part){
//***************************************************************

  Float_t length = 0.;

  //check if hit vector is empty
  if(part->Hits[2].empty()){
    //std::cout << "HitPosition vector empty! Returning -1" << std::endl;
    return -1;
  }

  // Initial position
  TVector3 disp(part->Hits[2][0].Position.X(),part->Hits[2][0].Position.Y(),part->Hits[2][0].Position.Z());

  // Add subsequent hits
  for(int i = 1; i < (int)part->Hits[2].size(); ++i){
    if (part->Hits[2][i].Position.X() == -999) break;
    TVector3 pos(part->Hits[2][i].Position.X(),part->Hits[2][i].Position.Y(),part->Hits[2][i].Position.Z());
    disp -= pos;
    length += disp.Mag();
    disp = pos;
  }
  
  return length;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> dEdx){
//***************************************************************

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(dEdx.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = dEdx.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + dEdx[i];
    counter ++;
  }
  
  return accumulated/counter;
}



//***************************************************************
Float_t pdAnaUtils::ComputeCalibrateddQdX(Float_t prim_dqdx, const TVector3& pos){
//***************************************************************

  Float_t hit_x = pos.X();
  Float_t hit_y = pos.Y();
  Float_t hit_z = pos.Z();
	
  if( hit_y < 0. || hit_y > 600. ) return prim_dqdx;
  if( hit_z < 0. || hit_z > 695. ) return prim_dqdx;
	
  Int_t X_bin = X_correction_hist->FindBin( hit_x );
  if (X_bin<1) X_bin = 1;
  if (X_bin>148) X_bin = 148;


  Float_t X_correction = X_correction_hist->GetBinContent(X_bin);
	
  double YZ_correction = (
                          ( hit_x < 0 )
                          ? YZ_neg->GetBinContent( YZ_neg->FindBin( hit_z, hit_y ) ) 
                          : YZ_pos->GetBinContent( YZ_pos->FindBin( hit_z, hit_y ) )  
                          );

  //  Float_t norm_factor = 0.983; // for plane 2
  Float_t norm_factor = 0.9947; // for plane 2
  //  double calib_factor =6.155e-3; //right cali constant for the run 5387. This converts from ADC to e
  double calib_factor =4.57e-3; //right cali constant for the run 5387. This converts from ADC to e
  
  if (debug){
    std::cout << " hit position: " << hit_x << " " << hit_y << " " << hit_z << std::endl;
    std::cout << " prim_dqdx,  X_correction , YZ_correction , norm_factor = "
              << prim_dqdx << " " <<  X_correction << " " <<  YZ_correction << " " <<  norm_factor << std::endl;	
  }
  
  Float_t corrected_dq_dx = prim_dqdx * X_correction * YZ_correction * norm_factor;	
  Float_t scaled_corrected_dq_dx = corrected_dq_dx / calib_factor;		

  
  return scaled_corrected_dq_dx;
}

//***************************************************************
Float_t pdAnaUtils::Compute3DWirePitch(Int_t planeKey, const TVector3& dir){
//***************************************************************

  std::map<int, double> fNormToWiresY;
  std::map<int, double> fNormToWiresZ;

  fNormToWiresY.clear();
  fNormToWiresZ.clear();

  int plane;

  // Numbers from PionAnalyzer_module output
  double dirY_0 = 0.812012;
  double dirZ_0 = 0.58364;

  int NTPC=12;
  plane=0;
  for (int t=0;t<NTPC;t++){
    for (int p=0;p<3;p++){

      double dirY=0;
      double dirZ=0;

      if (p==0){
        if (t%2 == 0) dirZ = -dirZ_0;
        else          dirZ =  dirZ_0;
        dirY = dirY_0;
      }
      else if (p==1){
        if (t%2 == 0) dirZ = -dirZ_0;
        else          dirZ =  dirZ_0;
        dirY = -dirY_0;
      }
      else if (p==2){
        if (t%2 == 0) dirY =  -1;
        else          dirY =   1;
        dirZ = 0;
      }
       
      fNormToWiresY.insert(std::make_pair(plane, -dirZ)); //y component of normal
      fNormToWiresZ.insert(std::make_pair(plane,  dirY)); //z component of normal


      //      std::cout << plane << " --> " << fNormToWiresY[plane] << " " << fNormToWiresZ[plane]  << std::endl;
      
      plane++;
    }
  }


  
  Float_t wirePitch = 0.4792;
      
  //Pitch to use in dEdx calculation
  Float_t yzPitch = wirePitch;   // TODO
  //      geom->WirePitch(hit->WireID().Plane,
  //                      hit->WireID().TPC); //pitch not taking into account angle of track or shower
  Float_t xComponent, pitch3D;
          
  //This assumes equal numbers of TPCs in each cryostat and equal numbers of planes in each TPC
  
  if (fNormToWiresY.count(planeKey) && fNormToWiresZ.count(planeKey)) {
    TVector3 normToWires(0.0, fNormToWiresY.at(planeKey), fNormToWiresZ.at(planeKey));
    yzPitch = wirePitch/fabs(dir.Dot(normToWires));
    //        geom->WirePitch(hit->WireID().Plane, hit->WireID().TPC) / fabs(dir.Dot(normToWires));
  }
  
  xComponent = yzPitch * dir[0] / sqrt(dir[1] * dir[1] + dir[2] * dir[2]);
  pitch3D = sqrt(xComponent * xComponent + yzPitch * yzPitch);///2254.*2742.;

  //  std::cout << yzPitch << " " << xComponent << " " << pitch3D << " " << dir.X() << " " << dir.Y() << " " << dir.Z() << std::endl;
  
  return pitch3D;
}
