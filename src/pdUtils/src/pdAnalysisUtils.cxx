#include "pdAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "AnalysisUtils.hxx"
#include <TH3F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/VavilovAccurate.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>
#include <TVectorD.h>
#include <TGraph2D.h>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <set>
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


TFile* dEdX_template_file = new TFile( (std::string(getenv("PDUTILSROOT"))+"/data/dEdxrestemplates.root").c_str(), "OPEN" );
std::map< int, TProfile* > templates;

TProfile* ProtonTemplate = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
TProfile* MuonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu" );
TProfile* KaonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_ka" );
/*
templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
*/

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
Float_t pdAnaUtils::ComputeKineticEnergy(const AnaParticlePD &part) {
//********************************************************************

  int plane = 2;

  if(part.Hits[plane].size() < 1)return -1;

  double kinetic = 0;
  double res     = 0;

  for(size_t i = 0; i < part.Hits[plane].size(); i++){
    if(part.Hits[plane][i].dEdx_calib > 1000. || part.Hits[plane][i].dEdx_calib==-999)continue;
    double dedxi = part.Hits[plane][i].dEdx_calib;
    double Residualrangei = part.Hits[plane][i].ResidualRange;
    kinetic = kinetic + dedxi * fabs(Residualrangei - res);
    res = Residualrangei;
  }

  // convert to GeV
  return kinetic/units::GeV;
}

//********************************************************************
void pdAnaUtils::ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance){
//********************************************************************

  distance.clear();
  if (!part) return;
  if (part->Daughters.empty()) return;

  for(size_t i = 0; i < part->Daughters.size(); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);

    //compute distance between part and daughter. Assume dau can be reconstructed backwards
    Float_t dis1 = sqrt(pow(part->PositionEnd[0]-dau->PositionStart[0],2)+pow(part->PositionEnd[1]-dau->PositionStart[1],2)+pow(part->PositionEnd[2]-dau->PositionStart[2],2));
    Float_t dis2 = sqrt(pow(part->PositionEnd[0]-dau->PositionEnd[0],2)+pow(part->PositionEnd[1]-dau->PositionEnd[1],2)+pow(part->PositionEnd[2]-dau->PositionEnd[2],2));

    //return lowest value
    if(dis1<dis2)distance.push_back(dis1);
    else distance.push_back(dis2);
  }
}

//********************************************************************
AnaTrueParticle* pdAnaUtils::GetBeamTrueParticle(const AnaSpillB& spill){
//********************************************************************

  AnaTrueParticle* beampart = NULL;

  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;

  int true_id = 0;
  if(beamPart)
    if(beamPart->TrueObject)
      true_id = static_cast<AnaTrueParticleB*>(beamPart->TrueObject)->ID;

  if(spill.TrueParticles.size() > 0){
    for(int i = 0; i < (int)spill.TrueParticles.size(); i++){
      if(true_id == spill.TrueParticles[i]->ID){
        beampart = static_cast<AnaTrueParticle*>(spill.TrueParticles[i]);
        break;
      }
    }
  }

  return beampart;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(AnaEventB* event, Int_t ID){
//*****************************************************************************

  // Get all reconstructed tracks in the event
  AnaTrueParticleB** trueParticles = event->TrueParticles;
  Int_t nTrueParts                 = event->nTrueParticles;

  for (Int_t i=0;i<nTrueParts;i++){
    if(trueParticles[i]->ID == ID){
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

//*****************************************************************************
AnaParticlePD* pdAnaUtils::GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID){
//*****************************************************************************

  // loop over reconstructed tracks
  for(UInt_t i = 0; i < particles.size(); i++){
    AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(particles[i]->TrueObject);
    if(!truepart)continue;
    if(truepart->ID == true_ID)
      return static_cast<AnaParticlePD*>(particles[i]);
  }

  return NULL;
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID(const AnaParticlePD& part, const int pdg ){
//*****************************************************************************

  double pid_chi2 = 0.;
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;

  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }

  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);

  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
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

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;
      dedx_res *= part.Hits[plane][i].dEdx;


      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);

  return std::make_pair(pid_chi2, npt);
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR){
//*****************************************************************************

  double pid_chi2 = 0.;
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;

  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }

  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);

  //check particles' length is at least as long as maximum RR
  // if(part.Length<RR)
  //   return std::make_pair(9999., -1);
  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    //break whenever above RR upper limit, 26 cm is maximum
    if(part.Hits[plane][i].ResidualRange > RR)continue;

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

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;
      dedx_res *= part.Hits[plane][i].dEdx;


      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

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
  for(size_t i = 1; i < part->Hits[2].size(); ++i){
    if (part->Hits[2][i].Position.X() == -999) break;
    TVector3 pos(part->Hits[2][i].Position.X(),part->Hits[2][i].Position.Y(),part->Hits[2][i].Position.Z());
    disp -= pos;
    length += disp.Mag();
    disp = pos;
  }

  return length;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part){
//***************************************************************

  if(!part)return -999;

  Float_t length = 0.;

  int ntps = part->TrjPoints.size();
  int i0 = 0;
  double x0 = 0,y0 = 0,z0 = 0,dx = 0,dy = 0,dz = 0;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      x0 = part->TrjPoints[itp].Position.X();
      y0 = part->TrjPoints[itp].Position.Y();
      z0 = part->TrjPoints[itp].Position.Z();
      i0 = itp;
      break;
    }
  }

  for(int itp = i0; itp < ntps; itp++){
    if(!part->TrjPoints[itp].IsValid())continue;
    dx = part->TrjPoints[itp].Position.X()-x0;
    dy = part->TrjPoints[itp].Position.Y()-y0;
    dz = part->TrjPoints[itp].Position.Z()-z0;

    length += sqrt(dx*dx+dy*dy+dz*dz);

    x0 = part->TrjPoints[itp].Position.X();
    y0 = part->TrjPoints[itp].Position.Y();
    z0 = part->TrjPoints[itp].Position.Z();
  }

  return length;
}

//***************************************************************
void pdAnaUtils::ComputeParticlePositionAndDirection(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  int ntps = part->TrjPoints.size();
  if(ntps<1)return;

  int ifirst = -1;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      ifirst = itp;
      break;
    }
  }
  if(ifirst != -1){
    part->PositionStart[0] = part->TrjPoints[ifirst].Position.X();
    part->PositionStart[1] = part->TrjPoints[ifirst].Position.Y();
    part->PositionStart[2] = part->TrjPoints[ifirst].Position.Z();
    part->DirectionStart[0] = part->TrjPoints[ifirst].Direction.X();
    part->DirectionStart[1] = part->TrjPoints[ifirst].Direction.Y();
    part->DirectionStart[2] = part->TrjPoints[ifirst].Direction.Z();
  }

  int ilast  = -1;
  for(int itp = 1; itp < ntps; itp++){
    if(part->TrjPoints[ntps-itp].IsValid()){
      ilast = ntps-itp;
      break;
    }
  }
  if(ilast != -1){
    part->PositionEnd[0] = part->TrjPoints[ilast].Position.X();
    part->PositionEnd[1] = part->TrjPoints[ilast].Position.Y();
    part->PositionEnd[2] = part->TrjPoints[ilast].Position.Z();
    part->DirectionEnd[0] = part->TrjPoints[ilast].Direction.X();
    part->DirectionEnd[1] = part->TrjPoints[ilast].Direction.Y();
    part->DirectionEnd[2] = part->TrjPoints[ilast].Direction.Z();
  }
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
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD> hits){
//***************************************************************

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(hits.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = hits.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + hits.at(i).dEdx;
    counter ++;
  }

  return accumulated/counter;
}

//***************************************************************
Float_t pdAnaUtils::ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeDistanceMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(mother->PositionEnd[i] - daughter->PositionStart[i],2);
  dis = sqrt(dis);
  return dis;
}

//***************************************************************
Float_t pdAnaUtils::ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeCosMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + mother->DirectionEnd[i] * daughter->DirectionStart[i];
  return cos;
}

//***************************************************************
Float_t pdAnaUtils::ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange){
//***************************************************************

  if(!part){
    std::cout << "ComputeAveragedEdxOverResRange: particle does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  if(part->Hits[2].empty()){
    //std::cout << "ComputeAveragedEdxOverResRange: has no hits" << std::endl;
    //std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double sumdedx = 0;
  int nhits      = 0;
  for(size_t i = 0; i < part->Hits[2].size(); i++){
    if(part->Hits[2][i].ResidualRange < maxresrange){
      if(part->Hits[2][i].dEdx_calib != -999 && part->Hits[2][i].dEdx_calib < 1000){
        sumdedx += part->Hits[2][i].dEdx_calib;
        nhits++;
      }
    }
  }

  return sumdedx/nhits;
}

//***************************************************************
bool pdAnaUtils::IsStoppingInFV(AnaParticlePD *part){
//***************************************************************

  if(!part)return false;

  bool ItIs = true;

  if((TMath::Abs(part->PositionStart[0])>350 ||
      part->PositionStart[1]<50 || part->PositionStart[1]>550 ||
      part->PositionStart[2]<50 || part->PositionStart[2]>645)
      &&
      (TMath::Abs(part->PositionEnd[0])>350 ||
       part->PositionEnd[1]<50 || part->PositionEnd[1]>550 ||
       part->PositionEnd[2]<50 || part->PositionEnd[2]>645))
    ItIs = false;

  return ItIs;
}

//***************************************************************
int pdAnaUtils::GetHitTPCid(AnaHitPD& hit){
//***************************************************************

  return GetPosTPCid(hit.Position);
}

//***************************************************************
int pdAnaUtils::GetPosTPCid(TVector3 pos){
//***************************************************************

  int TPCid = -1;

  if(pos.X() < 0){
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 1;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 5;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 9;
  }
  else{
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 2;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 6;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 10;
  }

  return TPCid;
}

//***************************************************************
void pdAnaUtils::EstimateHitsDirection(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  //loop over hits
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
    if(ihit < (int)part->Hits[2].size()-1){
      TVector3 dir = part->Hits[2][ihit+1].Position-part->Hits[2][ihit].Position;
      if(dir.Mag()!=0)dir.SetMag(1);
      part->Hits[2][ihit].Direction.SetXYZ(dir.X(),dir.Y(),dir.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(dir.X(),dir.Y(),dir.Z());
    }
    else{
      part->Hits[2][ihit].Direction.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
    }
  }
}

//***************************************************************
void pdAnaUtils::ComputeResidualRange(AnaParticlePD* part){
//***************************************************************

  if(!part)return;

  std::vector<double> delta; delta.clear();
  //loop over hits
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
    TVector3 diff = part->Hits[2][ihit+1].Position-part->Hits[2][ihit].Position;
    delta.push_back(diff.Mag());
  }

  //compute residual range
  std::vector<double> new_rr; new_rr.clear();
  new_rr.push_back(delta[0]/2);
  for(int i = 1; i < (int)delta.size(); i++)
    new_rr.push_back(new_rr[i-1]+delta[i-1]);

  //associate new rr to each hit
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    part->Hits[2][ihit].ResidualRange = new_rr[ihit];
}

//***************************************************************
Double_t pdAnaUtils::ComputeDepositedEnergy(AnaParticlePD* part){
//***************************************************************

  Double_t E = -999;
  if(!part)return E;

  int nhits = part->Hits[2].size();
  if(nhits <= 0)return E;

  E = 0;
  for(int ihit = 1; ihit < nhits-1; ihit++){
    if(part->Hits[2][ihit].dEdx > 1000. || part->Hits[2][ihit].dEdx==-999 || part->Hits[2][ihit].Pitch < 0)
      continue;
    E += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }

  return E;
}

//***************************************************************
Double_t pdAnaUtils::EstimateTrueMomAtAPABorder(AnaParticlePD* part){
//***************************************************************

  Double_t momf = -999;
  if(!part)return momf;

  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return momf;

  if(part->PositionStart[2]>220)return momf;

  int nhits = part->Hits[2].size();
  if(nhits <= 0)return momf;

  double depE = 0;
  for(int ihit = 0; ihit < nhits; ihit++){
    if(part->Hits[2][ihit].Position.Z()>220)break;
    depE += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }

  double mass;
  if(abs(truePart->PDG) == 211)mass       = 139.57;
  else if(abs(truePart->PDG) == 321)mass  = 493.7;
  else if(abs(truePart->PDG) == 2212)mass = 938.27;
  else return momf;

  double momi = truePart->Momentum*1000;
  momf = sqrt(pow(sqrt(momi*momi+mass*mass)-depE,2)-mass*mass);

  return momf;
}

//***************************************************************
Double_t pdAnaUtils::ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts){
//***************************************************************

  double distance = 9999;
  for(int ipart = 0; ipart < nparts; ipart++){
    AnaParticlePD* other = static_cast<AnaParticlePD*>(parts[ipart]);
    if(part->UniqueID == other->UniqueID)continue;
    double dis1 = 0;
    double dis2 = 0;
    for(int idis = 0; idis < 3; idis++){
      dis1 += pow(part->PositionEnd[idis]-other->PositionStart[idis],2);
      dis2 += pow(part->PositionEnd[idis]-other->PositionEnd[idis],2);
    }
    dis1 = sqrt(dis1);
    dis2 = sqrt(dis2);
    if(dis1 < distance)distance = dis1;
    if(dis2 < distance)distance = dis2;
  }

  return distance;
}


//***************************************************************
void pdAnaUtils::GetBeamQualityCuts(AnaEventPD* event,
				    double &mean_x, double &mean_y, double &mean_z,
				    double &sigma_x, double &sigma_y, double &sigma_z,
				    double &cos){
//***************************************************************

  //get nominal beam momentum. If none, set it to 1.
  AnaEventInfoPD* EventInfo = static_cast<AnaEventInfoPD*>(event->EventInfo);
  int NomBeamMom = (int)EventInfo->NominalBeamMom;
  if(NomBeamMom < 0 || NomBeamMom > 3)NomBeamMom = 1; //we still have no values for 6 and 7 GeV

  //get BQC parameters depending on beam mom and MC/data
  std::stringstream ssmom;
  ssmom << NomBeamMom;
  if(event->GetIsMC()){
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.MC."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
  else{
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.Data."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
}

//***************************************************************
double pdAnaUtils::GetDensityCorrection(double beta, double gamma){
//***************************************************************

  //Parameters for the density correction
  const double density_C  = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a  = 0.19559;
  const double density_k  = 3.0;

  //Estimate the density correction
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

//***************************************************************
double pdAnaUtils::GetdEdxBetheBloch(double KE, double mass){
//***************************************************************

  //Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar
  const double I   = 188.0e-6; // [MeV], mean excitation energy
  const double Me  = 0.511; // [Mev], mass of electron

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = GetDensityCorrection(beta, gamma);

  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}

//***************************************************************
double pdAnaUtils::GetWmax(double KE, double mass){
//***************************************************************

  const double Me  = 0.511; // [Mev], mass of electron

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));

  return Wmax;
}

//***************************************************************
double pdAnaUtils::GetLandauXi(double KE, double dx, double mass){
//***************************************************************

  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  return xi;
}

//***************************************************************
double pdAnaUtils::dEdxPDF(double *x, double *par){
//***************************************************************

  ROOT::Math::VavilovAccurate vav;

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
    this_vav = vav.Pdf(y, par[0], par[1]);
    this_vav = this_vav / a;
  }

  return this_vav;
}

//***************************************************************
Float_t pdAnaUtils::dEdxLikelihood(TGraph* tg, TGraph* tg_ke,
				   Float_t mass){
//***************************************************************

  double width = 0.65;
  TF1* pdf = new TF1("pdf", dEdxPDF, -10., 20., 5);
  double likelihood = 0;
  for(int i = 0; i < tg->GetN(); i++){
    double range = tg->GetPointX(i);//+L;
    double dEdx = tg->GetPointY(i);
    double ke = tg_ke->Eval(range);
    double gamma = (ke/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double xi = GetLandauXi(ke, width, mass);
    double Wmax = GetWmax(ke, mass);
    double kappa = xi / Wmax;
    double dEdx_BB = GetdEdxBetheBloch(ke, mass);
    double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
    pdf->SetParameters(par);
    if(pdf->Eval(dEdx) == 0)continue;
    else
      likelihood += log(pdf->Eval(dEdx));
  }
  delete pdf;
  return likelihood;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    if(part->Hits[2][ihit].ResidualRange > maxRR)continue;
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke,
					    Float_t mass){
//***************************************************************

  double width = 0.65;
  TF1* pdf = new TF1("pdf", dEdxPDF, -10., 20., 5);
  double L  = 0;
  double Lf = 10;
  double step = 0.1;
  std::vector<double> L_v,Likelihood_v;
  L_v.clear();
  Likelihood_v.clear();
  while(L<Lf){
    double likelihood = 0;
    for(int i = 0; i < tg->GetN(); i++){
      double range = tg->GetPointX(i)+L;
      double dEdx = tg->GetPointY(i);
      double ke = tg_ke->Eval(range);
      double gamma = (ke/mass)+1.0;
      double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
      double xi = GetLandauXi(ke, width, mass);
      double Wmax = GetWmax(ke, mass);
      double kappa = xi / Wmax;
      double dEdx_BB = GetdEdxBetheBloch(ke, mass);
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
  auto it = std::max_element(Likelihood_v.begin(), Likelihood_v.end());
  int index = std::distance(Likelihood_v.begin(), it);
  delete pdf;
  return std::make_pair(Likelihood_v[index],L_v[index]);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return std::make_pair(-999.,-999.);
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return std::make_pair(-999.,-999.);

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  std::pair<Float_t,Float_t>result = dEdxLikelihoodFreeRange(tg,tg_ke,mass);
  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return std::make_pair(-999.,-999.);
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return std::make_pair(-999.,-999.);

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    if(part->Hits[2][ihit].ResidualRange > maxRR)continue;
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  std::pair<Float_t,Float_t>result = dEdxLikelihoodFreeRange(tg,tg_ke,mass);
  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
// Usage example:
//   std::vector<AnaVertexPD*> vertices = pdAnaUtils::CreateReconstructedVertices(event, 5.0);
//   std::vector<AnaNeutralParticlePD*> neutralParticles = pdAnaUtils::CreateAnaNeutralParticles(event, vertices, 10.0, 2.0);
//   // This will find particles whose end positions are within 10cm of any vertex,
//   // extrapolate their trajectories, and create neutral particles for those with
//   // impact parameter >= 2.0cm (i.e., particles that don't pass too close to the vertex)
//   // A single particle can be assigned to multiple vertices if it meets the criteria for multiple vertices
//***************************************************************
std::vector<AnaNeutralParticlePD*> pdAnaUtils::CreateAnaNeutralParticles(AnaEventB& event, const std::vector<AnaVertexPD*>& vertices, double VertexRadius, double ImpactParameter){
//***************************************************************

  std::vector<AnaNeutralParticlePD*> neutralParticles;

  // Get all particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;
  int neutralParticleID = 0; // Counter for unique neutral particle IDs

  // Loop over all particles
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* particle = static_cast<AnaParticlePD*>(parts[i]);
    if(!particle) continue;

    // Check if particle has valid end position
    if (particle->PositionEnd[0] < -900 || particle->PositionEnd[1] < -900 || particle->PositionEnd[2] < -900) {
      continue; // Skip particles with invalid end positions
    }

    // Extrapolate the particle trajectory first
    double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

    // Check if this particle's end position is within any vertex sphere
    // and if impact parameter is sufficient

    // Loop over all vertices to check if particle end is within sphere
    for(size_t v = 0; v < vertices.size(); v++){
      AnaVertexPD* vertex = vertices[v];
      if(!vertex) continue;

      // Calculate distance from particle end position to vertex position
      TVector3 particleEnd(particle->PositionEnd[0], particle->PositionEnd[1], particle->PositionEnd[2]);
      TVector3 vertexPos(vertex->Position[0], vertex->Position[1], vertex->Position[2]);
      double distance = (particleEnd - vertexPos).Mag();

      // Check if particle end is within sphere
      if(distance < VertexRadius){
        std::vector<double> fitParams = {0,0,0,0,0,0};
        pdAnaUtils::ExtrapolateTrack(particle, fitParams, trackFitLength, false); // Use end position for parent particle
        // Check if extrapolation was successful
        double impactParameter = -999.0;
        if (fitParams[0] == -999.0) {
          // Extrapolation failed, set impact parameter to -999
          impactParameter = -999.0;
        }
        else{
          // Extrapolation successful, calculate impact parameter
          impactParameter = pdAnaUtils::CalculateImpactParameter(fitParams, vertexPos);
        }

        // Create neutral particle for this vertex
        AnaNeutralParticlePD* neutralParticle = new AnaNeutralParticlePD();
        neutralParticle->ImpactParameter = impactParameter;
        neutralParticle->UniqueID = neutralParticleID++;

        // Set vertex and parent
        neutralParticle->Vertex = vertex;
        neutralParticle->Parent = particle;

        neutralParticle->PositionStart[0] = neutralParticle->Parent->PositionEnd[0];
        neutralParticle->PositionStart[1] = neutralParticle->Parent->PositionEnd[1];
        neutralParticle->PositionStart[2] = neutralParticle->Parent->PositionEnd[2];
        neutralParticle->PositionStart[3] = neutralParticle->Parent->PositionEnd[3]; // Time component

        // Set position (vertex position) - using inherited PositionStart from AnaParticleB
        neutralParticle->PositionEnd[0] = vertex->Position[0];
        neutralParticle->PositionEnd[1] = vertex->Position[1];
        neutralParticle->PositionEnd[2] = vertex->Position[2];
        //TODO: Here we might add the time component of the composition of the particles in the vertex
        neutralParticle->PositionEnd[3] = -999; // Time component

        // Calculate direction as the vector from PositionStart to PositionEnd
        Float_t directionStart[3];
        directionStart[0] = neutralParticle->PositionEnd[0] - neutralParticle->PositionStart[0];
        directionStart[1] = neutralParticle->PositionEnd[1] - neutralParticle->PositionStart[1];
        directionStart[2] = neutralParticle->PositionEnd[2] - neutralParticle->PositionStart[2];

        // Normalize the direction vector
        Float_t norm = sqrt(directionStart[0]*directionStart[0] +
                            directionStart[1]*directionStart[1] +
                            directionStart[2]*directionStart[2]);
        if (norm > 0) {
          directionStart[0] /= norm;
          directionStart[1] /= norm;
          directionStart[2] /= norm;
        }

        // Set direction (from extrapolated line) - using inherited DirectionStart from AnaParticleB
        neutralParticle->DirectionStart[0] = directionStart[0];
        neutralParticle->DirectionStart[1] = directionStart[1];
        neutralParticle->DirectionStart[2] = directionStart[2];

        // Calculate DirectionEnd as the sum of DirectionStart of the two particles in the vertex
        if (vertex->NParticles >= 2) {
          AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
          AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

          if (part1 && part2) {
            // Sum the direction vectors of the two particles
            Float_t directionEnd[3];
            directionEnd[0] = part1->DirectionStart[0] + part2->DirectionStart[0];
            directionEnd[1] = part1->DirectionStart[1] + part2->DirectionStart[1];
            directionEnd[2] = part1->DirectionStart[2] + part2->DirectionStart[2];

            // Normalize the resulting direction vector
            Float_t normEnd = sqrt(directionEnd[0]*directionEnd[0] +
                                  directionEnd[1]*directionEnd[1] +
                                  directionEnd[2]*directionEnd[2]);
            if (normEnd > 0) {
              directionEnd[0] /= normEnd;
              directionEnd[1] /= normEnd;
              directionEnd[2] /= normEnd;
            }

            // Set the end direction
            neutralParticle->DirectionEnd[0] = directionEnd[0];
            neutralParticle->DirectionEnd[1] = directionEnd[1];
            neutralParticle->DirectionEnd[2] = directionEnd[2];

            if(part1->ParentID == part2->ParentID){
              // Check if the vertex particles have truth information
              AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(part1->TrueObject);
              AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(part2->TrueObject);

              if(trueParticle1 && trueParticle2) {
                // Use the true parent ID instead of reconstructed ParentID
                if(trueParticle1->ParentID == trueParticle2->ParentID) {
                  // Find the true object associated with the true parent ID
                  AnaTrueParticleB** trueParticles = event.TrueParticles;
                  Int_t nTrueParts = event.nTrueParticles;
                  AnaTrueParticleB* parentTrueObject = nullptr;

                  for (Int_t i = 0; i < nTrueParts; i++) {
                    if (trueParticles[i] && trueParticles[i]->ID == trueParticle1->ParentID) {
                      parentTrueObject = trueParticles[i];
                      break;
                    }
                  }
                  // Find the true particle that corresponds to the neutral particle
                  // It should be the true parent of BOTH vertex particles AND a daughter of the neutral particle's parent
                  AnaTrueParticleB* neutralTrueObject = nullptr;

                  // First, get the parent's true object if it exists
                  if (neutralParticle->Parent && neutralParticle->Parent->TrueObject) {
                    AnaTrueParticleB* parentTrue = static_cast<AnaTrueParticleB*>(neutralParticle->Parent->TrueObject);

                    // Look for a true particle that is both:
                    // 1. The true parent of BOTH vertex particles (ParentID corresponds to true particle ID)
                    // 2. A true daughter of the neutral particle's parent (ParentID corresponds to true particle ID)
                    if (trueParticle1->ParentID == trueParticle2->ParentID) {
                      for (Int_t i = 0; i < nTrueParts; i++) {
                        if (trueParticles[i] &&
                            trueParticles[i]->ID == trueParticle1->ParentID &&  // Is the parent of BOTH vertex particles (by ID)
                            trueParticles[i]->ParentID == parentTrue->ID) {    // Is a daughter of neutral particle's parent (by ID)
                          neutralTrueObject = trueParticles[i];
                          break;
                        }
                      }
                    }
                  }

                  neutralParticle->TrueObject = neutralTrueObject;
                } else {
                  neutralParticle->TrueObject = nullptr;
                }
              } else {
                neutralParticle->TrueObject = nullptr;
              }
            }
            else{
              // The two particles are not from the same parent, so the true object is not set
              neutralParticle->TrueObject = nullptr;
            }
          }
        }

        // Set other properties
        neutralParticle->Mass = -999; // Will be calculated based on particle type
        neutralParticle->Momentum = -999;
        neutralParticle->PDG = -999;
        neutralParticle->Lifetime = -999; // Will be calculated based on particle type
        neutralParticle->DecayLength = -999; // Use impact parameter as decay length

        neutralParticles.push_back(neutralParticle);

      }
    }
    // Note: We don't skip the particle if no vertices found, as it might be assigned to multiple vertices
  }

  return neutralParticles;

}

//********************************************************************
std::vector<AnaVertexPD*> pdAnaUtils::CreateReconstructedVertices(AnaEventB& event, double maxVertexRadius, double maxDaughterDistance){
//********************************************************************

  // Note: maxVertexRadius parameter is not currently used in this implementation
  (void)maxVertexRadius; // Suppress unused parameter warning

  // Get the array of particles from the event
  AnaParticleB** parts = event.Particles;
  int nParts = event.nParticles;

  // Create reconstructed vertices
  std::vector<AnaVertexPD*> reconstructedVertices;
  int vertexID = 0; // Counter for unique vertex IDs

  // Set to track which particle pairs have already been used to create vertices
  std::set<std::pair<AnaParticlePD*, AnaParticlePD*>> usedPairs;

  for(int i = 0; i < nParts; i++){
    AnaParticlePD* daughter1 = static_cast<AnaParticlePD*>(parts[i]);
    if(!daughter1) continue;

    // Check if daughter1 has valid start position
    if (daughter1->PositionStart[0] < -900 || daughter1->PositionStart[1] < -900 || daughter1->PositionStart[2] < -900) {
      continue; // Skip daughter1 with invalid start positions
    }
    // Check if daughter1 has valid end position
    if (daughter1->PositionEnd[0] < -900 || daughter1->PositionEnd[1] < -900 || daughter1->PositionEnd[2] < -900) {
      continue; // Skip daughter1 with invalid end positions
    }

    for(int j = 0; j < nParts; j++){
      AnaParticlePD* daughter2 = static_cast<AnaParticlePD*>(parts[j]);
      if(!daughter2) continue;

      // Check if daughter2 has valid start position
      if (daughter2->PositionStart[0] < -900 || daughter2->PositionStart[1] < -900 || daughter2->PositionStart[2] < -900) {
        continue; // Skip daughter2 with invalid start positions
      }
      // Check if daughter2 has valid end position
      if (daughter2->PositionEnd[0] < -900 || daughter2->PositionEnd[1] < -900 || daughter2->PositionEnd[2] < -900) {
        continue; // Skip daughter2 with invalid end positions
      }

      // Skip if daughter1 and daughter2 are the same particle
      if (daughter1 == daughter2) {
        continue;
      }

      // Create a pair to check for duplicates (order-independent)
      std::pair<AnaParticlePD*, AnaParticlePD*> currentPair;
      if (daughter1 < daughter2) {
        currentPair = std::make_pair(daughter1, daughter2);
      } else {
        currentPair = std::make_pair(daughter2, daughter1);
      }

      // Check if this pair has already been used to create a vertex
      if (usedPairs.find(currentPair) != usedPairs.end()) {
        continue; // Skip if this pair has already been used
      }

      // Check if daughter1 and daughter2 are close enough
      // Use the DefinePosition function to get the position of the particles, this can be changed by the user
      TVector3 pos1 = pdAnaUtils::DefinePosition(daughter1);
      TVector3 pos2 = pdAnaUtils::DefinePosition(daughter2);

      // Check if positions are valid
      if (pos1.X() < -900 || pos2.X() < -900) {
        continue; // Skip if positions are invalid
      }

      // Convert TVector3 to arrays for GetSeparationSquared
      Float_t pos1_array[3] = {static_cast<Float_t>(pos1.X()), static_cast<Float_t>(pos1.Y()), static_cast<Float_t>(pos1.Z())};
      Float_t pos2_array[3] = {static_cast<Float_t>(pos2.X()), static_cast<Float_t>(pos2.Y()), static_cast<Float_t>(pos2.Z())};

      float distance = sqrt(anaUtils::GetSeparationSquared(pos1_array, pos2_array));

      if(distance > maxDaughterDistance){
        continue; // Skip daughter1 and daughter2 if they are not close enough
      }

      // Create reconstructed vertex
      AnaVertexPD* reconstructedVertex = new AnaVertexPD();
      reconstructedVertex->OriginalDistance = distance;
      reconstructedVertex->UniqueID = vertexID++;
      reconstructedVertex->NParticles = 2;
      reconstructedVertex->Particles.push_back(daughter1);
      reconstructedVertex->Particles.push_back(daughter2);

      // Initialize vertex position and minimum distance to invalid values
      reconstructedVertex->Position[0] = -999.0;
      reconstructedVertex->Position[1] = -999.0;
      reconstructedVertex->Position[2] = -999.0;
      reconstructedVertex->MinimumDistance = -999.0;

      // Calculate the vertex position using fitted lines
      pdAnaUtils::FindVertexPosition(reconstructedVertex);

      // Only add vertex if line fitting was successful (position is valid)
      if (reconstructedVertex->Position[0] > -900 && reconstructedVertex->Position[1] > -900 && reconstructedVertex->Position[2] > -900) {
        // Set default values for other properties
        reconstructedVertex->Generation = -999;
        reconstructedVertex->Process = -999;

        reconstructedVertices.push_back(reconstructedVertex);

        // Mark this pair as used to prevent duplicate vertices
        usedPairs.insert(currentPair);
      } else {
        // Delete the vertex if line fitting failed
        delete reconstructedVertex;
      }
    }

  }

  // Return all created vertices without merging
  // The only constraint is that no two vertices can contain the same particle pair
  // (this is already enforced by the usedPairs set in the creation loop above)
  return reconstructedVertices;
}

//***************************************************************
void pdAnaUtils::FindVertexPosition(AnaVertexPD* vertex){
//***************************************************************

  if (!vertex || vertex->NParticles < 2) {
    return;
  }

  // Get the first two particles from the vertex
  AnaParticlePD* part1 = static_cast<AnaParticlePD*>(vertex->Particles[0]);
  AnaParticlePD* part2 = static_cast<AnaParticlePD*>(vertex->Particles[1]);

  if (!part1 || !part2) {
    return;
  }

  // Get track fit length from parameters
  double trackFitLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");

  // Fit lines to both particles
  std::vector<double> line1Params, line2Params;
  pdAnaUtils::ExtrapolateTrack(part1, line1Params, trackFitLength, true); // Use start position for vertex particles
  pdAnaUtils::ExtrapolateTrack(part2, line2Params, trackFitLength, true); // Use start position for vertex particles

  // Check if both fits were successful
  if (line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    return;
  }

  // Find the closest points between the two fitted lines
  TVector3 closestPoint1, closestPoint2;
  double minDistance = pdAnaUtils::FindClosestPointsBetweenLines(line1Params, line2Params,
                                                                 closestPoint1, closestPoint2);

  // Set the vertex position to the midpoint between the closest points
  TVector3 vertexPosition = 0.5 * (closestPoint1 + closestPoint2);

  vertex->Position[0] = vertexPosition.X();
  vertex->Position[1] = vertexPosition.Y();
  vertex->Position[2] = vertexPosition.Z();

  // Store the minimum distance between the fitted lines
  vertex->MinimumDistance = static_cast<Float_t>(minDistance);

  // Store the fitted line parameters for later use in event display
  vertex->FittedLineParams.clear();
  vertex->FittedLineParams.push_back(line1Params);
  vertex->FittedLineParams.push_back(line2Params);

}

//***************************************************************
// Usage example:
//   std::vector<double> fitParams;
//   pdAnaUtils::ExtrapolateTrack(particle, fitParams);
//   // fitParams: [x0, y0, z0, dx, dy, dz] for hits within 15 cm of DefinePosition
//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength, bool useStartPosition){
//***************************************************************

  // Initialize output vector with 6 parameters
  fitParams.clear();
  fitParams.resize(6, -999.0);

  if (!part) {
    return;
  }

  // Get the reference position for distance calculations
  TVector3 referencePos = pdAnaUtils::DefinePosition(part, useStartPosition);
  if (referencePos.X() < -900) {
    return; // Invalid reference position
  }

  // Collect all 3D hit positions from all planes with their distances from reference position
  std::vector<std::pair<TVector3, double>> hitPositionsWithDistance;

  for (int plane = 0; plane < 3; plane++) {
    for (size_t i = 0; i < part->Hits[plane].size(); i++) {
      AnaHitPD& hit = part->Hits[plane][i];
      TVector3 position;

      // Use SCE-corrected position to match event display coordinate system
      if (hit.Position.Z() != -999) {
        position = hit.Position;
      } else {
        continue; // Skip invalid hits
      }

      // Calculate distance from reference position
      double distance = (position - referencePos).Mag();

      hitPositionsWithDistance.push_back(std::make_pair(position, distance));
    }
  }

  // Need at least 2 points to fit a line
  if (hitPositionsWithDistance.size() < 2) {
    return;
  }

  // Fit line to hits within trackLength cm of reference position
  std::vector<TVector3> nearbyHits;
  for (const auto& hitPair : hitPositionsWithDistance) {
    if (hitPair.second <= trackLength) { // Within trackLength cm of reference position
      nearbyHits.push_back(hitPair.first);
    }
  }

  if (nearbyHits.size() >= 2) {
    // Check if ROOT fitting is enabled
    bool useROOTFitting = ND::params().GetParameterI("neutralKaonAnalysis.UseROOTFitting");

    if (useROOTFitting) {
      pdAnaUtils::FitLineToPoints(nearbyHits, fitParams);
    } else {
      // Use original PCA method
      pdAnaUtils::FitLineToPointsPCA(nearbyHits, fitParams);
    }
  }

}

//***************************************************************
// Helper function to fit a line to a set of 3D points using PCA (original method)
//***************************************************************
void pdAnaUtils::FitLineToPointsPCA(const std::vector<TVector3>& points, std::vector<double>& fitParams) {
//***************************************************************

  if (points.size() < 2) {
    return;
  }

  // Calculate centroid of points
  TVector3 centroid(0, 0, 0);
  for (const auto& pos : points) {
    centroid += pos;
  }
  centroid *= (1.0 / points.size());

  // Calculate covariance matrix
  double covXX = 0, covXY = 0, covXZ = 0;
  double covYY = 0, covYZ = 0, covZZ = 0;

  for (const auto& pos : points) {
    double dx = pos.X() - centroid.X();
    double dy = pos.Y() - centroid.Y();
    double dz = pos.Z() - centroid.Z();

    covXX += dx * dx;
    covXY += dx * dy;
    covXZ += dx * dz;
    covYY += dy * dy;
    covYZ += dy * dz;
    covZZ += dz * dz;
  }

  // Normalize by number of points
  int nPoints = points.size();
  covXX /= nPoints;
  covXY /= nPoints;
  covXZ /= nPoints;
  covYY /= nPoints;
  covYZ /= nPoints;
  covZZ /= nPoints;

  // Build covariance matrix
  TMatrixD covMatrix(3, 3);
  covMatrix(0, 0) = covXX; covMatrix(0, 1) = covXY; covMatrix(0, 2) = covXZ;
  covMatrix(1, 0) = covXY; covMatrix(1, 1) = covYY; covMatrix(1, 2) = covYZ;
  covMatrix(2, 0) = covXZ; covMatrix(2, 1) = covYZ; covMatrix(2, 2) = covZZ;

  // Find eigenvalues and eigenvectors
  TMatrixDEigen eigen(covMatrix);
  TVectorD eigenValues = eigen.GetEigenValuesRe();
  TMatrixD eigenVectors = eigen.GetEigenVectors();

  // Find the largest eigenvalue and corresponding eigenvector
  int maxEigenIndex = 0;
  double maxEigenValue = eigenValues[0];
  for (int i = 1; i < 3; i++) {
    if (eigenValues[i] > maxEigenValue) {
      maxEigenValue = eigenValues[i];
      maxEigenIndex = i;
    }
  }

  // The direction vector is the eigenvector corresponding to the largest eigenvalue
  TVector3 direction(eigenVectors(0, maxEigenIndex),
                    eigenVectors(1, maxEigenIndex),
                    eigenVectors(2, maxEigenIndex));
  direction = direction.Unit();

  // Set fit parameters: [x0, y0, z0, dx, dy, dz]
  fitParams[0] = centroid.X();  // x0
  fitParams[1] = centroid.Y();  // y0
  fitParams[2] = centroid.Z();  // z0
  fitParams[3] = direction.X(); // dx
  fitParams[4] = direction.Y(); // dy
  fitParams[5] = direction.Z(); // dz

}

//***************************************************************
// Functor class for 3D line fitting using ROOT
//***************************************************************
class Line3DFitFunction {
public:
  Line3DFitFunction(TGraph2D* graph) : _graph(graph) {}

  double operator()(const double* par) const {
    // Parameters: par[0]=x0, par[1]=y0, par[2]=z0, par[3]=dx, par[4]=dy, par[5]=dz
    TVector3 linePoint(par[0], par[1], par[2]);
    TVector3 lineDirection(par[3], par[4], par[5]);

    // Normalize direction vector
    lineDirection = lineDirection.Unit();

    double sumSquaredDistances = 0.0;

    for (int i = 0; i < _graph->GetN(); i++) {
      double x, y, z;
      _graph->GetPoint(i, x, y, z);
      TVector3 point(x, y, z);

      // Calculate distance from point to line
      // Distance = |(point - linePoint) × lineDirection| / |lineDirection|
      TVector3 diff = point - linePoint;
      TVector3 cross = diff.Cross(lineDirection);
      double distance = cross.Mag(); // lineDirection is already normalized

      sumSquaredDistances += distance * distance;
    }

    return sumSquaredDistances;
  }

private:
  TGraph2D* _graph;
};

//***************************************************************
// Helper function to fit a line to a set of 3D points using ROOT fitting
//***************************************************************
void pdAnaUtils::FitLineToPoints(const std::vector<TVector3>& points, std::vector<double>& fitParams) {
//***************************************************************

  if (points.size() < 2) {
    return;
  }

  // Create a TGraph2D with the 3D points
  TGraph2D* graph = new TGraph2D();
  for (size_t i = 0; i < points.size(); i++) {
    graph->SetPoint(i, points[i].X(), points[i].Y(), points[i].Z());
  }

  // Define the 3D line function: x = x0 + t*dx, y = y0 + t*dy, z = z0 + t*dz
  // We'll minimize the sum of squared distances from points to the line
  ROOT::Fit::Fitter fitter;

  // Create a functor for the distance minimization
  Line3DFitFunction lineFitFunc(graph);
  ROOT::Math::Functor fcn(lineFitFunc, 6); // 6 parameters: x0, y0, z0, dx, dy, dz

  // Set initial parameters using PCA as starting point
  TVector3 centroid(0, 0, 0);
  for (const auto& pos : points) {
    centroid += pos;
  }
  centroid *= (1.0 / points.size());

  // Calculate initial direction using PCA
  TMatrixD covMatrix(3, 3);
  double covXX = 0, covXY = 0, covXZ = 0, covYY = 0, covYZ = 0, covZZ = 0;
  for (const auto& pos : points) {
    double dx = pos.X() - centroid.X();
    double dy = pos.Y() - centroid.Y();
    double dz = pos.Z() - centroid.Z();
    covXX += dx * dx; covXY += dx * dy; covXZ += dx * dz;
    covYY += dy * dy; covYZ += dy * dz; covZZ += dz * dz;
  }
  int nPoints = points.size();
  covMatrix(0, 0) = covXX/nPoints; covMatrix(0, 1) = covXY/nPoints; covMatrix(0, 2) = covXZ/nPoints;
  covMatrix(1, 0) = covXY/nPoints; covMatrix(1, 1) = covYY/nPoints; covMatrix(1, 2) = covYZ/nPoints;
  covMatrix(2, 0) = covXZ/nPoints; covMatrix(2, 1) = covYZ/nPoints; covMatrix(2, 2) = covZZ/nPoints;

  TMatrixDEigen eigen(covMatrix);
  TVectorD eigenValues = eigen.GetEigenValuesRe();
  TMatrixD eigenVectors = eigen.GetEigenVectors();

  int maxEigenIndex = 0;
  double maxEigenValue = eigenValues[0];
  for (int i = 1; i < 3; i++) {
    if (eigenValues[i] > maxEigenValue) {
      maxEigenValue = eigenValues[i];
      maxEigenIndex = i;
    }
  }

  TVector3 initialDirection(eigenVectors(0, maxEigenIndex),
                           eigenVectors(1, maxEigenIndex),
                           eigenVectors(2, maxEigenIndex));
  initialDirection = initialDirection.Unit();

  // Set initial parameters
  double parStart[6] = {centroid.X(), centroid.Y(), centroid.Z(),
                       initialDirection.X(), initialDirection.Y(), initialDirection.Z()};

  fitter.SetFCN(fcn, parStart);

  // Set parameter names and limits
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("z0");
  fitter.Config().ParSettings(3).SetName("dx");
  fitter.Config().ParSettings(4).SetName("dy");
  fitter.Config().ParSettings(5).SetName("dz");

  // Set limits for direction vector components (unit vector)
  fitter.Config().ParSettings(3).SetLimits(-1, 1);
  fitter.Config().ParSettings(4).SetLimits(-1, 1);
  fitter.Config().ParSettings(5).SetLimits(-1, 1);

  // Perform the fit
  bool ok = fitter.FitFCN();

  if (ok) {
    const ROOT::Fit::FitResult& result = fitter.Result();
    const double* parFit = result.GetParams();

    // Normalize direction vector
    TVector3 direction(parFit[3], parFit[4], parFit[5]);
    direction = direction.Unit();

    // Set fit parameters: [x0, y0, z0, dx, dy, dz]
    fitParams[0] = parFit[0];  // x0
    fitParams[1] = parFit[1];  // y0
    fitParams[2] = parFit[2];  // z0
    fitParams[3] = direction.X(); // dx (normalized)
    fitParams[4] = direction.Y(); // dy (normalized)
    fitParams[5] = direction.Z(); // dz (normalized)

  } else {

    // Fallback to PCA method if ROOT fit fails
    TVector3 direction = initialDirection;
    fitParams[0] = centroid.X();  // x0
    fitParams[1] = centroid.Y();  // y0
    fitParams[2] = centroid.Z();  // z0
    fitParams[3] = direction.X(); // dx
    fitParams[4] = direction.Y(); // dy
    fitParams[5] = direction.Z(); // dz
  }

  delete graph;
}

//***************************************************************
// Define the position to use for calculations (distance, line fitting, etc.)
//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle, bool useStartPosition) {
//***************************************************************

  if (!particle) {
    return TVector3(-999, -999, -999);
  }

  // Choose between start and end position based on parameter
  if (useStartPosition) {
    // Use start position (for vertex particles)
    return TVector3(particle->PositionStart[0],
                    particle->PositionStart[1],
                    particle->PositionStart[2]);
  } else {
    // Use end position (for parent particles)
    return TVector3(particle->PositionEnd[0],
                    particle->PositionEnd[1],
                    particle->PositionEnd[2]);
  }

}

//***************************************************************
// Overloaded version for backward compatibility (defaults to start position)
//***************************************************************
TVector3 pdAnaUtils::DefinePosition(AnaParticlePD* particle) {
//***************************************************************
  return DefinePosition(particle, true); // Default to start position for backward compatibility
}

//***************************************************************
// Overloaded version of ExtrapolateTrack for backward compatibility (defaults to start position)
//***************************************************************
void pdAnaUtils::ExtrapolateTrack(AnaParticlePD* part, std::vector<double>& fitParams, double trackLength) {
//***************************************************************
  ExtrapolateTrack(part, fitParams, trackLength, true); // Default to start position for backward compatibility
}

//***************************************************************
// Find the closest points between two 3D lines
//***************************************************************
double pdAnaUtils::FindClosestPointsBetweenLines(const std::vector<double>& line1Params,
                                                const std::vector<double>& line2Params,
                                                TVector3& closestPoint1,
                                                TVector3& closestPoint2) {
//***************************************************************

  // Check if both lines have valid parameters
  if (line1Params.size() != 6 || line2Params.size() != 6 ||
      line1Params[0] == -999.0 || line2Params[0] == -999.0) {
    closestPoint1.SetXYZ(-999, -999, -999);
    closestPoint2.SetXYZ(-999, -999, -999);
    return -999.0;
  }

  // Extract line parameters
  // Line 1: P1(t1) = (x1, y1, z1) + t1 * (dx1, dy1, dz1)
  TVector3 point1(line1Params[0], line1Params[1], line1Params[2]);
  TVector3 dir1(line1Params[3], line1Params[4], line1Params[5]);

  // Line 2: P2(t2) = (x2, y2, z2) + t2 * (dx2, dy2, dz2)
  TVector3 point2(line2Params[0], line2Params[1], line2Params[2]);
  TVector3 dir2(line2Params[3], line2Params[4], line2Params[5]);

  // Normalize direction vectors
  dir1 = dir1.Unit();
  dir2 = dir2.Unit();

  // Vector between the two points
  TVector3 w0 = point1 - point2;

  // Calculate dot products
  double a = dir1.Dot(dir1);  // Should be 1 since dir1 is normalized
  double b = dir1.Dot(dir2);
  double c = dir2.Dot(dir2);  // Should be 1 since dir2 is normalized
  double d = dir1.Dot(w0);
  double e = dir2.Dot(w0);

  // Calculate denominator
  double denom = a * c - b * b;

  // Check if lines are parallel (denominator close to zero)
  if (fabs(denom) < 1e-10) {
    // Lines are parallel, use arbitrary point on line1 and project to line2
    closestPoint1 = point1;
    double t2 = e / c;  // Project point1 onto line2
    closestPoint2 = point2 + t2 * dir2;
  } else {
    // Lines are not parallel, find closest points
    double t1 = (b * e - c * d) / denom;
    double t2 = (a * e - b * d) / denom;

    closestPoint1 = point1 + t1 * dir1;
    closestPoint2 = point2 + t2 * dir2;
  }

  // Calculate the minimum distance
  double minDistance = (closestPoint1 - closestPoint2).Mag();

  return minDistance;
}

//********************************************************************
double pdAnaUtils::CalculateImpactParameter(const std::vector<double>& lineParams, const TVector3& point){
//********************************************************************

  if (lineParams.size() != 6 || lineParams[0] == -999.0) {
    return -999.0;
  }

  TVector3 linePoint(lineParams[0], lineParams[1], lineParams[2]);
  TVector3 lineDirection(lineParams[3], lineParams[4], lineParams[5]);
  lineDirection = lineDirection.Unit(); // Ensure it's normalized

  // Calculate distance from point to line
  TVector3 pointToLine = point - linePoint;
  TVector3 projection = (pointToLine.Dot(lineDirection)) * lineDirection;
  TVector3 perpendicular = pointToLine - projection;

  return perpendicular.Mag();
}