#include "pdAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include <TH3F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/VavilovAccurate.h>

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
