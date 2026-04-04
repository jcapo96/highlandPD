#include "pdUtilsTrack.hxx"
#include <cmath>
#include <iostream>
#include "TMath.h"

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromHitPosition(const AnaParticlePD* part){
//***************************************************************

  Float_t length = 0.;

  if(part->Hits[2].empty()){
    return -1;
  }

  TVector3 disp(part->Hits[2][0].Position.X(),part->Hits[2][0].Position.Y(),part->Hits[2][0].Position.Z());

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

  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    return -999;
  }

  if(dEdx.empty()){
    return -999;
  }

  int size   = dEdx.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

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

  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    return -999;
  }

  if(hits.empty()){
    return -999;
  }

  int size   = hits.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + hits.at(i).dEdx;
    counter ++;
  }

  return accumulated/counter;
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
    return -999.;
  }

  double sumdedx = 0;
  int nhits      = 0;
  for(size_t i = 0; i < part->Hits[2].size(); i++){
    if(part->Hits[2][i].ResidualRange < maxresrange){
      if(part->Hits[2][i].dEdx != -999 && part->Hits[2][i].dEdx < 1000){
        sumdedx += part->Hits[2][i].dEdx;
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
  const int nh = (int)part->Hits[2].size();
  if (nh < 2) return;
  for (int ihit = 0; ihit < nh - 1; ihit++) {
    TVector3 diff = part->Hits[2][ihit + 1].Position - part->Hits[2][ihit].Position;
    delta.push_back(diff.Mag());
  }

  std::vector<double> new_rr;
  new_rr.reserve(nh);
  new_rr.push_back(delta[0] / 2.);
  for (int ihit = 1; ihit < nh; ++ihit) {
    new_rr.push_back(new_rr[ihit - 1] + delta[ihit - 1]);
  }

  for (int ihit = 0; ihit < nh; ihit++) part->Hits[2][ihit].ResidualRange = new_rr[ihit];
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
