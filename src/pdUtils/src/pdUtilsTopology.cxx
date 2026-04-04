#include "pdUtilsTopology.hxx"
#include <iostream>
#include <cmath>

//********************************************************************
void pdAnaUtils::ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance){
//********************************************************************

  distance.clear();
  if (!part) return;
  if (part->Daughters.empty()) return;

  for(size_t i = 0; i < part->Daughters.size(); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);

    Float_t dis1 = sqrt(pow(part->PositionEnd[0]-dau->PositionStart[0],2)+pow(part->PositionEnd[1]-dau->PositionStart[1],2)+pow(part->PositionEnd[2]-dau->PositionStart[2],2));
    Float_t dis2 = sqrt(pow(part->PositionEnd[0]-dau->PositionEnd[0],2)+pow(part->PositionEnd[1]-dau->PositionEnd[1],2)+pow(part->PositionEnd[2]-dau->PositionEnd[2],2));

    if(dis1<dis2)distance.push_back(dis1);
    else distance.push_back(dis2);
  }
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
