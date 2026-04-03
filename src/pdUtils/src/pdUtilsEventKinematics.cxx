#include "pdUtilsEventKinematics.hxx"

#include "AnalysisUtils.hxx"
#include "Parameters.hxx"

#include <cmath>
#include <sstream>
#include <string>

//*****************************************************************************
Float_t pdAnaUtils::ComputeCalorimetricMomentum(AnaParticlePD* part, int pdg, bool includeDecayProducts){
//*****************************************************************************

  if(!part) return -999.0;

  const double pionMass = 0.13957;
  int plane = 2;

  if(part->Hits[plane].empty()) {
    return -999.0;
  }

  double kineticEnergy = 0.0;
  for(size_t i = 0; i < part->Hits[plane].size(); i++){
    if(part->Hits[plane][i].dEdx == -999 || part->Hits[plane][i].Pitch < 0) {
      continue;
    }

    double dedx = part->Hits[plane][i].dEdx;
    double pitch = part->Hits[plane][i].Pitch;
    kineticEnergy += dedx * pitch;
  }

  if(includeDecayProducts && !part->Daughters.empty()) {
    for(size_t i = 0; i < part->Daughters.size(); i++) {
      AnaParticlePD* daughter = static_cast<AnaParticlePD*>(part->Daughters[i]);
      if(!daughter) continue;

      if(!daughter->Hits[plane].empty()) {
        for(size_t j = 0; j < daughter->Hits[plane].size(); j++){
          if(daughter->Hits[plane][j].dEdx == -999 || daughter->Hits[plane][j].Pitch < 0) {
            continue;
          }

          double dedx = daughter->Hits[plane][j].dEdx;
          double pitch = daughter->Hits[plane][j].Pitch;
          kineticEnergy += dedx * pitch;
        }
      }
    }
  }

  double totalEnergy = kineticEnergy/1000 + pionMass;
  double momentumGeV = sqrt(totalEnergy*totalEnergy - pionMass*pionMass);

  double momentum = momentumGeV;
  return (Float_t)momentum;
}

//***************************************************************
void pdAnaUtils::GetBeamQualityCuts(AnaEventPD* event,
                                    double &mean_x, double &mean_y, double &mean_z,
                                    double &sigma_x, double &sigma_y, double &sigma_z,
                                    double &cos){
//***************************************************************

  AnaEventInfoPD* EventInfo = static_cast<AnaEventInfoPD*>(event->EventInfo);
  int NomBeamMom = (int)EventInfo->NominalBeamMom;
  if(NomBeamMom < 0 || NomBeamMom > 3)NomBeamMom = 1;

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

//*****************************************************************************
Float_t pdAnaUtils::ComputeTrueInvariantMass(const AnaTrueParticlePD& part1, const AnaTrueParticlePD& part2, Float_t mass1, Float_t mass2) {

  if (mass1 < 0. || mass2 < 0.)
    return -999.;

  Float_t E1 = sqrt(mass1*mass1 + part1.Momentum * part1.Momentum);
  Float_t E2 = sqrt(mass2*mass2 + part2.Momentum * part2.Momentum);

  Float_t p1x = part1.Momentum * part1.Direction[0];
  Float_t p1y = part1.Momentum * part1.Direction[1];
  Float_t p1z = part1.Momentum * part1.Direction[2];

  Float_t p2x = part2.Momentum * part2.Direction[0];
  Float_t p2y = part2.Momentum * part2.Direction[1];
  Float_t p2z = part2.Momentum * part2.Direction[2];

  Float_t totalE = E1 + E2;
  Float_t totalPx = p1x + p2x;
  Float_t totalPy = p1y + p2y;
  Float_t totalPz = p1z + p2z;
  Float_t totalP2 = totalPx * totalPx + totalPy * totalPy + totalPz * totalPz;

  Float_t M2 = totalE * totalE - totalP2;
  if (M2 > 0) {
    return sqrt(M2);
  }

  return -999.;
}
