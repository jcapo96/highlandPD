#include "pdUtilsTruthBeam.hxx"

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

  if(!event) return NULL;

  AnaTrueParticleB** trueParticles = event->TrueParticles;
  Int_t nTrueParts                 = event->nTrueParticles;

  if(!trueParticles || nTrueParts <= 0) return NULL;

  for (Int_t i=0;i<nTrueParts;i++){
    if(trueParticles[i] && trueParticles[i]->ID == ID){
      return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
    }
  }
  return NULL;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID){
//*****************************************************************************

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

  for(UInt_t i = 0; i < particles.size(); i++){
    AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(particles[i]->TrueObject);
    if(!truepart)continue;
    if(truepart->ID == true_ID)
      return static_cast<AnaParticlePD*>(particles[i]);
  }

  return NULL;
}

//*****************************************************************************
bool pdAnaUtils::isBeamLike(AnaParticlePD* part, AnaBeamPD* beam ){
//*****************************************************************************

  if (!part) return false;

  Float_t mccuts[7]  ={-3.,  7., -8.,  7., 27.5, 32.5, 0.93};
  Float_t datacuts[7]={ 0., 10., -5., 10., 30.,  35.0, 0.93};

  AnaParticlePD* beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if (!beampart) return false;

  Float_t beampos[3],beamdir[3], dist[3], dcos=0, cuts[7];

  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  if (trueBeamPart){
    for (int i=0;i<3;i++){
      beampos[i] = trueBeamPart->Position[i]-trueBeamPart->Position[2]*(trueBeamPart->Direction[i]/trueBeamPart->Direction[2]);
      beamdir[i] = trueBeamPart->Direction[i];
    }
    for (int i=0;i<7;i++) cuts[i] = mccuts[i];
  }
  else{
    if(beam->nMomenta != 1 || beam->nTracks != 1)return false;
    for (int i=0;i<3;i++){
      beampos[i] = beampart->PositionEnd[i];
      beamdir[i] = beampart->DirectionEnd[i];
    }
    for (int i=0;i<7;i++) cuts[i] = datacuts[i];
  }

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

  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<const AnaEventB*>(&event)->Beam);
  if (!beam->BeamParticle) return NULL;

  return static_cast<AnaParticlePD*>(beam->BeamParticle);
}

//***************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueBeamParticle(const AnaEventC& event){
//***************************************************************

  AnaParticlePD* beampart = GetBeamParticle(event);
  if (!beampart) return NULL;

  return static_cast<AnaTrueParticlePD*>(beampart->TrueObject);
}
