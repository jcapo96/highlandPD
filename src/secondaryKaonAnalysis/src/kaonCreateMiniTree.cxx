#include "kaonCreateMiniTree.hxx"
#include "InputManager.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include <iomanip>

//********************************************************************
kaonCreateMiniTree::kaonCreateMiniTree(int argc, char *argv[]):CreateMiniTree(argc, argv){
//********************************************************************

}

//********************************************************************
bool kaonCreateMiniTree::RecoCandidateExists(){
//********************************************************************

  bool itExists = false;
  
  //get the bunch
  AnaBunchPD* Bunch = static_cast<AnaBunchPD*>(_spill->Bunches[0]);

  //loop over particles
  for(size_t i = 0; i < Bunch->Particles.size(); i++){ 
    AnaParticlePD* Part = static_cast<AnaParticlePD*>(Bunch->Particles[i]);
    if(Part->isPandora)continue;
    if(Part->DaughtersIDs.size()==1 && Part->ParentID!=-1){
      itExists = true;
      break;
    }
  }
  
  return itExists;
}

//********************************************************************
bool kaonCreateMiniTree::TruthCandidateExists(){
//********************************************************************

  bool itExists = false;
  //loop over true particles
  for(size_t i = 1; i < _spill->TrueParticles.size(); i++){ //skip first particle because it is the primary
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(_spill->TrueParticles[i]);
    if(truePart->PDG==321){
      itExists = true;
      break;
    }
  }
  
  return itExists;
}

//********************************************************************
void kaonCreateMiniTree::DeleteUninterestingParticles(){
//********************************************************************

  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  _totalParticles += bunch->Particles.size();

  std::vector<AnaParticleB*> goodParticles;
  std::vector<AnaParticleB*> badParticles; 

  for(std::vector<AnaParticleB*>::iterator it = bunch->Particles.begin(); it != bunch->Particles.end(); it++){      
    AnaParticlePD* part = static_cast<AnaParticlePD*>(*it);
    //save the beam particle
    if(part->isPandora)goodParticles.push_back(static_cast<AnaParticleB*>(part));
    //if it not the beam particle, only save candidates and relatives
    else if(part->Daughters.size() == 1 && part->ParentID != -1){
      AnaParticleB* dau = static_cast<AnaParticleB*>(part->Daughters[0]);
      AnaParticleB* parent = (anaUtils::GetParticleByID(*bunch, part->ParentID));
      if(std::find(goodParticles.begin(), goodParticles.end(), part) == goodParticles.end())
        goodParticles.push_back(static_cast<AnaParticleB*>(part));
      if(dau)
        if(std::find(goodParticles.begin(), goodParticles.end(), dau) == goodParticles.end())
          goodParticles.push_back(dau);
      if(parent)
        if(std::find(goodParticles.begin(), goodParticles.end(), parent) == goodParticles.end())
          goodParticles.push_back(parent); 
    }
    else{
      if(std::find(goodParticles.begin(), goodParticles.end(), part) == goodParticles.end())
        badParticles.push_back(part);
    }
  }
  for(std::vector<AnaParticleB*>::iterator it = badParticles.begin(); it != badParticles.end(); it++){      
    if(std::find(goodParticles.begin(), goodParticles.end(), *it) == goodParticles.end())delete *it;
  }

  bunch->Particles = goodParticles;
  _savedParticles += bunch->Particles.size();
}

//********************************************************************
void kaonCreateMiniTree::DeleteUninterestingTrueParticles(){
//********************************************************************

  // Keep the true particles associated to reconstructed particles and all its descendants
  // TODO: in principle no need to delete the true particles that are not kept since they will be deleted at the end of the event
  //       Since this is the final spill and not the raw spill we cannot delete particles !!!!

  std::set<AnaTrueParticleB*> goodTrueParticles;
  std::vector<AnaTrueParticleB*> badTrueParticles;
  _totalTrueParticles += _spill->TrueParticles.size();


  AnaTrueParticle* beampart = pdAnaUtils::FindBeamTrueParticle(*_spill);
  if(beampart)goodTrueParticles.insert(beampart);
  
  // Loop over all reconstructed particles in the spill      
  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);    
  for(std::vector<AnaParticleB*>::iterator it = bunch->Particles.begin(); it != bunch->Particles.end(); it++){      
    AnaParticleB* part = *it;  
    if(part->TrueObject)
      if(std::find(goodTrueParticles.begin(), goodTrueParticles.end(), part->TrueObject) == goodTrueParticles.end())
        goodTrueParticles.insert(static_cast<AnaTrueParticleB*>(part->TrueObject));
  }
  
  // Loop over all true particles in the spill
  for(std::vector<AnaTrueParticleB*>::iterator it = _spill->TrueParticles.begin(); it != _spill->TrueParticles.end(); it++){      
    AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(*it);
    if(truepart->PDG==321){
      if(std::find(goodTrueParticles.begin(), goodTrueParticles.end(), truepart) == goodTrueParticles.end())
        goodTrueParticles.insert(static_cast<AnaTrueParticleB*>(truepart));
      for(int i = 0; i < (int)truepart->Daughters.size(); i++){
        AnaTrueParticlePD* truedau = pdAnaUtils::GetTrueParticle(_spill->TrueParticles,truepart->Daughters[i]);
        if(truedau)
          if(std::find(goodTrueParticles.begin(), goodTrueParticles.end(), truedau) == goodTrueParticles.end())
            goodTrueParticles.insert(static_cast<AnaTrueParticleB*>(truedau));
      }
      AnaTrueParticlePD* trueparent = pdAnaUtils::GetTrueParticle(_spill->TrueParticles,truepart->ParentID);
      if(trueparent)
        if(std::find(goodTrueParticles.begin(), goodTrueParticles.end(), trueparent) == goodTrueParticles.end())
          goodTrueParticles.insert(static_cast<AnaTrueParticleB*>(trueparent));
    }
  }
  
  // Transfer from the set to the vector
  _spill->TrueParticles.assign(goodTrueParticles.begin(), goodTrueParticles.end());
  _savedTrueParticles += _spill->TrueParticles.size();
}
