#include "CreateMiniTreeCandidates.hxx"
#include "InputManager.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include <iomanip>

//********************************************************************
CreateMiniTreeCandidates::CreateMiniTreeCandidates(int argc, char *argv[]):CreateMiniTreePD(argc, argv){
//********************************************************************

}

//********************************************************************
bool CreateMiniTreeCandidates::Process(){
//********************************************************************

  // Set the tree to fill
  SetFillSingleTree(minitree);
  SetCurrentTree(minitree);
  
  // Get the corrected spill and set the branch address
  _spill = static_cast<AnaSpill*>(&input().GetCorrectedSpill());
  GetTree(minitree)->SetBranchAddress("Spill",&_spill );  

  // The number of POT and Spills since the last saved Spill
  _POTSincePreviousSavedSpill    += _spill->Beam->POTSincePreviousSavedSpill; 
  _SpillsSincePreviousSavedSpill += _spill->Beam->SpillsSincePreviousSavedSpill; 

  _lastSpillSaved = false;

  // Save only candidates
  DeleteAllButCandidates();
  
  // Fill the minitree
  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  if(bunch->Particles.size()>0){
    FillMiniTree();
    // Mark this spill as saved
    _lastSpillSaved=true;
  } 
  
  /*// Fill the minitree
  if(SaveMiniTree()){
    FillMiniTree();
    // Mark this spill as saved
    _lastSpillSaved=true;
    } */ 
  
  return true;
}

//********************************************************************
void CreateMiniTreeCandidates::DeleteAllButCandidates(){
//********************************************************************

  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  
  std::vector<AnaParticleB*> goodParticles;
  std::vector<AnaParticleB*> badParticles; 
  std::set<AnaTrueParticleB*> goodTrueParticles;
  
  goodParticles.clear();
  badParticles.clear();
  goodTrueParticles.clear();

  //Save only selected candidates and its trueobjects, delete all the rest
  for(std::vector<AnaParticleB*>::iterator it = bunch->Particles.begin(); it != bunch->Particles.end(); it++){      
    AnaParticlePD* part = static_cast<AnaParticlePD*>(*it);
    if(IsCandidate(part)){
      goodParticles.push_back(part);
      if(part->TrueObject)goodTrueParticles.insert(static_cast<AnaTrueParticleB*>(part->TrueObject));
    }
    else badParticles.push_back(part);
  }
  
  for(std::vector<AnaParticleB*>::iterator it = badParticles.begin(); it != badParticles.end(); it++){      
    if(std::find(goodParticles.begin(), goodParticles.end(), *it) == goodParticles.end())delete *it;
  }
  
  bunch->Particles = goodParticles;
  _spill->TrueParticles.assign(goodTrueParticles.begin(), goodTrueParticles.end());
}

//********************************************************************
bool CreateMiniTreeCandidates::IsCandidate(AnaParticlePD* part){
//********************************************************************
  
  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);

  //beam filter
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(_spill->Beam);
  bool hadron = false;
  if(_spill->GetIsMC()){
    if(beam->BeamParticle){
      if(beam->BeamParticle->TrueObject){
	AnaTrueParticlePD* beamTruePart = static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);
        if(abs(beamTruePart->PDG)==13  || abs(beamTruePart->PDG)==211 || 
	   abs(beamTruePart->PDG)==321 || abs(beamTruePart->PDG)==2212)hadron = true;
      }
    }
  }
  else{
    for(int i = 0; i < (int)beam->PDGs.size(); i++){
      if(beam->PDGs[i] == 211 || beam->PDGs[i] == 321 || beam->PDGs[i] == 2212 || beam->PDGs[i] == 13){
	hadron = true;
	break;
      }
    }
  }
  if(!hadron)return false;

  //no beam particle
  if(part->isPandora)return false;

  //only one daughter
  if(part->DaughtersIDs.size()!=1)return false;

  //get daughter
  AnaParticlePD* dau = NULL;
  for(int i = 0; i < (int)bunch->Particles.size(); i++){
    if(bunch->Particles[i]->UniqueID == (int)part->DaughtersIDs[0]){
      dau = static_cast<AnaParticlePD*>(bunch->Particles[i]);
      break;
    }
  }
  if(!dau)return false;

  //daughter should be track like
  if(dau->Type != 2)return false;

  //daughter chi2 cut
  if(dau->Chi2Muon/dau->Chi2ndf <= 0.64 || dau->Chi2Muon/dau->Chi2ndf >= 5.84 || dau->Chi2Muon<=0)return false;

  //daughter range mom cut
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(mom <= 0.210 || mom >= 0.237)return false;

  //daughter CNN cut
  //if(abs(dau->CNNscore[0]-0.69)>=0.27)return false;

  //cos cut
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + part->DirectionEnd[i] * dau->DirectionStart[i];
  if(cos <= -0.98 || cos >= 0.64)return false;

  //candidate CNN cut
  //if(part->CNNscore[0]<=0.61)return false;

  //distance cut
  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(part->PositionEnd[i] - dau->PositionStart[i],2);
  dis = sqrt(dis);
  if(dis <= 0 || dis >= 9.9)return false;

  //if we are here we've passed the selection
  return true;
}

//********************************************************************
bool CreateMiniTreeCandidates::SaveMiniTree(){
//********************************************************************
  
  bool save = false;

  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  
  for(std::vector<AnaParticleB*>::iterator it = bunch->Particles.begin(); it != bunch->Particles.end(); it++){      
    AnaParticlePD* part = static_cast<AnaParticlePD*>(*it);
    if(IsCandidate(part)){
      part->IsCandidate = true;
      save = true;
    }
  }

  return save;
}
