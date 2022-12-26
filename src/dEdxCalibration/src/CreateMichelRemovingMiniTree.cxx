#include "CreateMichelRemovingMiniTree.hxx"
#include "InputManager.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include <iomanip>

//********************************************************************
CreateMichelRemovingMiniTree::CreateMichelRemovingMiniTree(int argc, char *argv[]):CreateMiniTreePD(argc, argv){
//********************************************************************

}

//********************************************************************
bool CreateMichelRemovingMiniTree::Process(){
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

  ClearUninterestingTracks();
  
  if(SaveMiniTree()){
    FillMiniTree();
    // Mark this spill as saved
    _lastSpillSaved=true;
  }
  
  return true;
}

//********************************************************************
bool CreateMichelRemovingMiniTree::SaveMiniTree(){
//********************************************************************
  
  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  if(bunch->Particles.size() > 0)return true;
  else return false;  
}

//********************************************************************
void CreateMichelRemovingMiniTree::ClearUninterestingTracks(){
//********************************************************************
  
  AnaBunchB* bunch = static_cast<AnaBunchB*>(_spill->Bunches[0]);
  std::vector<AnaParticleB*> goodParticles;
  std::vector<AnaParticleB*> badParticles;

  for(int itrk = 0; itrk < (int)bunch->Particles.size(); itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[itrk]);
    ClearUninterestingHits(part);
    if(part->Hits[2].size()>0)goodParticles.push_back(part);
    else                      badParticles.push_back(part);
  }

  for(int itrk = 0; itrk < (int)badParticles.size(); itrk++)
    delete badParticles[itrk];

  bunch->Particles = goodParticles;
}

//********************************************************************
void CreateMichelRemovingMiniTree::ClearUninterestingHits(AnaParticlePD* part){
//********************************************************************
  
  int ihit = 0;
  while(ihit < (int)part->Hits[2].size()){
    if(IsInterestingHit(part->Hits[2][ihit]))
      ihit++;
    else
      part->Hits[2].erase(part->Hits[2].begin()+ihit);
  }
  part->Hits[2].shrink_to_fit();
}

//********************************************************************
bool CreateMichelRemovingMiniTree::IsInterestingHit(AnaHitPD& hit){
//********************************************************************

  bool ItIs = false;

  if(hit.Position.X() > -360 && hit.Position.X() < 0 &&
     abs(hit.Position.Y()-300) < 300 &&
     abs(hit.Position.Z()-350) < 350)
    ItIs = true;

  return ItIs;

}
