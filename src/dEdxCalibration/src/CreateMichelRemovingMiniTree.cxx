#include "CreateMichelRemovingMiniTree.hxx"
#include "InputManager.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include <iomanip>

//********************************************************************
CreateMichelRemovingMiniTree::CreateMichelRemovingMiniTree(int argc, char *argv[]):CreateMiniTreePD(argc, argv){
//********************************************************************

  //initialize to default values
  _ThetaXZ_min = 0;
  _ThetaXZ_max = 60;
  _ThetaYZ_min = 0;
  _ThetaYZ_max = 80;
}

//********************************************************************
bool CreateMichelRemovingMiniTree::Initialize(){
//********************************************************************

  // call the base class
  CreateMiniTreePD::Initialize();
  
  _ThetaXZ_min = ND::params().GetParameterI("dEdxCalibration.MichelRemovingMiniTree.ThetaXZ_min");
  _ThetaXZ_max = ND::params().GetParameterI("dEdxCalibration.MichelRemovingMiniTree.ThetaXZ_max");
  _ThetaYZ_min = ND::params().GetParameterI("dEdxCalibration.MichelRemovingMiniTree.ThetaYZ_min");
  _ThetaYZ_max = ND::params().GetParameterI("dEdxCalibration.MichelRemovingMiniTree.ThetaYZ_max");
  
  std::cout << "ANGLE BINNING" << std::endl;
  std::cout << _ThetaXZ_min << " < thetaXZ <" << _ThetaXZ_max << std::endl;
  std::cout << 180-_ThetaXZ_min << " < thetaXZ <" << 180-_ThetaXZ_max << std::endl;
  std::cout << _ThetaYZ_min << " < thetaYZ <" << _ThetaYZ_max << std::endl;
  std::cout << 180-_ThetaYZ_min << " < thetaYZ <" << 180-_ThetaYZ_max << std::endl;

  return true;
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
    if(PartHasGoodAngles(part)){
      ClearUninterestingHits(part);
      if(part->Hits[2].size()>0)goodParticles.push_back(part);
      else                      badParticles.push_back(part);
    }
    else badParticles.push_back(part);;
  }

  for(int itrk = 0; itrk < (int)badParticles.size(); itrk++)
    delete badParticles[itrk];

  bunch->Particles = goodParticles;
}

//********************************************************************
bool CreateMichelRemovingMiniTree::PartHasGoodAngles(AnaParticlePD* part){
//********************************************************************
  
  double xz = abs(180/TMath::Pi()*part->ThetaXZ);
  double yz = abs(180/TMath::Pi()*part->ThetaYZ);
  if(((xz>_ThetaXZ_min && xz<_ThetaXZ_max) || 
      (xz>180-_ThetaXZ_max && xz<180-_ThetaXZ_min)) 
     &&
     ((yz>_ThetaYZ_min && yz<_ThetaYZ_max) || 
      (yz>180-_ThetaYZ_max && yz<180-_ThetaYZ_min)))
    return true;
  else 
    return false;
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
