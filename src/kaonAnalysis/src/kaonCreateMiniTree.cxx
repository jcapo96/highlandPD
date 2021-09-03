#include "kaonCreateMiniTree.hxx"
#include "InputManager.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include <iomanip>

//********************************************************************
kaonCreateMiniTree::kaonCreateMiniTree(int argc, char *argv[]):SimpleLoopBase(argc, argv){
//********************************************************************

  // Add the package version (not yet properly done)
  ND::versioning().AddPackage("kaonAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("KAONANALYSISROOT")));
  
  // add the different converters
  input().AddConverter("minitree", new HighlandMiniTreeConverter("highlandana/MiniTree",true));

}

//********************************************************************
bool kaonCreateMiniTree::Initialize(){
//********************************************************************
  
  _totalEntries=0;
  _savedEntries=0;
  
  _saveCandidates=true;

  return true;
}

//********************************************************************
void kaonCreateMiniTree::DefineOutputTree(){
//********************************************************************

  // Add a new tree called "MiniTree"
  AddTreeWithName(minitree,"MiniTree");

  SetCurrentTree(minitree);
  SetFillSingleTree(GetCurrentTree());

  // The full AnaSpill. _spill cannot be undefined here, otherwise it give problems in some platforms 
  _spill=NULL; 
  GetTree(minitree)->Branch("Spill", "AnaSpill", &_spill,64000,1);
}


//********************************************************************
bool kaonCreateMiniTree::Process(){
//********************************************************************

  // Set the tree to fill
  SetFillSingleTree(minitree);
  SetCurrentTree(minitree);
  
  // Get the corrected spill and set the branch address
  _spill = static_cast<AnaSpill*>(&input().GetCorrectedSpill());
  GetTree(minitree)->SetBranchAddress("Spill",&_spill );  

  // MC
  if(_spill->GetIsMC()){
    // Only considered beam hadrons
    AnaParticleB* BeamPart = static_cast<AnaBeam*>(_spill->Beam)->BeamParticle;
    AnaTrueParticlePD* TrueBeamPart = static_cast<AnaTrueParticlePD*>(BeamPart->TrueObject);
    if(!TrueBeamPart)return true;
    if(abs(TrueBeamPart->PDG)==13 || abs(TrueBeamPart->PDG)==211 || abs(TrueBeamPart->PDG)==2212 || abs(TrueBeamPart->PDG)==321){
      // The spill is saved when there are candidates or when there are true candidates on it (for eff calculation) 
      if(RecoCandidateExists() || TruthCandidateExists()){    
	// Fill the minitree
	FillMiniTree();
	_savedEntries++;
      }
    }
  }
  else{
    // Only beam hadrons considered
    AnaBeamPD* Beam = static_cast<AnaBeamPD*>(_spill->Beam);
    bool hadron = false;
    for(size_t i = 0; i < Beam->PDGs.size(); i++){
      if(Beam->PDGs[i] == 211 || Beam->PDGs[i] == 2212 || Beam->PDGs[i] == 321){
        hadron = true;
        break;
      }
    }
    if(hadron){
      // The spill is saved when there are candidates 
      if(RecoCandidateExists()){    
        // Fill the minitree
        FillMiniTree();
        _savedEntries++;
      }
    }
  }
  
  _totalEntries++;
  
  return true;
}

//********************************************************************
void kaonCreateMiniTree::Finalize(){
//********************************************************************

  // Delete Spills for the last entry
  input().DeleteSpill(); 

  std::cout << std::setprecision(2);
  std::cout << "Filtering factor " << 100*double(_savedEntries)/_totalEntries << "%" << std::endl; 
}

//********************************************************************
void kaonCreateMiniTree::FillMiniTree(){
//********************************************************************

  // Copy arrays into std::vectors
  _spill->CopyArraysIntoVectors();

  // Fill the tree
  FillTree(minitree);
}

//********************************************************************
bool kaonCreateMiniTree::RecoCandidateExists(){
//********************************************************************

  bool itExists = false;
  
  //get the bunch
  AnaBunchPD* Bunch = static_cast<AnaBunchPD*>(_spill->Bunches[0]);

  //loop over particles
  for(size_t i = 1; i < Bunch->Particles.size(); i++){ //skip first particle because it is the primary
    AnaParticlePD* Part = static_cast<AnaParticlePD*>(Bunch->Particles[i]);
    if(Part->DaughtersIDs.size()==1 && Part->BeamOrigin){
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
    if(truePart->PDG==321 || truePart->Origin==4){
      itExists = true;
      break;
    }
  }
  
  return itExists;
}
