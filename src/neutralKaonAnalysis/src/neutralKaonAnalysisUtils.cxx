#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

namespace {

bool GetRecoDaughterTrueObjects(AnaNeutralParticlePD* neutralParticle,
                                AnaTrueParticlePD*& trueDaughter1,
                                AnaTrueParticlePD*& trueDaughter2) {
  trueDaughter1 = nullptr;
  trueDaughter2 = nullptr;

  if(!neutralParticle || !neutralParticle->AnnihilationVertex ||
     neutralParticle->AnnihilationVertex->Particles.size() < 2) {
    return false;
  }

  AnaParticlePD* daughter1 = neutralParticle->AnnihilationVertex->Particles[0];
  AnaParticlePD* daughter2 = neutralParticle->AnnihilationVertex->Particles[1];
  if(!daughter1 || !daughter2) {
    return false;
  }

  trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
  return trueDaughter1 && trueDaughter2;
}

}

//********************************************************************
void neutralKaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddSignalCandidateCategory();
}



//********************************************************************
void neutralKaonAnaUtils::AddSignalCandidateCategory(){

  std::string part_types[] = {
    "signal",
    "background",
    NAMEOTHER};
  int part_codes[]         = {1, 2, CATOTHER};
  int part_colors[]        = {2, 6, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("signal", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
bool neutralKaonAnaUtils::IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    return false;
  }

  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if(!GetRecoDaughterTrueObjects(neutralParticle, trueDaughter1, trueDaughter2)) {
    return false;
  }

  if(trueDaughter1->ParentID <= 0 || trueDaughter1->ParentID != trueDaughter2->ParentID) {
    return false;
  }

  const bool arePions =
      ((trueDaughter1->PDG == 211 && trueDaughter2->PDG == -211) ||
       (trueDaughter1->PDG == -211 && trueDaughter2->PDG == 211));
  if(!arePions) {
    return false;
  }

  AnaTrueParticlePD* trueParent = pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter1->ParentID);
  if(!trueParent) {
    return false;
  }

  return (trueParent->PDG == 310 && trueParent->ProcessEnd == 2);
}

//********************************************************************
void neutralKaonAnaUtils::FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1);
    return;
  }

  const int signalCode = IsSignalCandidate(neutralParticle, event) ? 1 : 2;
  anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
}

