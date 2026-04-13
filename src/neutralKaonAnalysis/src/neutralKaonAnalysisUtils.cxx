#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

namespace {

bool GetRecoDaughterTrueObjectsFromVertex(AnaAnnihilationVertexPD* vertex,
                                          AnaTrueParticlePD*& trueDaughter1,
                                          AnaTrueParticlePD*& trueDaughter2) {
  trueDaughter1 = nullptr;
  trueDaughter2 = nullptr;
  if (!vertex || vertex->Particles.size() < 2) {
    return false;
  }
  AnaParticlePD* daughter1 = vertex->Particles[0];
  AnaParticlePD* daughter2 = vertex->Particles[1];
  if (!daughter1 || !daughter2) {
    return false;
  }
  trueDaughter1 = static_cast<AnaTrueParticlePD*>(daughter1->TrueObject);
  trueDaughter2 = static_cast<AnaTrueParticlePD*>(daughter2->TrueObject);
  return trueDaughter1 && trueDaughter2;
}

bool GetRecoDaughterTrueObjects(AnaNeutralParticlePD* neutralParticle,
                                AnaTrueParticlePD*& trueDaughter1,
                                AnaTrueParticlePD*& trueDaughter2) {
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) {
    trueDaughter1 = trueDaughter2 = nullptr;
    return false;
  }
  return GetRecoDaughterTrueObjectsFromVertex(neutralParticle->AnnihilationVertex, trueDaughter1, trueDaughter2);
}

AnaTrueParticlePD* GetSharedTrueParentFromVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event) {
  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if (!GetRecoDaughterTrueObjectsFromVertex(vertex, trueDaughter1, trueDaughter2)) {
    return nullptr;
  }
  if (trueDaughter1->ParentID <= 0 || trueDaughter1->ParentID != trueDaughter2->ParentID) {
    return nullptr;
  }
  return pdAnaUtils::GetTrueParticle(const_cast<AnaEventB*>(&event), trueDaughter1->ParentID);
}

bool IsStoppingEndProcess(Int_t proc) {
  return (proc == 2 || proc == 11);
}

int GetSignalStoppingSubtypeCodeFromVertex(AnaAnnihilationVertexPD* vertex) {
  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if (!GetRecoDaughterTrueObjectsFromVertex(vertex, trueDaughter1, trueDaughter2)) {
    return 2;
  }
  const int nStopping = static_cast<int>(IsStoppingEndProcess(trueDaughter1->ProcessEnd)) +
                        static_cast<int>(IsStoppingEndProcess(trueDaughter2->ProcessEnd));
  if (nStopping == 2) return 1;
  if (nStopping == 1) return 5;
  return 6;
}

int GetSignalStoppingSubtypeCode(AnaNeutralParticlePD* neutralParticle) {
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) {
    return 2;
  }
  return GetSignalStoppingSubtypeCodeFromVertex(neutralParticle->AnnihilationVertex);
}

AnaTrueParticlePD* GetSignalTrueParentFromVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event) {
  if (!vertex || vertex->Particles.size() < 2) {
    return nullptr;
  }
  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if (!GetRecoDaughterTrueObjectsFromVertex(vertex, trueDaughter1, trueDaughter2)) {
    return nullptr;
  }
  if (!(trueDaughter1->ParentID > 0 && trueDaughter1->ParentID == trueDaughter2->ParentID)) {
    return nullptr;
  }
  const bool arePions =
      ((trueDaughter1->PDG == 211 && trueDaughter2->PDG == -211) ||
       (trueDaughter1->PDG == -211 && trueDaughter2->PDG == 211));
  if (!arePions) {
    return nullptr;
  }
  AnaTrueParticlePD* trueParent = GetSharedTrueParentFromVertex(vertex, event);
  if (!trueParent) {
    return nullptr;
  }
  if (!(trueParent->PDG == 310 && trueParent->ProcessEnd == 2)) {
    return nullptr;
  }
  return trueParent;
}

bool IsLegitVertexFromTwoBodyDecayVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event) {
  AnaTrueParticlePD* trueParent = GetSharedTrueParentFromVertex(vertex, event);
  if (!trueParent) return false;
  return (trueParent->Daughters.size() == 2);
}

bool IsLegitVertexFromMultiBodyDecayVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event) {
  AnaTrueParticlePD* trueParent = GetSharedTrueParentFromVertex(vertex, event);
  if (!trueParent) return false;
  return (trueParent->Daughters.size() > 2);
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
    "two_stopping",
    "one_stopping",
    "interacting",
    "legit_vertex_2body",
    "legit_vertex_multibody",
    "background",
    NAMEOTHER};
  int part_codes[]         = {1, 5, 6, 3, 4, 2, CATOTHER};
  int part_colors[]        = {2, 8, 6, 4, 7, 46, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("signal", neutralKaonTree::nk0, "nk0",
              NPART, part_types, part_codes, part_colors,
              1, -1000);
}

//********************************************************************
bool neutralKaonAnaUtils::IsLegitVertexCandidate(AnaNeutralParticlePD* neutralParticle){
//********************************************************************

  if(!neutralParticle) {
    return false;
  }

  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if(!GetRecoDaughterTrueObjects(neutralParticle, trueDaughter1, trueDaughter2)) {
    return false;
  }

  return (trueDaughter1->ParentID > 0 && trueDaughter1->ParentID == trueDaughter2->ParentID);
}

//********************************************************************
AnaTrueParticlePD* neutralKaonAnaUtils::GetSignalTrueParent(AnaNeutralParticlePD* neutralParticle,
                                                            const AnaEventB& event){
//********************************************************************

  if (!neutralParticle || !neutralParticle->AnnihilationVertex) {
    return nullptr;
  }
  return GetSignalTrueParentFromVertex(neutralParticle->AnnihilationVertex, event);
}

//********************************************************************
bool neutralKaonAnaUtils::IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  return GetSignalTrueParent(neutralParticle, event) != nullptr;
}

//********************************************************************
bool neutralKaonAnaUtils::IsLegitVertexFromTwoBodyDecay(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) return false;
  return IsLegitVertexFromTwoBodyDecayVertex(neutralParticle->AnnihilationVertex, event);
}

//********************************************************************
bool neutralKaonAnaUtils::IsLegitVertexFromMultiBodyDecay(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) return false;
  return IsLegitVertexFromMultiBodyDecayVertex(neutralParticle->AnnihilationVertex, event);
}

//********************************************************************
void neutralKaonAnaUtils::FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1);
    return;
  }

  int signalCode = 2;
  if(IsSignalCandidate(neutralParticle, event)) {
    signalCode = GetSignalStoppingSubtypeCode(neutralParticle);
  }
  anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
}

//********************************************************************
int neutralKaonAnaUtils::GetSignalCategoryCodeForAnnihilationVertex(AnaAnnihilationVertexPD* vertex,
                                                                   const AnaEventB& event){
//********************************************************************
  if (!vertex) {
    return 2;
  }
  int signalCode = 2;
  if (GetSignalTrueParentFromVertex(vertex, event) != nullptr) {
    signalCode = GetSignalStoppingSubtypeCodeFromVertex(vertex);
  } else if (IsLegitVertexFromTwoBodyDecayVertex(vertex, event)) {
    signalCode = 3;
  } else if (IsLegitVertexFromMultiBodyDecayVertex(vertex, event)) {
    signalCode = 4;
  }
  return signalCode;
}

