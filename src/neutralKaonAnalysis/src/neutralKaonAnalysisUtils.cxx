#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"
#include "TVector3.h"
#include <cmath>

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

// True K0S -> pi+ pi- whose K0S is not in the stored truth (the minitree truth keeps a K0S only when it descends from a
// K+; 96 % of the stored ones do): opposite-charge pions with the same true ParentID and a common true start, whose true
// invariant mass is the K0S mass. The two-body decay fixes it exactly (stored K0S: within 0.003 MeV of 497.611 MeV),
// while pi+ pi- pairs from hadronic interactions populate it at about 25 pairs/MeV in the 12 MC files (about 2.6 pairs
// expected inside the window, against 849 K0S).
const double kK0SMassGeV = 0.497611;
const double kPionMassGeV = 0.13957;
const double kK0SPairMassToleranceGeV = 5e-5;
const double kCommonTrueStartCm = 0.01;

bool IsK0SPionPairByTrueMass(const AnaTrueParticlePD* t1, const AnaTrueParticlePD* t2) {
  if (!(t1->Position[0] > -900.f && t2->Position[0] > -900.f && t1->Momentum > 0.f && t2->Momentum > 0.f)) {
    return false;
  }
  const TVector3 x1(t1->Position[0], t1->Position[1], t1->Position[2]);
  const TVector3 x2(t2->Position[0], t2->Position[1], t2->Position[2]);
  if ((x1 - x2).Mag() > kCommonTrueStartCm) return false;
  TVector3 u1(t1->Direction[0], t1->Direction[1], t1->Direction[2]);
  TVector3 u2(t2->Direction[0], t2->Direction[1], t2->Direction[2]);
  if (u1.Mag2() <= 0. || u2.Mag2() <= 0.) return false;
  const double p1 = t1->Momentum, p2 = t2->Momentum;
  const double e1 = std::sqrt(p1 * p1 + kPionMassGeV * kPionMassGeV);
  const double e2 = std::sqrt(p2 * p2 + kPionMassGeV * kPionMassGeV);
  const double m2 = 2. * kPionMassGeV * kPionMassGeV + 2. * (e1 * e2 - p1 * p2 * u1.Unit().Dot(u2.Unit()));
  return m2 > 0. && std::fabs(std::sqrt(m2) - kK0SMassGeV) < kK0SPairMassToleranceGeV;
}

// Signal decay vertex: its two reco particles are the pi+ and pi- of one true K0S -> pi+ pi- decay, the K0S being
// stored in the truth (GetSignalTrueParentFromVertex) or identified by the true pair mass when it is not.
bool IsSignalPairFromVertex(AnaAnnihilationVertexPD* vertex, const AnaEventB& event) {
  if (GetSignalTrueParentFromVertex(vertex, event)) return true;
  AnaTrueParticlePD* trueDaughter1 = nullptr;
  AnaTrueParticlePD* trueDaughter2 = nullptr;
  if (!GetRecoDaughterTrueObjectsFromVertex(vertex, trueDaughter1, trueDaughter2)) return false;
  if (!(trueDaughter1->ParentID > 0 && trueDaughter1->ParentID == trueDaughter2->ParentID)) return false;
  if (!((trueDaughter1->PDG == 211 && trueDaughter2->PDG == -211) ||
        (trueDaughter1->PDG == -211 && trueDaughter2->PDG == 211))) {
    return false;
  }
  if (GetSharedTrueParentFromVertex(vertex, event)) return false;  // stored parent that is not a decaying K0S
  return IsK0SPionPairByTrueMass(trueDaughter1, trueDaughter2);
}

bool HasNonK0TrueObject(AnaNeutralParticlePD* neutralParticle) {
  if (!neutralParticle) return false;
  AnaTrueParticlePD* neutralTrue = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  if (!neutralTrue) return false;
  return (neutralTrue->PDG != 310 && neutralTrue->PDG != 130);
}

bool HasMatchingRecoParentTruthToVertex(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event) {
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) return false;
  if (!neutralParticle->Parent || !neutralParticle->Parent->TrueObject) return false;

  AnaTrueParticlePD* vertexParent = GetSharedTrueParentFromVertex(neutralParticle->AnnihilationVertex, event);
  if (!vertexParent) return false;

  AnaTrueParticlePD* neutralParentTrue = static_cast<AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
  if (!neutralParentTrue) return false;

  return (vertexParent->ID == neutralParentTrue->ID);
}

bool IsVertexSignalFromAnnihilationVertex(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event) {
  if (!neutralParticle || !neutralParticle->AnnihilationVertex) {
    return false;
  }
  return GetSignalTrueParentFromVertex(neutralParticle->AnnihilationVertex, event) != nullptr;
}

}

//********************************************************************
void neutralKaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddSignalCandidateCategory();
  AddVertexCandidateCategory();
}



//********************************************************************
void neutralKaonAnaUtils::AddSignalCandidateCategory(){

  std::string part_types[] = {
    "two_stopping",
    "one_stopping",
    "interacting",
    "other_true_neutral",
    "interaction_vertex",
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
void neutralKaonAnaUtils::AddVertexCandidateCategory(){

  std::string part_types[] = {
    "signal",
    "background",
    NAMEOTHER};
  int part_codes[]         = {1, 2, CATOTHER};
  int part_colors[]        = {2, 46, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("vertex", neutralKaonTree::nk0, "nk0",
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
bool neutralKaonAnaUtils::IsParentMatchedSignalCandidate(AnaNeutralParticlePD* neutralParticle,
                                                         const AnaEventB& event){
//********************************************************************

  AnaTrueParticlePD* vertexParent = GetSignalTrueParent(neutralParticle, event);
  if (!vertexParent) {
    return false;
  }
  const AnaTrueParticlePD* neutralTrue = static_cast<AnaTrueParticlePD*>(neutralParticle->TrueObject);
  return neutralTrue && neutralTrue->ID == vertexParent->ID;
}

//********************************************************************
bool neutralKaonAnaUtils::IsSignalCandidate(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  return neutralParticle && IsSignalPairFromVertex(neutralParticle->AnnihilationVertex, event);
}

//********************************************************************
void neutralKaonAnaUtils::FillSignalCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  if(!neutralParticle) {
    anaUtils::_categ->SetObjectCode("signal", 2, CATOTHER, -1);
    return;
  }

  int signalCode = 2;

  const bool hasNeutralTrueObject = (neutralParticle->TrueObject != nullptr);

  if (IsSignalCandidate(neutralParticle, event)) {
    signalCode = GetSignalStoppingSubtypeCode(neutralParticle);
  } else if (HasNonK0TrueObject(neutralParticle)) {
    signalCode = 3;
  } else if (!hasNeutralTrueObject && HasMatchingRecoParentTruthToVertex(neutralParticle, event)) {
    signalCode = 4;
  }

  anaUtils::_categ->SetObjectCode("signal", signalCode, CATOTHER, -1);
}

//********************************************************************
void neutralKaonAnaUtils::FillVertexCandidateCategory(AnaNeutralParticlePD* neutralParticle, const AnaEventB& event){
//********************************************************************

  int vertexCode = 2;
  if (IsVertexSignalFromAnnihilationVertex(neutralParticle, event)) {
    vertexCode = 1;
  }

  anaUtils::_categ->SetObjectCode("vertex", vertexCode, CATOTHER, -1);
}

//********************************************************************
int neutralKaonAnaUtils::GetSignalCategoryCodeForAnnihilationVertex(AnaAnnihilationVertexPD* vertex,
                                                                   const AnaEventB& event){
//********************************************************************
  if (!vertex) {
    return 2;
  }
  int signalCode = 2;
  if (IsSignalPairFromVertex(vertex, event)) {
    signalCode = GetSignalStoppingSubtypeCodeFromVertex(vertex);
  }
  return signalCode;
}

