#include "neutralKaonAnalysisUtils.hxx"
#include "neutralKaonAnalysis.hxx"
#include "neutralKaonTree.hxx"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void neutralKaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddNeutralParticleSignalBackgroundCategory();
}

//********************************************************************
void neutralKaonAnaUtils::AddNeutralParticleSignalBackgroundCategory(){
//********************************************************************

  std::string part_types[] = {"signal", "background", NAMEOTHER};
  int part_codes[]         = {1        , 0           , CATOTHER};
  int part_colors[]        = {2        , 1           , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  // name to be called in terminal, counter to follow, variable name
  anaUtils::_categ->AddObjectCategory("isk0", neutralKaonTree::nk0, "nk0",
				      NPART, part_types, part_codes, part_colors,
				      1, -100);
}

//********************************************************************
void neutralKaonAnaUtils::FillNeutralParticleSignalBackgroundCategory(AnaNeutralParticlePD* neutralParticle){
//********************************************************************

  if(!neutralParticle) return;

  // Check if we have a vertex and particles
  if(!neutralParticle->Vertex || neutralParticle->Vertex->Particles.size() < 2) {
    anaUtils::_categ->SetObjectCode("isk0", CATNOTRUTH, CATOTHER, -1);
    return;
  }

  // Get the two particles from the vertex
  AnaParticlePD* particle1 = neutralParticle->Vertex->Particles[0];
  AnaParticlePD* particle2 = neutralParticle->Vertex->Particles[1];

  if(!particle1 || !particle2) {
    anaUtils::_categ->SetObjectCode("isk0", CATNOTRUTH, CATOTHER, -1);
    return;
  }

  // Get true particles
  AnaTrueParticlePD* trueParticle1 = static_cast<AnaTrueParticlePD*>(particle1->TrueObject);
  AnaTrueParticlePD* trueParticle2 = static_cast<AnaTrueParticlePD*>(particle2->TrueObject);

  if(!trueParticle1 || !trueParticle2) {
    anaUtils::_categ->SetObjectCode("isk0", CATNOTRUTH, CATOTHER, -1);
    return;
  }

  // Check if both particles are pi+ and pi- (PDG codes 211 and -211)
  bool isPiPlus1 = (trueParticle1->PDG == 211);
  bool isPiMinus1 = (trueParticle1->PDG == -211);
  bool isPiPlus2 = (trueParticle2->PDG == 211);
  bool isPiMinus2 = (trueParticle2->PDG == -211);

  // Check if we have one pi+ and one pi-
  bool hasPiPlusPiMinus = (isPiPlus1 && isPiMinus2) || (isPiMinus1 && isPiPlus2);

  if(!hasPiPlusPiMinus) {
    anaUtils::_categ->SetObjectCode("isk0", 0, CATOTHER, -1); // background
    return;
  } else {
    anaUtils::_categ->SetObjectCode("isk0", 1, CATOTHER, -1); // signal
  }

//   // Now check if both particles have the same parent K0 (PDG code 310)
//   if(trueParticle1->ParentPDG == 310 && trueParticle2->ParentPDG == 310) {
//     anaUtils::_categ->SetObjectCode("isk0", 1, CATOTHER, -1); // signal
//   } else {
//     anaUtils::_categ->SetObjectCode("isk0", 0, CATOTHER, -1); // background
//   }
}