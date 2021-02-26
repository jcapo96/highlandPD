#include "kaonAnalysisUtils.hxx"
#include "kaonAnalysis.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void kaonAnaUtils::AddCustomCategories(){
//********************************************************************

  // --------- daughter category ----------------
  
  std::string part_types[] = {"#mu^{-}", "e^{-}", "#pi^{0}", "k^{-}", "#mu^{+}", "e^{+}", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {13       , 11     , 111      , -321   , -13      , -11    , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {2        , 3      , 4        ,  94    , 7        , 6      , 31       ,  92     ,   8  , COLOTHER};

  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("daughter", standardPDTree::seltrk_ndau, "seltrk_ndau", -100, 
                                      NPART, part_types, part_codes, part_colors); 
}

//********************************************************************
void kaonAnaUtils::FillCustomCategories(AnaEventB* event, AnaParticle* part){
//********************************************************************

  if (!part->TrueObject) return;  
  AnaTrueParticlePD*  truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);  
  if (!truePart) return;


  // --------- daughter category ------------------------
  
  for (Int_t i = 0; i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,(Int_t)part->Daughters.size()); i++){
    kaonAnaUtils::FillDaughterCategory(event, static_cast<AnaTrueParticlePD*>(part->Daughters[i]->TrueObject));
  }
}

//********************************************************************
void kaonAnaUtils::FillDaughterCategory(AnaEventB* event, AnaTrueParticlePD* truePart){
//********************************************************************

  if (!truePart) return; 
  
  if      (truePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("daughter", 13  );
  else if (truePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("daughter", 11  );
  else if (truePart->PDG == 111 ) anaUtils::_categ->SetObjectCode("daughter", 111 );
  else if (truePart->PDG == -321) anaUtils::_categ->SetObjectCode("daughter", -321);
  else if (truePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("daughter", -13 );
  else if (truePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("daughter", -11 );
  else if (truePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("daughter", 211 );
  else if (truePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("daughter", 321 );
  else if (truePart->PDG == 2212) anaUtils::_categ->SetObjectCode("daughter", 2212);
  else                            anaUtils::_categ->SetObjectCode("daughter", 10  );

}
