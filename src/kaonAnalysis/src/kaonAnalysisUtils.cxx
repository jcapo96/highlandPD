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
  
  // --------- gdaughter category ----------------
  anaUtils::_categ->AddObjectCategory("gdaughter", standardPDTree::seltrk_ngdau, "seltrk_ngdau", -100, 
                                      NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void kaonAnaUtils::FillCustomCategories(AnaEventB* event, AnaParticle* part){
//********************************************************************
  
  // --------- daughter category ------------------------
  
  for (Int_t i = 0; i<std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDDAUGHTERS,(Int_t)part->Daughters.size()); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);
    kaonAnaUtils::FillDaughterCategory(event, dau);

    // --------- gdaughter category ------------------------
    Int_t ngdau = std::min((Int_t)kaonAnalysisConstants::NMAXSAVEDGDAUGHTERS,(Int_t)dau->Daughters.size());
    for (Int_t j = 0; j < ngdau; j++){
      AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[j]);
      kaonAnaUtils::FillGDaughterCategory(event, gdau);
    }
    if(ngdau < (Int_t)kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS){
      AnaParticlePD* dummy = new AnaParticlePD();
      for(Int_t j = ngdau; j < (Int_t)kaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS; j++){
	kaonAnaUtils::FillGDaughterCategory(event,dummy);
      }
      delete dummy;
    }
  }
}

//********************************************************************
void kaonAnaUtils::FillDaughterCategory(AnaEventB* event, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;  
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  
  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ndau counter and
  //the seltrk_ndau counter is incremented despite of the trueObject existing or not
  if(!dauTruePart)                     anaUtils::_categ->SetObjectCode("daughter", 10  );
  else{
    if      (dauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("daughter", 13  );
    else if (dauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("daughter", 11  );
    else if (dauTruePart->PDG == 111 ) anaUtils::_categ->SetObjectCode("daughter", 111 );
    else if (dauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("daughter", -321);
    else if (dauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("daughter", -13 );
    else if (dauTruePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("daughter", -11 );
    else if (dauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("daughter", 211 );
    else if (dauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("daughter", 321 );
    else if (dauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("daughter", 2212);
    else                               anaUtils::_categ->SetObjectCode("daughter", 10  );
  }
    
}

//********************************************************************
void kaonAnaUtils::FillGDaughterCategory(AnaEventB* event, AnaParticlePD* gdau){
//********************************************************************
  
  if(!gdau) return;  
  AnaTrueParticle* gdauTruePart = static_cast<AnaTrueParticle*>(gdau->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ngdau counter and
  //the seltrk_ngdau counter is incremented despite of the trueObject existing or not
  if(!gdauTruePart)                     anaUtils::_categ->SetObjectCode("gdaughter", 10  );
  else{
    if      (gdauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("gdaughter", 13  );
    else if (gdauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("gdaughter", 11  );
    else if (gdauTruePart->PDG == 111 ) anaUtils::_categ->SetObjectCode("gdaughter", 111 );
    else if (gdauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("gdaughter", -321);
    else if (gdauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("gdaughter", -13 );
    else if (gdauTruePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("gdaughter", -11 );
    else if (gdauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("gdaughter", 211 );
    else if (gdauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("gdaughter", 321 );
    else if (gdauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("gdaughter", 2212);
    else                                anaUtils::_categ->SetObjectCode("gdaughter", 10  );
  }

}
