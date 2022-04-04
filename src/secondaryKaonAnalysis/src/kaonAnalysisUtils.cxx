#include "kaonAnalysisUtils.hxx"
#include "secondaryKaonAnalysis.hxx"
#include "kaonTree.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void kaonAnaUtils::AddCustomCategories(){
//********************************************************************

  AddCandidateDaughterMuonCategory();
}

//********************************************************************
void kaonAnaUtils::AddCandidateDaughterMuonCategory(){
//********************************************************************

  // --------- candidate's daughter category ----------------
  std::string part_types[] = {"#mu^{-}", "#pi^{-}", "#mu^{+}(k^{+})", "#mu^{+}(#pi^{+})", "#mu^{+}(other)",  "e^{+}       ", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {13       ,  -211    , -130            , -131              , -132            ,  -11           , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {2        ,  4       ,  1               , 72               , 7               ,  6             , 31       ,  92     ,   8  , COLOTHER};

  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("candidatedaumuon", kaonTree::ncandidates, "ncandidates", 
				      NPART, part_types, part_codes, part_colors,
				      1, -100);
}

//********************************************************************
void kaonAnaUtils::FillCandidateDaughterMuonCategory(AnaParticlePD* parent, AnaParticlePD* daughter){
//********************************************************************
  
  if(!parent || !daughter) return;  
  AnaTrueParticle* parentTruePart   = static_cast<AnaTrueParticle*>(parent->TrueObject);
  AnaTrueParticle* daughterTruePart = static_cast<AnaTrueParticle*>(daughter->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the ncandidates counter
  if(!daughterTruePart)                    anaUtils::_categ->SetObjectCode("candidatedaumuon", CATNOTRUTH, CATOTHER, -1);
  else{												             
    if      (daughterTruePart->PDG == 13  )anaUtils::_categ->SetObjectCode("candidatedaumuon", 13,         CATOTHER, -1);
    else if (daughterTruePart->PDG == -211)anaUtils::_categ->SetObjectCode("candidatedaumuon", -211,       CATOTHER, -1);
    else if (daughterTruePart->PDG == -13 ){
      if (!parentTruePart)                 anaUtils::_categ->SetObjectCode("candidatedaumuon", -132,       CATOTHER, -1);
      else if (parentTruePart->PDG == 321 )anaUtils::_categ->SetObjectCode("candidatedaumuon", -130,       CATOTHER, -1);
      else if (parentTruePart->PDG == 211 )anaUtils::_categ->SetObjectCode("candidatedaumuon", -131,       CATOTHER, -1);
      else                                 anaUtils::_categ->SetObjectCode("candidatedaumuon", -132,       CATOTHER, -1);
    }												             
    else if (daughterTruePart->PDG == -11 )anaUtils::_categ->SetObjectCode("candidatedaumuon", -11,        CATOTHER, -1);
    else if (daughterTruePart->PDG == 211 )anaUtils::_categ->SetObjectCode("candidatedaumuon", 211,        CATOTHER, -1);
    else if (daughterTruePart->PDG == 321 )anaUtils::_categ->SetObjectCode("candidatedaumuon", 321,        CATOTHER, -1);
    else if (daughterTruePart->PDG == 2212)anaUtils::_categ->SetObjectCode("candidatedaumuon", 2212,       CATOTHER, -1);
    else                                   anaUtils::_categ->SetObjectCode("candidatedaumuon", CATOTHER,   CATOTHER, -1);
  }
}
