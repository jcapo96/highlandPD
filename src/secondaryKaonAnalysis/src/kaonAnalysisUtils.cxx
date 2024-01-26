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

  AddBeamParticleReducedCategory();
  AddBestCandidateParticleReducedCategory();
  AddCandidateParticleReducedCategory();
  AddCandidateDaughterParticleReducedCategory();
  AddCandidateDaughterMuonCategory();
  AddCandidateDaughterMuonReducedCategory();
}

//********************************************************************
void kaonAnaUtils::AddBeamParticleReducedCategory(){
//********************************************************************

  std::string part_types[] = {"#mu^{+}", "e^{+}", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {-13      , -11    , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {7        , 6      , 31       ,  92     ,   8  , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddCategory("beamparticlereduced", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void kaonAnaUtils::AddBestCandidateParticleReducedCategory(){
//********************************************************************

  std::string part_types[] = {"#pi^{-}", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {-211     , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {4        , 31       ,  92     ,   8  , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddCategory("bestcandidateparticlereduced", NPART, part_types, part_codes, part_colors);
}

//********************************************************************
void kaonAnaUtils::AddCandidateParticleReducedCategory(){
//********************************************************************

  std::string part_types[] = {"#pi^{-}", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {-211     , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {4        , 31       ,  92     ,   8  , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("candidateparticlereduced", kaonTree::ncandidates, "ncandidates", 
				      NPART, part_types, part_codes, part_colors, 
				      1, -100);
}

//********************************************************************
void kaonAnaUtils::AddCandidateDaughterParticleReducedCategory(){
//********************************************************************

  std::string part_types[] = {"#pi^{-}","#mu^{+}", "e^{+}", "#pi^{+}","p"  , NAMEOTHER};
  int part_codes[]         = {-211     ,-13      , -11    , 211      , 2212, CATOTHER};
  int part_colors[]        = {4        ,7        , 6      , 31       ,  8  , COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);

  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("candidatedauparticlereduced", kaonTree::ncandidates, "ncandidates", 
				      NPART, part_types, part_codes, part_colors, 
				      1, -100);
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
void kaonAnaUtils::AddCandidateDaughterMuonReducedCategory(){
//********************************************************************

  // --------- candidate's daughter category ----------------
  std::string part_types[] = {"#pi^{-}", "#mu^{+}(k^{+})" , "#mu^{+}(other)", "#pi^{+}", "p"  , NAMEOTHER};
  int part_codes[]         = { -211    , -130             , -132            , 211      ,  2212, CATOTHER};
  int part_colors[]        = { 4       ,  1               , 7               , 31       ,   8  , COLOTHER};

  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("candidatedaumuonreduced", kaonTree::ncandidates, "ncandidates", 
				      NPART, part_types, part_codes, part_colors,
				      1, -100);
}

//********************************************************************
void kaonAnaUtils::FillBeamParticleReducedCategory(AnaParticlePD* beampart){
//********************************************************************
  
  if(!beampart) return;  
  AnaTrueParticle* truepart = static_cast<AnaTrueParticle*>(beampart->TrueObject);

  if(!truepart)                    anaUtils::_categ->SetCode("beamparticlereduced", CATNOTRUTH, CATNOTRUTH);
  else{												          
    if      (truepart->PDG == -13 )anaUtils::_categ->SetCode("beamparticlereduced", -13,        CATOTHER);
    else if (truepart->PDG == -11 )anaUtils::_categ->SetCode("beamparticlereduced", -11,        CATOTHER);
    else if (truepart->PDG == 211 )anaUtils::_categ->SetCode("beamparticlereduced", 211,        CATOTHER);
    else if (truepart->PDG == 321 )anaUtils::_categ->SetCode("beamparticlereduced", 321,        CATOTHER);
    else if (truepart->PDG == 2212)anaUtils::_categ->SetCode("beamparticlereduced", 2212,       CATOTHER);
    else                           anaUtils::_categ->SetCode("beamparticlereduced", CATOTHER,   CATOTHER);
  }
}

//********************************************************************
void kaonAnaUtils::FillBestCandidateParticleReducedCategory(AnaParticlePD* part){
//********************************************************************
  
  if(!part) return;  
  AnaTrueParticle* truepart = static_cast<AnaTrueParticle*>(part->TrueObject);

  if(!truepart)                    anaUtils::_categ->SetCode("bestcandidateparticlereduced", CATNOTRUTH, CATNOTRUTH);
  else{
    if      (truepart->PDG == -211)anaUtils::_categ->SetCode("bestcandidateparticlereduced", -211,       CATOTHER);
    else if (truepart->PDG == 211 )anaUtils::_categ->SetCode("bestcandidateparticlereduced", 211,        CATOTHER);
    else if (truepart->PDG == 321 )anaUtils::_categ->SetCode("bestcandidateparticlereduced", 321,        CATOTHER);
    else if (truepart->PDG == 2212)anaUtils::_categ->SetCode("bestcandidateparticlereduced", 2212,       CATOTHER);
    else                           anaUtils::_categ->SetCode("bestcandidateparticlereduced", CATOTHER,   CATOTHER);
  }
}

//********************************************************************
void kaonAnaUtils::FillCandidateParticleReducedCategory(AnaParticlePD* part){
//********************************************************************
  
  if(!part) return;  
  AnaTrueParticle* truepart = static_cast<AnaTrueParticle*>(part->TrueObject);

  if(!truepart)                    anaUtils::_categ->SetObjectCode("candidateparticlereduced", CATNOTRUTH, CATNOTRUTH, -1);
  else{
         if (truepart->PDG == -211)anaUtils::_categ->SetObjectCode("candidateparticlereduced", -211,       CATOTHER, -1);
    else if (truepart->PDG == 211 )anaUtils::_categ->SetObjectCode("candidateparticlereduced", 211,        CATOTHER, -1);
    else if (truepart->PDG == 321 )anaUtils::_categ->SetObjectCode("candidateparticlereduced", 321,        CATOTHER, -1);
    else if (truepart->PDG == 2212)anaUtils::_categ->SetObjectCode("candidateparticlereduced", 2212,       CATOTHER, -1);
    else                           anaUtils::_categ->SetObjectCode("candidateparticlereduced", CATOTHER,   CATOTHER, -1);
  }
}

//********************************************************************
void kaonAnaUtils::FillCandidateDaughterParticleReducedCategory(AnaParticlePD* part){
//********************************************************************
  
  if(!part) return;  
  AnaTrueParticle* truepart = static_cast<AnaTrueParticle*>(part->TrueObject);

  if(!truepart)                    anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", CATNOTRUTH, CATNOTRUTH, -1);
  else{
         if (truepart->PDG == -211)anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", -211,       CATOTHER, -1);
    else if (truepart->PDG == -13 )anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", -13,        CATOTHER, -1);
    else if (truepart->PDG == -11 )anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", -11,        CATOTHER, -1);
    else if (truepart->PDG == 211 )anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", 211,        CATOTHER, -1);
    else if (truepart->PDG == 2212)anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", 2212,       CATOTHER, -1);
    else                           anaUtils::_categ->SetObjectCode("candidatedauparticlereduced", CATOTHER,   CATOTHER, -1);
  }
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

//********************************************************************
void kaonAnaUtils::FillCandidateDaughterMuonReducedCategory(AnaParticlePD* parent, AnaParticlePD* daughter){
//********************************************************************
  
  if(!parent || !daughter) return;  
  AnaTrueParticle* parentTruePart   = static_cast<AnaTrueParticle*>(parent->TrueObject);
  AnaTrueParticle* daughterTruePart = static_cast<AnaTrueParticle*>(daughter->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the ncandidates counter
  if(!daughterTruePart)                    anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", CATNOTRUTH, CATOTHER, -1);
  else{					
         if (daughterTruePart->PDG == -211)anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", -211,       CATOTHER, -1);
    else if (daughterTruePart->PDG == -13 ){
      if (!parentTruePart)                 anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", -132,       CATOTHER, -1);
      else if (parentTruePart->PDG == 321 )anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", -130,       CATOTHER, -1);
      else                                 anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", -132,       CATOTHER, -1);
    }
    else if (daughterTruePart->PDG == 211 )anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", 211,        CATOTHER, -1);
    else if (daughterTruePart->PDG == 2212)anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", 2212,       CATOTHER, -1);
    else                                   anaUtils::_categ->SetObjectCode("candidatedaumuonreduced", CATOTHER,   CATOTHER, -1);
  }
}
