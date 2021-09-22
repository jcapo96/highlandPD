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

  // --------- daughter category ----------------
  
  std::string part_types[] = {"#mu^{-}", "e^{-}", "#pi^{-}" , "k^{-}", "#mu^{+}", "e^{+}", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes[]         = {13       , 11     , -211      , -321   , -13      , -11    , 211      ,  321    ,  2212, CATOTHER};
  int part_colors[]        = {2        , 3      , 4         ,  1     , 7        , 6      , 31       ,  92     ,   8  , COLOTHER};

  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("daughter", standardPDTree::seltrk_ndau, "seltrk_ndau", 
                                      NPART, part_types, part_codes, part_colors,
				      1, -100);
  
  // --------- gdaughter category ----------------
  anaUtils::_categ->AddObjectCategory("gdaughter", standardPDTree::seltrk_ndau, "seltrk_ndau", 
				      NPART, part_types, part_codes, part_colors,
				      2, -100, secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);

  // --------- ggdaughter category ----------------
  anaUtils::_categ->AddObjectCategory("ggdaughter", standardPDTree::seltrk_ndau, "seltrk_ndau", 
				      NPART, part_types, part_codes, part_colors,
				      3, -100, secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS, secondaryKaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS);

  // --------- true secondary daughter category ----------------
  
  std::string part_types_ts[] = {"#mu^{-}", "e^{-}", "#pi^{-}" , "k^{-}", "#mu^{+}", "e^{+}", "#pi^{+}", "TS k^{+}" , "non TS k^{+}" , "p"  , NAMEOTHER};
  int part_codes_ts[]         = {13       , 11     , -211      , -321   , -13      , -11    , 211      ,  3210      ,  3211          ,  2212, CATOTHER};
  int part_colors_ts[]        = {2        , 3      , 4         ,  92    , 7        , 6      , 31       ,  1         ,  11            ,   8  , COLOTHER};

  const int NPART_ts = sizeof(part_types_ts)/sizeof(part_types_ts[0]);
  
  std::reverse(part_types_ts,  part_types_ts  + NPART_ts);
  std::reverse(part_codes_ts,  part_codes_ts  + NPART_ts);
  std::reverse(part_colors_ts, part_colors_ts + NPART_ts);

  anaUtils::_categ->AddObjectCategory("daughterts", standardPDTree::seltrk_ndau, "seltrk_ndau", 
                                      NPART_ts, part_types_ts, part_codes_ts, part_colors_ts,
				      1, -100);

  // --------- gdaughter kaon category ----------------
  std::string part_types_gdkaon[] = {"#mu^{-}", "e^{-}", "#pi^{-}", "k^{-}", "#mu^{+}(k^{+})", "#mu^{+}(#pi^{+})", "#mu^{+}(other)", "e^{+}(k^{+})", "e^{+}       ", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes_gdkaon[]         = {13       , 11     , -211      , -321   , -130            , -131             , -132            , -111          , -11           , 211      ,  321    ,  2212, CATOTHER};
  int part_colors_gdkaon[]        = {2        , 3      , 4         ,  94    , 1               , 72               , 7               , 12            , 6             , 31       ,  92     ,   8  , COLOTHER};

  const int NPART_gdkaon = sizeof(part_types_gdkaon)/sizeof(part_types_gdkaon[0]);
  
  std::reverse(part_types_gdkaon,  part_types_gdkaon  + NPART_gdkaon);
  std::reverse(part_codes_gdkaon,  part_codes_gdkaon  + NPART_gdkaon);
  std::reverse(part_colors_gdkaon, part_colors_gdkaon + NPART_gdkaon);

  anaUtils::_categ->AddObjectCategory("gdaughterkaon", standardPDTree::seltrk_ndau, "seltrk_ndau", 
				      NPART_gdkaon, part_types_gdkaon, part_codes_gdkaon, part_colors_gdkaon,
				      2, -100, secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS);

  anaUtils::_categ->AddObjectCategory("candidatedaumuon", kaonTree::ncandidates, "ncandidates", 
				      NPART_gdkaon, part_types_gdkaon, part_codes_gdkaon, part_colors_gdkaon,
				      1, -100);

  // --------- ggdaughter kaon category ----------------
  std::string part_types_ggdkaon[] = {"#mu^{-}", "e^{-}", "#pi^{-}", "k^{-}" , "#mu^{+}", "e^{+}(k^{+} #rightarrow #mu^{+})", "e^{+}(#pi^{+} #rightarrow #mu^{+})", "e^{+}(other)", "#pi^{+}", "k^{+}" , "p"  , NAMEOTHER};
  int part_codes_ggdkaon[]         = {13       , 11     , -211      , -321   , -13      , -111                              , -112                                , -11           , 211      ,  321    ,  2212, CATOTHER};
  int part_colors_ggdkaon[]        = {2        , 3      , 4         ,  94    , 7        , 1                                 , 72                                  , 6             , 31       ,  92     ,   8  , COLOTHER};

  const int NPART_ggdkaon = sizeof(part_types_ggdkaon)/sizeof(part_types_ggdkaon[0]);
  
  std::reverse(part_types_ggdkaon,  part_types_ggdkaon  + NPART_ggdkaon);
  std::reverse(part_codes_ggdkaon,  part_codes_ggdkaon  + NPART_ggdkaon);
  std::reverse(part_colors_ggdkaon, part_colors_ggdkaon + NPART_ggdkaon);

  anaUtils::_categ->AddObjectCategory("ggdaughterkaon", standardPDTree::seltrk_ndau, "seltrk_ndau", 
				      NPART_ggdkaon, part_types_ggdkaon, part_codes_ggdkaon, part_colors_ggdkaon,
				      3, -100, secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS, secondaryKaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS);

}

//********************************************************************
void kaonAnaUtils::FillCustomCategories(AnaEventB* event, AnaParticle* part){
//********************************************************************

  //dummy anaparticle to fill length-fixed categories
  AnaParticlePD* dummy = new AnaParticlePD();

  // --------- daughter category ------------------------
  
  for (Int_t idau = 0; idau<std::min((Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDDAUGHTERS,(Int_t)part->Daughters.size()); idau++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[idau]);
    kaonAnaUtils::FillDaughterCategory(event, dau);
    kaonAnaUtils::FillDaughterTSCategory(event, static_cast<AnaParticlePD*>(part), dau);

    // --------- gdaughter category ------------------------
    Int_t ngdau = std::min((Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS,(Int_t)dau->Daughters.size());
    for (Int_t jgdau = 0; jgdau < ngdau; jgdau++){
      AnaParticlePD* gdau = static_cast<AnaParticlePD*>(dau->Daughters[jgdau]);
      kaonAnaUtils::FillGDaughterCategory(event, gdau, idau, jgdau);
      kaonAnaUtils::FillGDaughterKaonCategory("gdaughterkaon",event, dau, gdau, idau, jgdau);

      // --------- ggdaughter category ------------------------
      Int_t nggdau = std::min((Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS,(Int_t)gdau->Daughters.size());
      for (Int_t kggdau = 0; kggdau < nggdau; kggdau++){
	AnaParticlePD* ggdau = static_cast<AnaParticlePD*>(gdau->Daughters[kggdau]);
	kaonAnaUtils::FillGGDaughterCategory(event, ggdau, idau, jgdau, kggdau);
	kaonAnaUtils::FillGGDaughterKaonCategory(event, dau, gdau, ggdau, idau, jgdau, kggdau);
      }

      // as the variables in the three have a fixed length NMAXSAVEDGGDAUGHTERS, the category has to be filled even if the
      // ggdaughter does not exist
      for (Int_t kggdau = nggdau; kggdau < (Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS; kggdau++){
	kaonAnaUtils::FillGGDaughterCategory(event, dummy, idau, jgdau, kggdau);
	kaonAnaUtils::FillGGDaughterKaonCategory(event, dummy, dummy, dummy, idau, jgdau, kggdau);
      }
    }
    
    // as the variables in the three have a fixed length NMAXSAVEDGDAUGHTERS, the category has to be filled even if the
    // gdaughter does not exist
    for (Int_t jgdau = ngdau; jgdau < (Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDGDAUGHTERS; jgdau++){
      kaonAnaUtils::FillGDaughterCategory(event, dummy, idau, jgdau);
      kaonAnaUtils::FillGDaughterKaonCategory("gdaughterkaon",event, dummy, dummy, idau, jgdau);
      for (Int_t kggdau = 0; kggdau < (Int_t)secondaryKaonAnalysisConstants::NMAXSAVEDGGDAUGHTERS; kggdau++){
	kaonAnaUtils::FillGGDaughterCategory(event, dummy, idau, jgdau, kggdau);
	kaonAnaUtils::FillGGDaughterKaonCategory(event, dummy, dummy, dummy, idau, jgdau, kggdau);
      }
    }
  }

  delete dummy;
}

//********************************************************************
void kaonAnaUtils::FillDaughterCategory(AnaEventB* event, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;  
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  
  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ndau counter and
  //the seltrk_ndau counter is incremented despite of the trueObject existing or not
  if(!dauTruePart)                     anaUtils::_categ->SetObjectCode("daughter", CATNOTRUTH, CATOTHER);
  else{
    if      (dauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("daughter", 13,         CATOTHER);
    else if (dauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("daughter", 11,         CATOTHER);
    else if (dauTruePart->PDG == -211) anaUtils::_categ->SetObjectCode("daughter", -211,       CATOTHER);
    else if (dauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("daughter", -321,       CATOTHER);
    else if (dauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("daughter", -13,        CATOTHER);
    else if (dauTruePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("daughter", -11,        CATOTHER);
    else if (dauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("daughter", 211,        CATOTHER);
    else if (dauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("daughter", 321,        CATOTHER);
    else if (dauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("daughter", 2212,       CATOTHER);
    else                               anaUtils::_categ->SetObjectCode("daughter", CATOTHER,   CATOTHER);
  }
    
}

//********************************************************************
void kaonAnaUtils::FillGDaughterCategory(AnaEventB* event, AnaParticlePD* gdau, Int_t indx1, Int_t indx2){
//********************************************************************
  
  if(!gdau) return;  
  AnaTrueParticle* gdauTruePart = static_cast<AnaTrueParticle*>(gdau->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ngdau counter and
  //the seltrk_ngdau counter is incremented despite of the trueObject existing or not
  if(!gdauTruePart)                     anaUtils::_categ->SetObjectCode("gdaughter", CATNOTRUTH, CATOTHER, indx2);
  else{												         
    if      (gdauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("gdaughter", 13,         CATOTHER, indx2);
    else if (gdauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("gdaughter", 11,         CATOTHER, indx2);
    else if (gdauTruePart->PDG == -211) anaUtils::_categ->SetObjectCode("gdaughter", -211,       CATOTHER, indx2);
    else if (gdauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("gdaughter", -321,       CATOTHER, indx2);
    else if (gdauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("gdaughter", -13,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("gdaughter", -11,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("gdaughter", 211,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("gdaughter", 321,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("gdaughter", 2212,       CATOTHER, indx2);
    else                                anaUtils::_categ->SetObjectCode("gdaughter", CATOTHER,   CATOTHER, indx2);
  }

}

//********************************************************************
void kaonAnaUtils::FillGGDaughterCategory(AnaEventB* event, AnaParticlePD* ggdau, Int_t indx1, Int_t indx2, Int_t indx3){
//********************************************************************
  
  if(!ggdau) return;  
  AnaTrueParticle* ggdauTruePart = static_cast<AnaTrueParticle*>(ggdau->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ngdau counter and
  //the seltrk_ngdau counter is incremented despite of the trueObject existing or not
  if(!ggdauTruePart)                     anaUtils::_categ->SetObjectCode("ggdaughter", CATNOTRUTH, CATOTHER, indx2, indx3);
  else{										       		              
    if      (ggdauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("ggdaughter", 13,         CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("ggdaughter", 11,         CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -211) anaUtils::_categ->SetObjectCode("ggdaughter", -211,       CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("ggdaughter", -321,       CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("ggdaughter", -13,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -11 ) anaUtils::_categ->SetObjectCode("ggdaughter", -11,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("ggdaughter", 211,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("ggdaughter", 321,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("ggdaughter", 2212,       CATOTHER, indx2, indx3);
    else                                 anaUtils::_categ->SetObjectCode("ggdaughter", CATOTHER,   CATOTHER, indx2, indx3);
  }
}

//********************************************************************
void kaonAnaUtils::FillDaughterTSCategory(AnaEventB* event, AnaParticlePD* part, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;  
  AnaTrueParticle* truePart    = static_cast<AnaTrueParticle*>(part->TrueObject);
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  
  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ndau counter and
  //the seltrk_ndau counter is incremented despite of the trueObject existing or not
  if(!dauTruePart)                                 anaUtils::_categ->SetObjectCode("daughterts", CATNOTRUTH, CATOTHER);
  else{
    if      (dauTruePart->PDG == 13  )             anaUtils::_categ->SetObjectCode("daughterts", 13,         CATOTHER);
    else if (dauTruePart->PDG == 11  )             anaUtils::_categ->SetObjectCode("daughterts", 11,         CATOTHER);
    else if (dauTruePart->PDG == -211)             anaUtils::_categ->SetObjectCode("daughterts", -211,       CATOTHER);
    else if (dauTruePart->PDG == -321)             anaUtils::_categ->SetObjectCode("daughterts", -321,       CATOTHER);
    else if (dauTruePart->PDG == -13 )             anaUtils::_categ->SetObjectCode("daughterts", -13,        CATOTHER);
    else if (dauTruePart->PDG == -11 )             anaUtils::_categ->SetObjectCode("daughterts", -11,        CATOTHER);
    else if (dauTruePart->PDG == 211 )             anaUtils::_categ->SetObjectCode("daughterts", 211,        CATOTHER);
    else if (dauTruePart->PDG == 321 ){
      if(!truePart)anaUtils::_categ->SetObjectCode("daughterts", 3211, CATOTHER);
      else if (truePart->ParentID==0 && 
	       truePart->ID==dauTruePart->ParentID)anaUtils::_categ->SetObjectCode("daughterts", 3210,       CATOTHER);
      else                                         anaUtils::_categ->SetObjectCode("daughterts", 3211,       CATOTHER);
    }
    else if (dauTruePart->PDG == 2212)             anaUtils::_categ->SetObjectCode("daughterts", 2212,       CATOTHER);
    else                                           anaUtils::_categ->SetObjectCode("daughterts", CATOTHER,   CATOTHER);
  }
    
}


//********************************************************************
void kaonAnaUtils::FillGDaughterKaonCategory(std::string catname, AnaEventB* event, AnaParticlePD* dau, AnaParticlePD* gdau, Int_t indx1, Int_t indx2){
//********************************************************************
  
  if(!dau || !gdau) return;  
  AnaTrueParticle* dauTruePart  = static_cast<AnaTrueParticle*>(dau->TrueObject);
  AnaTrueParticle* gdauTruePart = static_cast<AnaTrueParticle*>(gdau->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ngdau counter and
  //the seltrk_ngdau counter is incremented despite of the trueObject existing or not
  if(!gdauTruePart)                     anaUtils::_categ->SetObjectCode(catname, CATNOTRUTH, CATOTHER, indx2);
  else{												             
    if      (gdauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode(catname, 13,         CATOTHER, indx2);
    else if (gdauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode(catname, 11,         CATOTHER, indx2);
    else if (gdauTruePart->PDG == -211) anaUtils::_categ->SetObjectCode(catname, -211,       CATOTHER, indx2);
    else if (gdauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode(catname, -321,       CATOTHER, indx2);
    else if (gdauTruePart->PDG == -13 ){
      if (!dauTruePart)                 anaUtils::_categ->SetObjectCode(catname, -132,       CATOTHER, indx2);
      else if (dauTruePart->PDG == 321) anaUtils::_categ->SetObjectCode(catname, -130,       CATOTHER, indx2);
      else if (dauTruePart->PDG == 211) anaUtils::_categ->SetObjectCode(catname, -131,       CATOTHER, indx2);
      else                              anaUtils::_categ->SetObjectCode(catname, -132,       CATOTHER, indx2);
    }												             
    else if (gdauTruePart->PDG == -11 ){
      if (!dauTruePart)                 anaUtils::_categ->SetObjectCode(catname, -11,        CATOTHER, indx2);
      else if (dauTruePart->PDG == 321) anaUtils::_categ->SetObjectCode(catname, -111,       CATOTHER, indx2);
      else                              anaUtils::_categ->SetObjectCode(catname, -11,        CATOTHER, indx2);
    }
    else if (gdauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode(catname, 211,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode(catname, 321,        CATOTHER, indx2);
    else if (gdauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode(catname, 2212,       CATOTHER, indx2);
    else                                anaUtils::_categ->SetObjectCode(catname, CATOTHER,   CATOTHER, indx2);
  }
}

//********************************************************************
void kaonAnaUtils::FillGGDaughterKaonCategory(AnaEventB* event, AnaParticlePD* dau, AnaParticlePD* gdau, AnaParticlePD* ggdau, Int_t indx1, Int_t indx2, Int_t indx3){
//********************************************************************
  
  if(!dau || !gdau || !ggdau) return;  
  AnaTrueParticle* dauTruePart   = static_cast<AnaTrueParticle*>(dau->TrueObject);
  AnaTrueParticle* gdauTruePart  = static_cast<AnaTrueParticle*>(gdau->TrueObject);
  AnaTrueParticle* ggdauTruePart = static_cast<AnaTrueParticle*>(ggdau->TrueObject);

  //if there is no trueObject associated to the object we can't avoid adding the
  //object to the category, since the category follows the seltrk_ngdau counter and
  //the seltrk_ngdau counter is incremented despite of the trueObject existing or not
  if(!ggdauTruePart)                     anaUtils::_categ->SetObjectCode("ggdaughterkaon", CATNOTRUTH, CATOTHER, indx2, indx3);
  else{
    if      (ggdauTruePart->PDG == 13  ) anaUtils::_categ->SetObjectCode("ggdaughterkaon", 13,         CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 11  ) anaUtils::_categ->SetObjectCode("ggdaughterkaon", 11,         CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -211) anaUtils::_categ->SetObjectCode("ggdaughterkaon", -211,       CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -321) anaUtils::_categ->SetObjectCode("ggdaughterkaon", -321,       CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -13 ) anaUtils::_categ->SetObjectCode("ggdaughterkaon", -13,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == -11 ){							      
      if(!dauTruePart || !gdauTruePart)  anaUtils::_categ->SetObjectCode("ggdaughterkaon", -11,        CATOTHER, indx2, indx3);
      else if(dauTruePart->PDG  == 321 &&							      
	      gdauTruePart->PDG == -13)  anaUtils::_categ->SetObjectCode("ggdaughterkaon", -111,       CATOTHER, indx2, indx3);
      else if(dauTruePart->PDG  == 211 &&							      
	      gdauTruePart->PDG == -13)  anaUtils::_categ->SetObjectCode("ggdaughterkaon", -112,       CATOTHER, indx2, indx3);
      else                               anaUtils::_categ->SetObjectCode("ggdaughterkaon", -11,        CATOTHER, indx2, indx3);
    }												      
    else if (ggdauTruePart->PDG == 211 ) anaUtils::_categ->SetObjectCode("ggdaughterkaon", 211,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 321 ) anaUtils::_categ->SetObjectCode("ggdaughterkaon", 321,        CATOTHER, indx2, indx3);
    else if (ggdauTruePart->PDG == 2212) anaUtils::_categ->SetObjectCode("ggdaughterkaon", 2212,       CATOTHER, indx2, indx3);
    else                                 anaUtils::_categ->SetObjectCode("ggdaughterkaon", CATOTHER,   CATOTHER, indx2, indx3);
  }
}

//********************************************************************
std::pair<Int_t, Int_t> kaonAnaUtils::GetKaonDecayMode(AnaTrueParticleB* trueKaon, std::vector<AnaTrueParticleB*> trueParts){
//********************************************************************

  bool muon = false;
  bool pion = false;
  int ID    = -999;
  int dmode = 0;

  //loop over daughters
  for(int i = 0; i < (int)static_cast<AnaTrueParticle*>(trueKaon)->Daughters.size(); i++){
    for(int j = 0; j < (int)trueParts.size(); j++){
      if(trueParts[j]->ID == static_cast<AnaTrueParticle*>(trueKaon)->Daughters[i]){
	if(abs(trueParts[j]->PDG) == 13){
	  muon = true;
	  ID = trueParts[j]->ID;
	}
	else if(abs(trueParts[j]->PDG) == 211 || abs(trueParts[j]->PDG) == 111)pion = true;
      }
    }
  }

  if(muon)      dmode = 1;
  else if(pion) dmode = 2;
  else          dmode = 3;

  std::pair<Int_t, Int_t> result = std::make_pair(dmode,ID);

  return result;
}

//********************************************************************
std::pair<bool, Int_t> kaonAnaUtils::MuonFromKaonChain(AnaTrueParticleB* trueKaon, std::vector<AnaTrueParticleB*> trueParts){
//********************************************************************

  bool muon = false;
  Int_t ID  = -999;
  std::pair<bool, Int_t> dummy;

  //loop over daughters
  for(int i = 0; i < (int)static_cast<AnaTrueParticle*>(trueKaon)->Daughters.size(); i++){
    for(int j = 0; j < (int)trueParts.size(); j++){
      if(trueParts[j]->ID == static_cast<AnaTrueParticle*>(trueKaon)->Daughters[i]){
	if(trueParts[j]->PDG == 321){
	  dummy = MuonFromKaonChain(trueParts[j],trueParts);
	  muon  = dummy.first;
	  ID    = dummy.second;
	}
	if(trueParts[j]->PDG == -13){
	  muon = true; 
	  ID = trueParts[j]->ID;
	  break;
	}
      }
      if(muon)break;
    }
    if(muon)break;
  } 
  
  std::pair<bool, Int_t> result = std::make_pair(muon,ID);
  
  return result;
}
