#include "pionAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void pionAnaUtils::AddCustomCategories(){
//********************************************************************

  // ------  WORK IN PROGRESS: a category similar to the one by Jake/Francesca -----------


  // --------- daupionana category ----------------
  
  std::string part_types[] = {"Cosmic", "Nuc"     , "Self", "Gamma"  , "pi0Gamma", "Other", "GDaughter", "GGDaughter" , "Michel", "Pi" , "Proton" , "mu"};
  int part_codes[]         = {1       , 2         , 3     , 4        ,  5        ,  6     , 7          , 8            ,  9      , 10   , 11      , 12   };
  int part_colors[]        = {kRed+2  , kOrange+10, kBlue , kOrange+1, kSpring-8 ,  kBlack, kMagenta   , kOrange      ,  kPink  , kTeal, kViolet-3, kRed};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddObjectCategory("daupionana", standardPDTree::seltrk_ndau, "seltrk_ndau", -100, 
                                      NPART, part_types, part_codes, part_colors);


  /*
  struct categ_type_def{
    categ_def(const std::string& name, Int_t code, Int_t color){_name=name;_code=code;_color=color;}
    std::string _name;
    Int_t _code;
    Int_t _color;
  };

  categ_def pion_categ[] = {categ_type_def("Primary #pi^{+}"         ,  1, 1 ),                                                           
                            categ_type_def("Primary #mu"             ,  2, 2 ),                                                               
                            categ_type_def("Cosmic"                  ,  3, 3 ),                                                                    
                            categ_type_def("Extra Primary"           ,  4, 4 ),                                                             
                            categ_type_def("#pi^{+} From Interaction",  5, 5 ),                                                  
                            categ_type_def("#pi^{-} From Interaction",  6, 6 ),                                                  
                            categ_type_def("p From Interaction"      ,  7, 7 ),                                                        
                            categ_type_def("K From Interaction"      ,  8, 8 ),                                                        
                            categ_type_def("Nuc From Interaction"    ,  9, 9 ),                                                      
                            categ_type_def("Decay Product"           , 10,10 ),                                                             
                            categ_type_def("Other"                   , 11,11 )};
  */  

  // --------- pionana2 category ----------------

  std::string pionana2_types[] = {"Primary #pi^{+}",
                                  "Primary #mu",
                                  "Cosmic",
                                  "Extra Primary",
                                  "#pi^{+} From Interaction",
                                  "#pi^{-} From Interaction",
                                  "p From Interaction",
                                  "K From Interaction",
                                  "Nuc From Interaction",
                                  "Decay Product",
                                  "Other"};                                                                    

  
  int pionana2_codes[]         = {1      , 2        , 3         , 4     , 5         , 6      , 7     , 8        ,  9          , 10       , 11      };
  int pionana2_colors[]        = {kBlue-4, kOrange+1, kSpring-8 , kRed+2, kViolet-3 , kCyan-2, kTeal , kYellow-6, kOrange+10  , kOrange-7, kMagenta-10};
  const int NPIONANA2 = sizeof(pionana2_types)/sizeof(pionana2_types[0]);

  std::reverse(pionana2_types,  pionana2_types  + NPIONANA2);
  std::reverse(pionana2_codes,  pionana2_codes  + NPIONANA2);
  std::reverse(pionana2_colors, pionana2_colors + NPIONANA2);
  
  anaUtils::_categ->AddCategory("pionana2",     NPIONANA2, pionana2_types,  pionana2_codes, pionana2_colors);

  
  // --------- pionana3 category ----------------
  
  std::string pionana3_types[] = {"Background","True Abs Signal", "True Cex Signal", "True n-#pi^{0} Signal"};
  int pionana3_codes[]         = {1           , 2               , 3                , 4};                      
  int pionana3_colors[]        = {kBlue       , kGreen          , kMagenta         , kTeal};
  const int NPIONANA3 = sizeof(pionana3_types)/sizeof(pionana3_types[0]);

  std::reverse(pionana3_types,  pionana3_types  + NPIONANA3);
  std::reverse(pionana3_codes,  pionana3_codes  + NPIONANA3);
  std::reverse(pionana3_colors, pionana3_colors + NPIONANA3);
  
  anaUtils::_categ->AddCategory("pionana3",     NPIONANA3, pionana3_types,  pionana3_codes, pionana3_colors);

  
}

//********************************************************************
void pionAnaUtils::FillCustomCategories(AnaEventB* event, AnaParticle* part, PDCounters& counters){
//********************************************************************

  AnaBeamPD* beam  = static_cast<AnaBeamPD*>(event->Beam); 
  
  if (!part->TrueObject) return;  
  AnaTrueParticlePD*  truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);  
  if (!truePart) return;


  // --------- daupionana category ------------------------
  
  AnaParticleB* beamPart = beam->BeamParticle;
  AnaTrueParticlePD* trueBeamPart = static_cast<AnaTrueParticlePD*>(beamPart->TrueObject);  

  for (Int_t i=0;i<std::min((Int_t)100,(Int_t)part->Daughters.size()); ++i){
    pionAnaUtils::FillDaupionanaCategory(event, static_cast<AnaTrueParticlePD*>(part->Daughters[i]->TrueObject), trueBeamPart);
  }

  
  // --------- pionana2 category ------------------------

  if (truePart->ProcessStart ==  AnaTrueParticle::primary && truePart->Origin ==  4){
    if      (truePart->Matched && truePart->PDG==211)   anaUtils::_categ->SetCode("pionana2", 1);  // Primary pion
    else if (truePart->Matched && truePart->PDG==-13)   anaUtils::_categ->SetCode("pionana2", 2);  // Primary muon
    else if (!truePart->Matched)                        anaUtils::_categ->SetCode("pionana2", 4);  // Extra primary
  }
  else if (truePart->Origin == 2)                       anaUtils::_categ->SetCode("pionana2", 3);  // Cosmic
  else if (truePart->Origin == 4 && 
           (truePart->ProcessStart ==  AnaTrueParticle::neutronInelastic ||
            truePart->ProcessStart ==  AnaTrueParticle::protonInelastic ||
            truePart->ProcessStart ==  AnaTrueParticle::piminusInelastic ||
            truePart->ProcessStart ==  AnaTrueParticle::piplusInelastic ||
            truePart->ProcessStart ==  AnaTrueParticle::hadElastic)){
    if      (truePart->PDG==211)                        anaUtils::_categ->SetCode("pionana2", 5);  // #pi^{+} From Interaction
    else if (truePart->PDG==-211)                       anaUtils::_categ->SetCode("pionana2", 6);  // #pi^{-} From Interaction
    else if (truePart->PDG==2212)                       anaUtils::_categ->SetCode("pionana2", 7);  // p From Interaction
    else if (abs(truePart->PDG)==321)                   anaUtils::_categ->SetCode("pionana2", 8);  // k From Interaction
    else if (truePart->PDG>2212)                        anaUtils::_categ->SetCode("pionana2", 9);  // Nuc From Interaction
  }
  else if (truePart->ProcessStart ==  AnaTrueParticle::Decay && truePart->Origin ==  4)  anaUtils::_categ->SetCode("pionana2", 10);  // Decay
  else                                                  anaUtils::_categ->SetCode("pionana2", 11);  // Other
   
  // --------- pionana3 category ------------------------  

  if (trueBeamPart->ProcessEnd == AnaTrueParticle::piplusInelastic && counters.ntrue_beamdaughter_piplus+counters.ntrue_beamdaughter_piminus==0){
    if      (counters.ntrue_beamdaughter_pi0==0)    anaUtils::_categ->SetCode("pionana3", 2);  // Abs
    else if (counters.ntrue_beamdaughter_pi0==1)    anaUtils::_categ->SetCode("pionana3", 3);  // Cex
    else                                            anaUtils::_categ->SetCode("pionana3", 4);  // Npi0
  } 
  else                                              anaUtils::_categ->SetCode("pionana3", 1);  // Background

}

//********************************************************************
void pionAnaUtils::FillDaupionanaCategory(AnaEventB* event, AnaTrueParticlePD* truePart, AnaTrueParticlePD*  trueBeamPart){
//********************************************************************

  if (!trueBeamPart || !truePart) return; 

  bool true_daughter = (trueBeamPart->ID == truePart->ParentID);
  bool is_self       = (trueBeamPart->ID == truePart->ID);

  bool true_grand_daughter = false;
  bool true_great_grand_daughter = false;
  for(UInt_t i = 0; i < trueBeamPart->Daughters.size(); i++){
    AnaTrueParticlePD* trueBeamDaughter = pdAnaUtils::GetTrueParticle(event, trueBeamPart->Daughters[i]);
    if (!trueBeamDaughter) continue;
    for(UInt_t j = 0; j < trueBeamDaughter->Daughters.size(); j++){
      AnaTrueParticlePD* trueBeamGrandDaughter = pdAnaUtils::GetTrueParticle(event, trueBeamDaughter->Daughters[j]);
      if (trueBeamGrandDaughter->ID == truePart->ID){
        true_grand_daughter=true;
        break;
      }
      else if (trueBeamGrandDaughter->ID == truePart->ParentID){
        true_great_grand_daughter=true;
        break;
      }

      /*      else{ 
        for(UInt_t k = 0; k < trueBeamGrandDaughter->Daughters.size(); k++){
          AnaTrueParticlePD* trueBeamGreatGrandDaughter = pdAnaUtils::GetTrueParticle(spill->TrueParticles, trueBeamGrandDaughter->Daughters[k]);
          if (trueBeamGreatGrandDaughter->ID == truePart->ID){
            true_great_grand_daughter=true;
            break;
          }
        }
      }
      */
    }
    if (true_grand_daughter) break;
  }
  
  bool pi0gamma=false;
  if (truePart->PDG==22){
    for(UInt_t i = 0; i < trueBeamPart->Daughters.size(); i++){
      AnaTrueParticlePD* trueBeamDaughter = pdAnaUtils::GetTrueParticle(event, trueBeamPart->Daughters[i]);
      if (!trueBeamDaughter) continue;
      if(trueBeamDaughter->PDG==111){
        std::vector<Int_t> miau = trueBeamPart->Pi0_decay_ID;
        auto pos = std::find(miau.begin(),miau.end(),truePart->ID);
        if(pos != miau.end()){
          pi0gamma=true;                                    
          break;
        }
      }
    }
  }
  
  if      (truePart->Origin  == 2)                                   anaUtils::_categ->SetObjectCode("daupionana", 1);
  else if (is_self)                                                  anaUtils::_categ->SetObjectCode("daupionana", 3);
  else if (pi0gamma)                                                 anaUtils::_categ->SetObjectCode("daupionana", 5);
  else if (true_grand_daughter)                                      anaUtils::_categ->SetObjectCode("daupionana", 7);
  else if (true_great_grand_daughter && abs(truePart->PDG)==11)      anaUtils::_categ->SetObjectCode("daupionana", 9);
  else if (true_great_grand_daughter)                                anaUtils::_categ->SetObjectCode("daupionana", 8);
  else if (true_daughter){
    if      (abs(truePart->PDG)==211)                                anaUtils::_categ->SetObjectCode("daupionana", 10);
    else if (abs(truePart->PDG)==2212)                               anaUtils::_categ->SetObjectCode("daupionana", 11);
    else if (abs(truePart->PDG)==13)                                 anaUtils::_categ->SetObjectCode("daupionana", 12);
    else if (abs(truePart->PDG)>2212)                                anaUtils::_categ->SetObjectCode("daupionana", 2);
    else if (abs(truePart->PDG)==22)                                 anaUtils::_categ->SetObjectCode("daupionana", 4);
    else                                                             anaUtils::_categ->SetObjectCode("daupionana", 6);  }

  else                                                               anaUtils::_categ->SetObjectCode("daupionana", 6);
}
