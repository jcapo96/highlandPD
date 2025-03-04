#include "QPIXUtils.hxx"

#include "QPIXDataClasses.hxx"
#include "CategoryManager.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void QPIXUtils::AddCustomCategories(){
//********************************************************************

  AddSingleRadiogenicCategory();
  AddChainRadiogenicCategory();
}

//********************************************************************
void QPIXUtils::AddSingleRadiogenicCategory(){
//********************************************************************

  std::string types[] = {"B8", "hep", "del #gamma", 
			 "Ar39", 
			 "Ar42", "K42",
			 "Kr85",
			 "Rn222", "Po218", "Pb214", "At218", "Rn218", "Bi214", "Po214", "Tl210", "Pb209", "Bi209", "Pb210", "Hg206", "Bi210", "Tl206", "Po210",
			 "cavern #gamma", 
			 NAMEOTHER};
  int codes[]         = {8, 9, 2200, 
			 1000180390, 
			 1000180420, 1000190420, 
			 1000360850,
			 1000862220, 1000842180, 1000822140, 1000852180, 1000862180, 1000832140, 1000842140, 1000812100, 1000822090, 1000832090, 1000822100, 1000802060, 1000832100, 1000812060, 1000842100, 
			 22, 
			 CATOTHER};
  int colors[]        = {2, 208, 224, 
			 8,
			 9, 216,
			 218,
			 13, 16, 28, 37, 48, 111, 121, 133, 142, 145, 148, 177, 105, 205, 158,  
			 70, 
			 COLOTHER};
  const int NPART = sizeof(types)/sizeof(types[0]);
  
  std::reverse(types,  types  + NPART);
  std::reverse(codes,  codes  + NPART);
  std::reverse(colors, colors + NPART);

  anaUtils::_categ->AddCategory("qpixsingle", NPART, types, codes, colors);
}

//********************************************************************
void QPIXUtils::AddChainRadiogenicCategory(){
//********************************************************************

  std::string types[] = {"B8", "hep", "del #gamma", "Ar39 chain", "Ar42 chain", "Kr85 chain", "Rn222 chain", "cavern #gamma", "neutrons" ,NAMEOTHER};
  int codes[]         = {8   , 9    , 2200        , 39          , 42          , 85          , 222          , 22             , 2112       ,CATOTHER};
  int colors[]        = {2   , 208  , 224         , 8           , 9           , 90          , 226          , 70             , 38         ,COLOTHER};
  const int NPART = sizeof(types)/sizeof(types[0]);
  
  std::reverse(types,  types  + NPART);
  std::reverse(codes,  codes  + NPART);
  std::reverse(colors, colors + NPART);

  anaUtils::_categ->AddCategory("qpixchain", NPART, types, codes, colors);
}


//********************************************************************
void QPIXUtils::FillCustomCategories(AnaEventPD* event){
//********************************************************************
  
  FillSingleRadiogenicCategory(event);
  FillChainRadiogenicCategory(event);
}

//********************************************************************
void QPIXUtils::FillSingleRadiogenicCategory(AnaEventPD* event){
//********************************************************************

  AnaEventInfoQPIX* EventInfo = static_cast<AnaEventInfoQPIX*>(event->EventInfo);
  AnaTrueParticleB* TruePart = event->TrueParticles[0];
  if(!EventInfo->IsBackground){
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::B8)
      anaUtils::_categ->SetCode("qpixsingle", 8, CATOTHER);
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::hep)
      anaUtils::_categ->SetCode("qpixsingle", 9, CATOTHER);
    if(TruePart->PDG == 22)
      anaUtils::_categ->SetCode("qpixsingle", 2200, CATOTHER);
  }
  else
    anaUtils::_categ->SetCode("qpixsingle", TruePart->PDG, CATOTHER);
}

//********************************************************************
void QPIXUtils::FillChainRadiogenicCategory(AnaEventPD* event){
//********************************************************************

  AnaEventInfoQPIX* EventInfo = static_cast<AnaEventInfoQPIX*>(event->EventInfo);
  AnaTrueParticleB* TruePart = event->TrueParticles[0];
  if(!EventInfo->IsBackground){
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::B8)
      anaUtils::_categ->SetCode("qpixchain", 8, CATOTHER);
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::hep)
      anaUtils::_categ->SetCode("qpixchain", 9, CATOTHER);
    if(TruePart->PDG == 22)
      anaUtils::_categ->SetCode("qpixchain", 2200, CATOTHER);
  }
  else{
    if(TruePart->PDG == 1000180390)
      anaUtils::_categ->SetCode("qpixchain", 39, CATOTHER);
    else if(TruePart->PDG == 1000180420 || TruePart->PDG == 1000190420)
      anaUtils::_categ->SetCode("qpixchain", 42, CATOTHER);
    else if(TruePart->PDG == 1000360850)
      anaUtils::_categ->SetCode("qpixchain", 85, CATOTHER);
    else
      anaUtils::_categ->SetCode("qpixchain", 222, CATOTHER);
  }
  //TOFIX
  if(TruePart->PDG == 2112)anaUtils::_categ->SetCode("qpixchain", 2112, CATOTHER);
}

//********************************************************************
void QPIXUtils::AddQPIXVariables(OutputManager& output){
//********************************************************************
  
  AddVarI(    output, is_background,    "is this a background event?");
  AddVarI(    output, has_potasium ,    "has this event potasium?");
  AddVarI(    output, has_gamma    ,    "has this event a delayed gamma?");
  AddVarF(    output, nu_E         ,    "neutrino energy");
  AddVar4VF(  output, nu_pos       ,    "neutrino position");
  AddVarF(    output, e_E          ,    "electron energy");
  AddVar4VF(  output, gamma_pos    ,    "gamma position");
  AddVarFixVI(output, nphotons     ,    "detected photons per plane", 7);
  AddVarFixVI(output, wf_total     ,    "PDS total wf"              , 1024);
  AddVarFixVI(output, wf_plane_0   ,    "PDS plane 0 wf"            , 1024);
  AddVarFixVI(output, wf_plane_1   ,    "PDS plane 1 wf"            , 1024);
  AddVarFixVI(output, wf_plane_2   ,    "PDS plane 2 wf"            , 1024);
  AddVarFixVI(output, wf_plane_3   ,    "PDS plane 3 wf"            , 1024);
  AddVarFixVI(output, wf_plane_4   ,    "PDS plane 4 wf"            , 1024);
  AddVarFixVI(output, wf_plane_5   ,    "PDS plane 5 wf"            , 1024);
}

//********************************************************************
void QPIXUtils::FillQPIXVariables(OutputManager& output, AnaEventPD* event){
//********************************************************************

  AnaEventInfoQPIX* EventInfo = static_cast<AnaEventInfoQPIX*>(event->EventInfo);
  
  output.FillVar(               is_background, (int)EventInfo->IsBackground);
  output.FillVar(               has_potasium , (int)EventInfo->HasPotasium);
  output.FillVar(               has_gamma    , (int)EventInfo->HasGamma);
  output.FillVectorVarFromArray(nphotons     , EventInfo->DetectedPhotons, 7);
  for(int i = 0; i < 1024; i++){
    output.FillVectorVar(       wf_total     , EventInfo->PDSwf[0][i]    , i);
    output.FillVectorVar(       wf_plane_0   , EventInfo->PDSwf[1][i]    , i);
    output.FillVectorVar(       wf_plane_1   , EventInfo->PDSwf[2][i]    , i);
    output.FillVectorVar(       wf_plane_2   , EventInfo->PDSwf[3][i]    , i);
    output.FillVectorVar(       wf_plane_3   , EventInfo->PDSwf[4][i]    , i);
    output.FillVectorVar(       wf_plane_4   , EventInfo->PDSwf[5][i]    , i);
    output.FillVectorVar(       wf_plane_5   , EventInfo->PDSwf[6][i]    , i);
  }
  //if(EventInfo->IsBackground)return;
  AnaTrueVertex* vertex = static_cast<AnaTrueVertex*>(event->TrueVertices[0]);
  if(!vertex)return;
  output.FillVar(               nu_E         , vertex->NuEnergy);
  output.FillVectorVarFromArray(nu_pos       , vertex->Position, 4);
  output.FillVar(               e_E          , vertex->LeptonMom);
  
  //gamma information
  AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(event->TrueParticles[0]);
  if(!truepart)return;
  if(truepart->PDG!=22)return;
  output.FillVectorVarFromArray(gamma_pos    , truepart->Position, 4);
}

