//header
#include "BrokenTracks.hxx"

//highland
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"
#include "dEdxUtils.hxx"

//converters
#include "HighlandMiniTreeConverter.hxx"
#include "MichelRemovingTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//selections
#include "BrokenTracksSelection.hxx"

//corrections
#include "ParticlePositionSCECorrection.hxx"
#include "BeamMomSmearingCorrection.hxx"

//********************************************************************
BrokenTracks::BrokenTracks(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************
  
}

//********************************************************************
bool BrokenTracks::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("BrokenTracks.MinAccumLevelToSave"));
                                                       
  // Define standard categories for color drawing
  anaUtils::AddStandardCategories();        // This is for the candidate particle
  anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle

  AddBrokenCategory();

  return true;
}

//********************************************************************
void BrokenTracks::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("PDSPAnalyzer",   new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void BrokenTracks::DefineSelections(){
//********************************************************************

  sel().AddSelection("Broken tracks Selection", "Broken tracks on APAs", new BrokenTracksSelection(false));  
}

//********************************************************************
void BrokenTracks::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
  corr().AddCorrection(0, "SCE correction"  , new ParticlePositionSCECorrection());
  corr().AddCorrection(1, "Beam smearing correction"  , new BeamMomSmearingCorrection());
}

//********************************************************************
void BrokenTracks::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();
}

//********************************************************************
void BrokenTracks::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();
}

//********************************************************************
void BrokenTracks::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  standardPDTree::AddStandardVariables_EventInfo(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleReco(output());
  standardPDTree::AddStandardVariables_BeamParticleTrue(output());
  standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),5);
  standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),5);

  AddVarI(output(), seltrk_ndau_APAs, "seltrk ndaughters in APA region");
}

//********************************************************************
void BrokenTracks::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();
}

//********************************************************************
void BrokenTracks::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if(addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  standardPDTree::FillStandardVariables_EventInfo(        output(), static_cast<AnaEventInfoPD*>(GetEvent().EventInfo));
  standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);
  standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack);
  standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);    
  
  int ndau_APAs = 0;
  int ntracks = std::min((int)box().MainTrack->DaughtersIDs.size(),5);
  for(int idau = 0; idau < ntracks; idau++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(box().MainTrack->Daughters[idau]);
    standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), dau);
    standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), dau);
    output().IncrementCounter(standardPDTree::seltrk_ndau);
    double Z = dau->PositionEnd[2];
    if(Z>205 && Z < 249)ndau_APAs++;
  }

  output().FillVar(seltrk_ndau_APAs, ndau_APAs);
}

//********************************************************************
void BrokenTracks::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  // if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);

}

//********************************************************************
void BrokenTracks::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  // baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void BrokenTracks::FillCategories(){
//********************************************************************

  if(box().MainTrack){
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils
    FillBrokenCategory(box().MainTrack);
  }

  // For the beam 
  AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
  if (beamPart){
    anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
  }

}

//********************************************************************
void BrokenTracks::AddBrokenCategory(){
//********************************************************************

  std::string part_types[] = {"Broken", "Not broken", "No Truth", NAMEOTHER};
  int part_codes[]         = {1       , 2           , 3, CATOTHER};
  int part_colors[]        = {4       , 2           , COLNOTRUTH, COLOTHER};
  const int NPART = sizeof(part_types)/sizeof(part_types[0]);
  
  std::reverse(part_types,  part_types  + NPART);
  std::reverse(part_codes,  part_codes  + NPART);
  std::reverse(part_colors, part_colors + NPART);

  anaUtils::_categ->AddCategory("broken", NPART, part_types, part_codes, part_colors);
  // anaUtils::_categ->AddObjectCategory("brokendaughter", standardPDTree::seltrk_ndau, "seltrk_ndau", 
  // 				      NPART, part_types, part_codes, part_colors, 
  // 				      1, -100);
}

//********************************************************************
void BrokenTracks::FillBrokenCategory(AnaParticlePD* part){
//********************************************************************

  if(!part)return;
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return;

  // bool broken = false;
  // int ndau = part->Daughters.size();
  // int broken_dau = -1;
  // for(int idau = 0; idau < ndau; idau++){
  //   AnaTrueParticlePD* dauTruePart = static_cast<AnaTrueParticlePD*>(part->Daughters[idau]->TrueObject);
  //   if(!dauTruePart)continue;
  //   if(truePart->ID == dauTruePart->ID){
  //     broken = true;
  //     broken_dau = idau;
  //     break;
  //   }
  // }

  // //beam part category
  // if(broken)anaUtils::_categ->SetCode("broken", 1, CATOTHER);
  // else      anaUtils::_categ->SetCode("broken", 2, CATOTHER);

  // //daughters category
  // for(int idau = 0; idau < ndau; idau++){
  //   if(idau == broken_dau)anaUtils::_categ->SetObjectCode("brokendaughter", 1, CATOTHER, -1);
  //   else                  anaUtils::_categ->SetObjectCode("brokendaughter", 2, CATOTHER, -1);
  // }

  if(!truePart->Matched)anaUtils::_categ->SetCode("broken", 3, CATOTHER);
  else{
    double Z = truePart->PositionEnd[2];
    if(Z>215.0 && Z<234.0)anaUtils::_categ->SetCode("broken", 2, CATOTHER);
    else if (Z>234)       anaUtils::_categ->SetCode("broken", 1, CATOTHER);
    else                  anaUtils::_categ->SetCode("broken", 4, CATOTHER);
  }
}
