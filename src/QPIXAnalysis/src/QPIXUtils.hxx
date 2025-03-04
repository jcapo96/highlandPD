#ifndef QPIXUtils_h
#define QPIXUtils_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdDataClasses.hxx"

namespace QPIXUtils{

  // Add custom categories
  void AddCustomCategories();
  void AddSingleRadiogenicCategory();
  void AddChainRadiogenicCategory();
  // Fill Custom categories
  void FillCustomCategories(AnaEventPD* event);
  void FillSingleRadiogenicCategory(AnaEventPD* event);
  void FillChainRadiogenicCategory(AnaEventPD* event);

  //Add QPIX variables
  void AddQPIXVariables(OutputManager& output);
  //Fill QPIX variables
  void FillQPIXVariables(OutputManager& output, AnaEventPD* event);


  // Enum with unique indexes for output tree variables  
  enum enumQPIXTree{

    // selected track (beam particle) true info
    is_background = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    has_potasium,
    has_gamma,
    nu_E, 
    nu_pos,
    e_E,
    gamma_pos,
    nphotons,
    wf_total,
    wf_plane_0,
    wf_plane_1,
    wf_plane_2,
    wf_plane_3,
    wf_plane_4,
    wf_plane_5
  };
}

#endif
