#ifndef pionMomentumTree_h
#define pionMomentumTree_h

#include "OutputManager.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include <vector>

namespace pionMomentumTree {

  void AddPionMomentumVariables_BeamParticleReco(OutputManager& output);
  void FillPionMomentumVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beampart,
                                                  const std::vector<double>& mainthetascatter,
                                                  const std::vector<double>& mainsegmentlength);

  enum enumPionMomentumMicroTrees {
    mainispandora = standardPDTree::enumStandardMicroTreesLast_standardPDTree + 1,
    maintruepdg,
    maintruestartmomentum,
    maintrueendmomentum,
    maintrueid,
    beamtrueid,
    mainmcsntriplets,
    mainmcsdeltatheta,
    mainmcssegmentlength,
    mainnhitscollection,
    mainntrjpoints,
    mainnvalidrrhits,
    mainnvalidxyzhits,
    mainnvalidxyztrj,
  };

}  // namespace pionMomentumTree

#endif
