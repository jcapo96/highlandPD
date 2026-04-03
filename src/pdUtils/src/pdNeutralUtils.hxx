#ifndef pdNeutralUtils_h
#define pdNeutralUtils_h

#include "pdDataClasses.hxx"
#include "BaseDataClasses.hxx"
#include <vector>

namespace pdNeutralUtils {

  // Create neutral candidates from annihilation vertices.
  std::vector<AnaNeutralParticlePD*> CreateNeutralsFromAnnihilationVertices(
      AnaEventB& event,
      const std::vector<AnaAnnihilationVertexPD*>& annihilationVertices);

}

#endif

