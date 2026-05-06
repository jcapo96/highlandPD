#ifndef pdMomReconstructionFromParams_h
#define pdMomReconstructionFromParams_h

#include "pdMomReconstructionJointK0s.hxx"

namespace pdMomReconstruction {

void FillPionTLEFitConfig_FromPionMomentumParams(PionTLEFitConfig& cfg);
void FillPionMCSConfig_FromPionMomentumParams(PionMCSConfig& cfg);

void FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams(JointK0sTwoPionGridFitConfig& cfg);

}  // namespace pdMomReconstruction

#endif
