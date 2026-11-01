#ifndef pdMomReconstructionFromParams_h
#define pdMomReconstructionFromParams_h

#include "pdMomReconstructionJointK0s.hxx"
#include <string>

namespace pdMomReconstruction {

void FillPionTLEFitConfig_FromPionMomentumParams(PionTLEFitConfig& cfg);
void FillPionMCSConfig_FromPionMomentumParams(PionMCSConfig& cfg);
void FillPionMomentumRecoConfig_FromPionMomentumParams(PionMomentumRecoConfig& cfg);

void FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams(JointK0sTwoPionGridFitConfig& cfg);

/// Sets the dE/dx likelihood model (pdDedxPdfModel) from <prefix>.DedxPdfModel* parameters; the dE/dx window is
/// <prefix>.FreeRangeDedxMin/MaxMeVcm. Missing parameters keep the defaults (model disabled: legacy behaviour).
void SetDedxPdfModel_FromParams(const std::string& prefix);

}  // namespace pdMomReconstruction

#endif
