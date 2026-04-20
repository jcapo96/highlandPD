#ifndef pdUtilsNeutralKaonSignalDiagnostics_h
#define pdUtilsNeutralKaonSignalDiagnostics_h

#include "neutralKaonTree.hxx"

namespace neutralKaonTreeDiagnostics {
  void MaybeAccumulateSignalPionDedxMultiGraphs(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                                size_t candidateIndex);
  void MaybeAccumulateJointMomentumLogLikelihoodGraphs(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                                       size_t candidateIndex);
  void MaybeAccumulateJointObjective2DHeatmaps(AnaNeutralParticlePD* candidate, const AnaEventB& event,
                                               size_t candidateIndex);
  void ResetSignalK0DiagnosticsIfNewEvent(const AnaEventB& event);
  bool RegisterSignalTrueK0Id(Int_t trueK0Id);
  void WriteSignalPionDedxDiagnostics(OutputManager& output);
}

#endif
