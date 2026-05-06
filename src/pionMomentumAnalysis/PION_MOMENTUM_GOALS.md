# Pion Momentum Analysis Goals

Canonical goals for the pion momentum analysis line.

## Primary objective

Evaluate the performance of momentum reconstruction algorithms for pions that come from the beam particle sample.

## Algorithms in scope

- `TLEFit`
- `Calorimetry`
- `MCS`
- `Joint fit`
- Any additional momentum estimators explicitly added for comparison

## Analysis outcomes

1. Quantify reconstruction performance for each algorithm (resolution, bias, stability, and relevant diagnostics).
2. Extract calibration/correction parameters from the beam-pion control sample.
3. Produce validated corrections/parameterizations that can be propagated to `neutralKaonAnalysis`.

## Downstream purpose

Use the extracted corrections to improve momentum reconstruction of pions associated with the annihilation vertex in `neutralKaonAnalysis`.

## Parity contract (pionMomentum vs neutralKaon joint momentum)

Implementation keeps one TLE stack (`EstimatePionMomentumTLEFit` → `BuildPionFreeRangeLogLikelihoodVsMomentumCurve` on the track) and one neutral-kaon joint recipe (`FitJointK0sTwoPionMomentaOnGrid` in `pdMomReconstructionJointK0s.cxx`), with annihilation and diagnostics calling the same parameter loaders.

| Role | Code |
|------|------|
| Decode `pionMomentumAnalysis.*` into TLE/MCS structs | `pdMomReconstruction::FillPionTLEFitConfig_FromPionMomentumParams`, `FillPionMCSConfig_FromPionMomentumParams` (`pdMomReconstructionFromParams.cxx`) |
| Decode `neutralKaonAnalysis.FreeRange*`, `JointK0s*`, MCS keys into joint config | `FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams` (same file) |
| Per-daughter diagnostic TLE+MCS curves (no third copy of FreeRange/MCS parsing) | `BuildNeutralKaonJointDiagnosticsCurvesForDaughter` |
| Pion MCS likelihood surface (segmentation, `ComputePionMcsNLLFromScatteringSamples`, `BuildPionMCSLogLikelihoodVsMomentumCurve`, …) | `pdMomReconstruction` in `pdMomReconstructionMCS.hxx` / `.cxx` |
| Two-body grid + mass term on pre-built log L curves | `pdJointK0sPionMomentum::FitJointMomentaOnGrid` and helpers in `pdJointK0sPionMomentum.cxx` |

**Parameter alignment:** comparable physics between the beam-pion line and annihilation daughters requires matching numeric values for the shared quantities (scan limits, interior-hit cuts, dE/dx gates, PDF path, fine scan, low-*p* refine, and joint MCS likelihood knobs). Keys live under `pionMomentumAnalysis.*` vs `neutralKaonAnalysis.FreeRange*` / `JointK0s*` by design; the loaders above are the checklist for what must be kept in sync manually.

**Intentional MCS difference:** standalone MCS (`EstimatePionMomentumMCSReco` / `MCSTrue`) honors `MCSDropFirstNSegments` and `MCSDropLastNSegments`. The joint MCS axis uses `pdMomReconstruction::BuildPionMCSLogLikelihoodVsMomentumCurve`, which does not apply those drops (the neutral-kaon joint loader sets drop counts to zero). Do not expect identical MCS curves between standalone and joint terms unless drops are zero or the joint API is extended.

## Main/beam variable naming convention

- For variables from `box().MainTrack`, use the `main` prefix.
- For variables from `beam->BeamParticle`, use the `beam` prefix.
- Do not use `_` in variable names.
- Do not use capital letters in variable names.
- Use straightforward lowercase names (examples: `mainispandora`, `beamtrueid`).

## Notes

- This file stores stable analysis goals and should be updated only when analysis objectives change.
- Implementation details, workflows, and run procedures should be documented separately.
