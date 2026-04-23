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

## Main/beam variable naming convention

- For variables from `box().MainTrack`, use the `main` prefix.
- For variables from `beam->BeamParticle`, use the `beam` prefix.
- Do not use `_` in variable names.
- Do not use capital letters in variable names.
- Use straightforward lowercase names (examples: `mainispandora`, `beamtrueid`).

## Notes

- This file stores stable analysis goals and should be updated only when analysis objectives change.
- Implementation details, workflows, and run procedures should be documented separately.
