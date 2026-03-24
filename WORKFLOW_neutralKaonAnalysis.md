# neutralKaonAnalysis workflow (current)

Technical reference for the current implementation in `highlandPD`.

## Scope

- Main entry point and execution flow
- Selection and candidate-building chain
- Parameters and output path
- Current risks/mismatches

## Key files

- `src/neutralKaonAnalysis/app/RunNeutralKaonAnalysis.cxx`: binary entry point (`main`).
- `src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`: analysis lifecycle (`Initialize`, microtree/truth filling, categories).
- `src/neutralKaonAnalysis/src/neutralKaonSelection.cxx`: ordered selection steps and cuts.
- `src/neutralKaonAnalysis/src/neutralKaonTree.cxx`: neutral kaon microtree variable definitions/filling.
- `src/neutralKaonAnalysis/src/neutralKaonTruthTree.cxx`: truth-tree variable handling.
- `src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.cxx`: signal/category logic.
- `src/neutralKaonAnalysis/src/EventDisplayDataTree.cxx`: event-display wrapper bridge.
- `src/neutralKaonAnalysis/src/neutralKaonEventDisplay.cxx`: analysis-specific display content hooks.
- `src/pdUtils/src/pdAnnihilationUtils.cxx`: annihilation vertex creation/scoring.
- `src/pdUtils/src/pdNeutralUtils.cxx`: neutral candidate construction/scoring.
- `src/neutralKaonAnalysis/parameters/neutralKaonAnalysis.parameters.dat`: analysis parameters.
- `CMakeLists.txt`: executable target declaration (`neutralKaon.exe`).

## Startup and execution flow

1. `RunNeutralKaonAnalysis.cxx` creates `neutralKaonAnalysis` and passes it to `AnalysisLoop`.
2. `AnalysisLoop` parses command-line options (`-o`, `-p`, etc.) and loads parameters.
3. `neutralKaonAnalysis::Initialize()`:
   - reads neutral-kaon parameters,
   - sets minimum accumulation level,
   - initializes event display tree wrapper,
   - registers custom categories.
4. Input converters are registered (notably `minitreefiltered` and `PDSPAnalyzerTree` paths).
5. `DefineSelections()` registers `neutralKaonSelection`.
6. Per-event selection steps run in sequence:
   - beam track preselection,
   - particle finding,
   - neutral-candidate creation,
   - K0 quality/geometry cuts.
7. Candidate building uses:
   - `pdAnnihilationUtils::CreateVertices(...)`,
   - `pdCreationUtils::CreateCreationVertices(...)`,
   - `pdNeutralUtils::CreateNeutrals(...)`.
8. Output fill phase:
   - analysis/microtree variables via `neutralKaonTree`,
   - event display payload (`EventDisplayData`),
   - truth/config/header trees through the loop lifecycle.

## Selection chain (neutralKaonSelection)

Current selection path is implemented as a step chain including:

- `FindBeamTrackAction`, `BeamTrackExistsCut`
- `FindAllParticlesAction`, `HasEnoughParticlesCut`
- `FindNeutralCandidatesAction`, `HasNeutralCandidatesCut`
- K0 quality cuts such as:
  - `K0StartEndDirCut`
  - `K0VtxOpeningCut`
  - `K0ParentTruePDGCut`
  - `K0LengthCut`
  - daughter length cuts

## Parameters and configuration flow

- Parameter override file can be passed with `-p` (processed by `AnalysisLoop`).
- Neutral-kaon-specific configuration comes from:
  - `src/neutralKaonAnalysis/parameters/neutralKaonAnalysis.parameters.dat`
- Representative parameter groups currently used:
  - accumulation level and SCE toggles,
  - annihilation/creation vertex geometry thresholds,
  - creation-annihilation separation,
  - filtering and track-fit controls.

## Output model

- Output file is controlled by `AnalysisLoop` (`-o`).
- Main produced content:
  - analysis trees (with neutral kaon variables),
  - `truth`,
  - `config`,
  - `header`,
  - `EventDisplayData`.
- Operational convention from run docs:
  - inputs under `/DATA`,
  - outputs under `/microTrees`.

## Detailed reconstruction workflow (candidate building)

This section describes the reconstruction logic from raw event particles to selected neutral-kaon candidates.

1. **Beam and particle preconditions**
   - `FindBeamTrackAction` + `BeamTrackExistsCut` ensure a beam track in TPC is available.
   - `FindAllParticlesAction` counts particles with valid start positions.
   - `HasEnoughParticlesCut` requires at least two valid particles.

2. **Annihilation vertex construction**
   - `pdAnnihilationUtils::CreateVertices(...)` scans particle pairs and applies:
     - start/end position validity checks,
     - same-parent prefilter for the two daughter tracks,
     - PID rejection against proton/kaon-like tracks (`Chi2PID` thresholds),
     - geometric proximity using fitted line closest approach.
   - Vertex position is computed with one of three methods (`VertexFindingMethod`):
     - geometric,
     - TMinuit fit,
     - Kalman fit.
   - Vertex validation includes score sanity and max closest-approach requirement.
   - Vertices are filtered so one daughter track is not reused in multiple annihilation vertices.

3. **Creation vertex construction**
   - `pdCreationUtils::CreateCreationVertices(...)` builds creation vertices around beam/secondary topology.
   - Current call in selection uses creation radius parameter and an empty exclusion list.
   - Additional filtering consistency is applied later when combining with annihilation vertices.

4. **Neutral candidate generation (`AnaNeutralParticlePD`)**
   - `pdNeutralUtils::CreateNeutrals(...)` forms combinations of creation and annihilation vertices.
   - Pair-level guards:
     - rejects combinations reusing annihilation-vertex particles at creation side,
     - enforces minimum creation-annihilation separation,
     - rejects beam-only creation vertices (`DistanceScore < 0.1`).
   - For each surviving pair:
     - assigns creation/annihilation vertices and parent,
     - computes start/end position, directions, length,
     - computes impact parameter via parent-track extrapolation,
     - computes neutral score metrics and true-equivalent associations,
     - computes invariant mass from annihilation daughters (pion-mass hypothesis).

5. **Direction and score filtering**
   - Candidates are first filtered to positive alignment along beam direction.
   - Scores are normalized using alignment relative to beam direction.
   - Final pruning keeps only one neutral candidate per annihilation vertex (best score).

6. **Selection-level candidate cuts**
   - Event passes if at least one candidate satisfies each event-level cut:
     - start/end direction consistency,
     - annihilation opening-angle threshold,
     - true parent PDG requirement,
     - K0 and daughter length thresholds.

7. **Tree filling**
   - For events that pass selection, all candidates currently in `ToyBoxNeutralKaon` are written to microtree variables.
   - `EventDisplayData` is also filled per selected event.

## Current branch/commit (template)

- Branch: `<fill>`
- Commit: `<fill>`
- Date checked: `<fill>`

## Known gaps/uncertainties

- Executable generated by compilation is `neutralKaon.exe` (from `highlandPD/CMakeLists.txt`).
- Event display in `highland` is generic base infrastructure; analysis-specific behavior is implemented in `highlandPD`.
- Some documented run commands are placeholders; manifest entries should capture full command lines per production output.

## Potential errors or bugs to check

1. **Truth-fill coverage inconsistency**
   - `CheckFillTruthTreePD` currently accepts only `PDG == 310`, while `FillTruthTree` contains logic for `310`, `130`, and `311`.
   - Risk: `K0L`/`K0` branches in `FillTruthTree` may never execute in practice.

2. **Event-level cuts vs candidate-level output**
   - Selection cuts are event-level (“at least one candidate passes”), but microtree filling writes all candidates in the box.
   - Risk: output includes candidates that would fail intended per-candidate quality requirements.

3. **Hardcoded debug print for one event**
   - `FindNeutralCandidatesAction` includes explicit diagnostics for run/subrun/event `22591250/40/4959`.
   - Risk: unnecessary stdout spam and non-general behavior in production processing.

4. **Potential null-pointer assumptions in candidate creation**
   - Candidate filling uses creation-vertex beam members (e.g., timestamp field) with implicit availability.
   - Risk: crash if a malformed creation vertex lacks expected beam pointer state.

5. **Score meaning can be overwritten**
   - `CalculateNeutralScore` computes detailed multi-component score, but later normalization rewrites `NeutralScore` from beam-alignment ranking.
   - Risk: ambiguity between “physics score” and “final ranking score”, making tuning/debugging confusing.

## Next verification actions

1. Confirm whether truth tree should include only `PDG=310` or also `130/311`, then align `CheckFillTruthTreePD` and documentation.
2. Decide whether candidate cuts should be applied per candidate before writing microtrees.
3. Guard or remove hardcoded event diagnostics for production runs.
4. Record one full production command in `RUN_MANIFEST.md` with concrete file names and parameter hash.
5. Validate that one produced file in `/microTrees` can be fully reproduced from manifest fields.
