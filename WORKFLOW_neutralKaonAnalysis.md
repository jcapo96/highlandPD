# neutralKaon framework workflow (detailed chain)

Detailed reconstruction reference for neutral-particle building in `highlandPD`, including where each part of the chain lives.

## 1) Entry points and framework plumbing

- **Binary entrypoint**: `highlandPD/src/neutralKaonAnalysis/app/RunNeutralKaonAnalysis.cxx`
  - Creates `neutralKaonAnalysis` and runs it through `AnalysisLoop`.
- **Build target**: `highlandPD/CMakeLists.txt`
  - Produces `neutralKaon.exe`.
- **Core loop (framework side)**: `highland/src/highland2/highlandCore/src/AnalysisLoop.cxx`
  - Handles options (`-o`, `-n`, `-s`, `-p`), event loop, tree lifecycle.
- **Analysis class**: `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`
  - `Initialize()`, `DefineSelections()`, `FillMicroTrees()`, `FillToyVarsInMicroTrees()`, truth filling.

## 2) Runtime configuration and parameters

- **Default parameter file**:
  - `highlandPD/src/neutralKaonAnalysis/parameters/neutralKaonAnalysis.parameters.dat`
- **Override mechanism**:
  - `neutralKaon.exe -p <override.parameters.dat>`
- **Typical parameter groups used by reconstruction**:
  - vertex radii and fit controls (annihilation/creation),
  - candidate pairing constraints (creation-annihilation separation),
  - filtering/selection thresholds,
  - debug and event-display toggles.

## 3) Selection chain where reconstruction is triggered

- **Selection definition**: `highlandPD/src/neutralKaonAnalysis/src/neutralKaonSelection.cxx`
  - `neutralKaonSelection::DefineSteps()`
  - Ordered chain:
    1. beam track actions/cuts,
    2. particle preconditions,
    3. `FindNeutralCandidatesAction`,
    4. event-level K0 quality cuts.
- **Candidate builder action**:
  - `FindNeutralCandidatesAction::Apply(...)` in the same file.
  - This is the main orchestration point for neutral reconstruction.

## 4) Reconstruction chain (neutral particles) with exact code locations

### 4.1 Input objects and preconditions

- **File**: `highlandPD/src/neutralKaonAnalysis/src/neutralKaonSelection.cxx`
- **Function**: `FindNeutralCandidatesAction::Apply(...)`
- Reads event particles, beam context, and parameter values (e.g., daughter/creation radii).

### 4.2 Build annihilation vertices (V-like decay side)

- **Orchestrator call** (from selection):
  - `pdAnnihilationUtils::CreateVertices(event, maxDaughterDistance)`
- **Implementation file**:
  - `highlandPD/src/pdUtils/src/pdAnnihilationUtils.cxx`
- **Main functions**:
  - `CreateVertices(...)`
  - `CreateVerticesCommon(...)`
  - vertex-position strategies:
    - `FindVertexPositionGeometric(...)`
    - `FindVertexPositionWithFit(...)`
    - `FindVertexPositionKalman(...)`
  - fit utility:
    - `FindVertexPositionFit(...)`
- **What happens conceptually**:
  - particle-pair scan,
  - quality and topology constraints,
  - reconstructed annihilation vertex position/direction/fit products,
  - internal filtering and scoring.

### 4.3 Build creation vertices (beam + partner side)

- **Orchestrator call** (from selection):
  - `pdCreationUtils::CreateCreationVertices(event, creationRadius, excludeIDs)`
- **Implementation file**:
  - `highlandPD/src/pdUtils/src/pdCreationUtils.cxx`
- **Main functions**:
  - `CreateCreationVertices(...)`
  - `CalculateCreationVertexScores(...)`
  - `FilterCreationVerticesByScore(...)`
- **Helper math (line closest approach / Pandora vertex helper)**:
  - `highlandPD/src/pdUtils/src/pdNeutralHelpers.cxx`
  - `CalculatePandoraVertexPosition(...)`
  - `CalculateLineMinDistanceMidpoint(...)`

### 4.4 Combine creation + annihilation into neutral candidates

- **Orchestrator call** (from selection):
  - `pdNeutralUtils::CreateNeutrals(event, creationVertices, annihilationVertices)`
- **Implementation file**:
  - `highlandPD/src/pdUtils/src/pdNeutralUtils.cxx`
- **Main functions**:
  - `CreateNeutrals(...)`
  - `CalculateNeutralScore(...)`
  - `NormalizeNeutralParticleScores(...)`
  - `FilterNeutralsByScore(...)`
- **What is built per candidate (`AnaNeutralParticlePD`)**:
  - creation/annihilation vertex links,
  - parent association,
  - neutral trajectory (start/end/direction/length),
  - impact parameter and derived kinematics,
  - score components and final ranking,
  - truth-equivalent association for validation studies.

## 5) Event-level neutral-K0 cuts after candidate creation

- **File**: `highlandPD/src/neutralKaonAnalysis/src/neutralKaonSelection.cxx`
- **Representative cuts in chain**:
  - `K0StartEndDirCut`
  - `K0VtxOpeningCut`
  - `K0ParentTruePDGCut`
  - `K0LengthCut`
  - daughter-length cuts
- **Important behavior**:
  - these are applied as selection steps at event level; event may pass if at least one candidate satisfies the condition.

## 6) Output writing: where reconstructed content goes

### 6.1 Microtree (main analysis output)

- **Variables declaration/fill**:
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
- **Called from**:
  - `neutralKaonAnalysis::FillMicroTrees()` in
    `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`

### 6.2 Truth-validation tree

- **Implementation**:
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTruthTree.cxx`
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTruthTree.hxx`
- **Dispatch logic**:
  - `CheckFillTruthTreePD(...)`, `FillTruthTree(...)` in `neutralKaonAnalysis.cxx`

### 6.3 Event display payload

- **Event display tree bridge**:
  - `highlandPD/src/neutralKaonAnalysis/src/EventDisplayDataTree.cxx`
- **Analysis-specific event display class**:
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonEventDisplay.cxx`
  - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonEventDisplay.hxx`
- **Generic display base (framework side)**:
  - `highland/src/highland2/highlandTools/src/EventDisplayBase.cxx/.hxx`
  - `highland/src/highland2/highlandTools/src/DrawingTools.cxx/.hxx`

## 7) Data model classes used by reconstruction

- **Main ProtoDUNE data classes**:
  - `highlandPD/src/pdEventModel/src/pdDataClasses.hxx`
  - `highlandPD/src/pdEventModel/src/pdDataClasses.cxx`
- **Backward-compat neutral include wrapper**:
  - `highlandPD/src/pdEventModel/src/pdNeutralDataClasses.hxx`
  - `highlandPD/src/pdEventModel/src/pdNeutralDataClasses.cxx`

## 8) Minimal chain you can follow to audit one event

1. Start in `RunNeutralKaonAnalysis.cxx` (`main`).
2. Jump to `neutralKaonAnalysis::Initialize()` and `DefineSelections()` in `neutralKaonAnalysis.cxx`.
3. Open `neutralKaonSelection::DefineSteps()` and `FindNeutralCandidatesAction::Apply(...)`.
4. Follow the three reconstruction calls in order:
   - `pdAnnihilationUtils::CreateVertices(...)`
   - `pdCreationUtils::CreateCreationVertices(...)`
   - `pdNeutralUtils::CreateNeutrals(...)`
5. Follow candidate filters/scores in `pdNeutralUtils.cxx`.
6. Inspect `FillMicroTrees()` in `neutralKaonAnalysis.cxx`.
7. Inspect per-variable filling in `neutralKaonTree.cxx`.
8. If validating visually, inspect event-display fill path:
   - `EventDisplayDataTree.cxx` -> `neutralKaonEventDisplay.cxx`.

## 9) Practical run/check commands

From `/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW`:

```bash
source highland/scripts/setup.sh
source highlandPD/scripts/setup.sh
neutralKaon.exe -n 20 -o /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/microTrees/nk_check.root /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA/6GeV_prod4a_00_minitree_2023-01-27.root
```

Event display listing:

```bash
root -l -b <<'EOF'
gSystem->Load("libhighland");
gSystem->Load("libhighlandPD");
DrawingTools d("/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/microTrees/nk_check.root");
d.ListEvtDisplay();
.q
EOF
```
