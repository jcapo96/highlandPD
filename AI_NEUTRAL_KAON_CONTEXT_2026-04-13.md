# Neutral Kaon Reconstruction Context (AI Handoff)

Last updated: 2026-04-13
Primary branch context: `neutralKaon/creation_vertex/parent`

---

## 1) Purpose of this document

This file is a **current-state handoff** for another AI agent that needs to continue work on neutral kaon reconstruction in `highlandPD` without re-discovering the full codebase.

It captures:
- how neutral candidates are currently reconstructed,
- which parts are reco-only vs truth-assisted,
- what changed recently,
- what analysis goals are driving current development,
- where to edit next.

---

## 2) Current reconstruction workflow (what actually runs now)

### 2.1 Selection step that builds neutral candidates

File: `src/neutralKaonAnalysis/src/neutralKaonSelection.cxx`
Function: `FindNeutralCandidatesAction::Apply(...)`

Current flow:
1. Clear previous candidate containers in `ToyBoxNeutralKaon`.
2. Build annihilation vertices with:
   - `pdAnnihilationUtils::CreateVertices(event, maxDaughterDistance, ...)`
   - `maxDaughterDistance` from `neutralKaonAnalysis.AnnihilationVertexRadius`.
3. Build neutral candidates with:
   - `pdNeutralUtils::CreateNeutralsFromAnnihilationVertices(event, annihilationVertices)`.

Important: current neutral construction is **annihilation-vertex-driven** (one neutral wrapper per annihilation vertex), not a separate creation+annihilation combinator.

---

### 2.2 Annihilation-vertex building logic

File: `src/pdNeutralUtils/src/pdAnnihilationUtils.cxx`
Function: `CreateVerticesCommon(...)`

Core pairing constraints are reco-side:
- both daughters must exist,
- both must have valid start position,
- both must pass minimum collection-plane hit requirement,
- **must share the same reco parent ID** (`daughter1->ParentID == daughter2->ParentID`),
- daughter start-position distance must be within the annihilation radius.

So the “same parent” condition for reconstruction is already on the **reco chain**.

---

### 2.3 Neutral candidate creation from annihilation vertices

File: `src/pdNeutralUtils/src/pdNeutralUtils.cxx`
Function: `CreateNeutralsFromAnnihilationVertices(...)`

Per annihilation vertex:
1. Infer reco parent from daughter `ParentID` lookup in event particles.
2. **Hard guard (recent)**: no neutral candidate is created if:
   - reco parent is null, or
   - reco parent end position is invalid/sentinel (`<= -900` in any component).
3. Build `AnaCreationVertexPD` as beam-parent anchored object.
4. Compute projected parent-end (`pdCreationUtils::ProjectBeamTailOntoStartDirection`) and use it when valid.
5. Set creation coordinates:
   - `CreationVertex->PositionPandora` = selected parent end (currently projected if valid, raw fallback),
   - `CreationVertex->Position` = same selected value.
6. Create `AnaNeutralParticlePD`, link:
   - `CreationVertex`,
   - `AnnihilationVertex`,
   - `Parent`.
7. Fill neutral geometry (`PositionStart`, `PositionEnd`, length/alignment helpers).
8. Attach truth parent only for validation (`neutralParticle->TrueObject = GetCommonTrueParent(...)`).

---

## 3) Signal category behavior (current)

File: `src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.cxx`

### 3.1 Category definition
`AddSignalCandidateCategory()` still defines labels/codes:
- `two_stopping` (1)
- `one_stopping` (5)
- `interacting` (6)
- `legit_vertex_2body` (3)
- `legit_vertex_multibody` (4)
- `background` (2)

### 3.2 Actual fill behavior used for candidates
`FillSignalCandidateCategory(AnaNeutralParticlePD*, ...)` now does:
- `signalCode = 2` by default,
- if `IsSignalCandidate(...)` then use stopping subtype (1/5/6),
- **no longer promotes to 3 or 4 in candidate fill fallback**.

Meaning: the category assignment in analysis output is now effectively **neutral-candidate-based signal vs background/stopping subtype**, matching the requested “no neutral => no signal” policy.

### 3.3 Vertex helper still exists
`GetSignalCategoryCodeForAnnihilationVertex(...)` still has 3/4 fallback logic and is used in annihilation utilities (e.g., momentum-method gating).
That is a helper path and not the direct candidate-category fill.

---

## 4) Tree writing workflow

Files:
- `src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`
- `src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
- `src/neutralKaonAnalysis/src/neutralKaonTree.cxx`

`DefineMicroTrees()` registers blocks in this order:
1. `AddNeutralKaonVariables_K0Particle(...)`
2. `AddNeutralKaonVariables_K0Parent(...)`
3. `AddNeutralKaonVariables_K0CreationVtx(...)`  ← recent
4. `AddNeutralKaonVariables_K0Vtx(...)`

`FillMicroTrees()` loops candidates and calls:
- `FillNeutralKaonVariables(...)`, which dispatches to:
  - particle block,
  - parent block,
  - creation-vtx block,
  - annihilation-vtx block.

---

## 5) New creation-vertex residual variables (recent)

Added in enum and tree add/fill:
- `k0cvtxpandoraresidual`
- `k0cvtxfitresidual`
- `k0cvtxpandoraresidualx/y/z`
- `k0cvtxfitresidualx/y/z`

Definition implemented in `FillNeutralKaonVariables_K0CreationVtx(...)`:
- component residuals are `reco - true`,
- magnitude residual is Euclidean norm.

Current source of reco positions:
- **Pandora residual:** raw reco parent end (`recoParent->PositionEnd`) (uncorrected endpoint).
- **Fit residual:** projected creation endpoint (`creationVtx->Position`) (tail-projected endpoint).

Current source of true creation reference:
1. preferred: `candidate->TrueObject->Position` (true K0 start),
2. fallback: true reco-parent end (`recoParent->TrueObject->PositionEnd`) when needed.

---

## 6) Key physics/analysis intent behind recent requests

Recent user intent has been:
1. Keep neutral reconstruction rooted in reco information.
2. Do not create neutral wrappers when parent endpoint quality is invalid.
3. Avoid classifying signal through vertex-only fallbacks in candidate category fill.
4. Quantify creation-vertex quality explicitly (raw vs projected endpoint residuals to truth).

---

## 7) Important caveats for future edits

1. `CreationVertex->PositionPandora` is currently populated with the selected end value (projected-if-valid fallback logic in `pdNeutralUtils`).
   - The new tree residual computation intentionally uses raw parent end for `k0cvtxpandora*` to match user semantics.
2. `GetSignalCategoryCodeForAnnihilationVertex(...)` still carries vertex-level subtype codes (3/4).
   - If strict policy becomes “signal logic only via neutral candidates everywhere”, this helper usage in `pdAnnihilationUtils` should be revisited.
3. The legacy workflow document `WORKFLOW_neutralKaonAnalysis.md` contains broader historical architecture and may not match this branch’s exact current builder path.

---

## 8) Practical edit map (where to touch for common requests)

- Change candidate-building constraints:
  `src/pdNeutralUtils/src/pdNeutralUtils.cxx`

- Change annihilation pairing/vertex construction:
  `src/pdNeutralUtils/src/pdAnnihilationUtils.cxx`

- Change signal definitions/categories:
  `src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.cxx`

- Add/remove microtree variables:
  `src/neutralKaonAnalysis/src/neutralKaonTree.hxx/.cxx` and registration call in `neutralKaonAnalysis.cxx`

- Event-display behavior:
  `src/neutralKaonAnalysis/src/neutralKaonEventDisplay.cxx`

---

## 9) Current near-term goals (recommended backlog)

1. **Unify semantics of creation reco positions**
   - Decide whether `CreationVertex->PositionPandora` should always be raw parent end and `Position` always projected, then enforce explicitly at construction.

2. **Signal-policy cleanup**
   - If desired, remove or isolate vertex-only fallback (`3/4`) from helper paths used outside category fill.

3. **Reco-only purity audit**
   - Keep truth usage exclusively for monitoring/labels, not for reco candidate existence or reconstruction cuts.

4. **Validation plots for new creation residuals**
   - Compare `k0cvtxpandora*` vs `k0cvtxfit*` to quantify impact of projection correction.

5. **Documentation harmonization**
   - Update/replace `WORKFLOW_neutralKaonAnalysis.md` with branch-accurate flow once this path is finalized.

---

## 10) Build / run quick checks

From repository root (`highlandPD` build already configured):

```bash
cd build
make install
```

Then run analysis as usual (example):

```bash
cd ../..
source setup.sh
cd DATA
neutralKaon.exe -o test.root minitree_pure_signal.root
root -l test.root
```

---

## 11) Short status snapshot for next AI

- Branch is active and up to date: `neutralKaon/creation_vertex/parent`.
- Neutral creation now guarded against invalid/null parent-end conditions.
- Candidate `signal` category fill is neutral-based (no vertex fallback in fill function).
- Creation-vertex residual microtree block is implemented and connected.
- Remaining cleanup is mostly semantic consistency and validation, not missing plumbing.
