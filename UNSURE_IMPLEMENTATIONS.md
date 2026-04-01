# Unsure Implementations Log

This file tracks code changes that are experimental or not fully validated yet.
If needed, you can ask to remove any item by its `ID`.

---

## ID: NK-MOMDIAG-003

- **Title**: Temporary daughter momentum-method diagnostics in K0 vertex microtree
- **Date**: 2026-03-31
- **Status**: Active (tentative)
- **Intent**: Diagnose why extension-fit momentum assignment is rare by exposing per-daughter decision and fit-quality diagnostics.

### What it does

Adds temporary per-vertex/per-daughter diagnostic branches:

- Existing method flags:
  - `k0vtxdaughter1momentummethod`
  - `k0vtxdaughter2momentummethod`
- New diagnostics:
  - `k0vtxdau1haspreexistingmomentum`, `k0vtxdau2haspreexistingmomentum`
  - `k0vtxdau1extensionattempted`, `k0vtxdau2extensionattempted`
  - `k0vtxdau1extensionvalid`, `k0vtxdau2extensionvalid`
  - `k0vtxdau1extensionchi2ndf`, `k0vtxdau2extensionchi2ndf`
  - `k0vtxdau1extensionnvalidhits`, `k0vtxdau2extensionnvalidhits`

### Where implemented

1. `highlandPD/src/pdUtils/src/pdAnnihilationUtils.cxx`
   - Captures and stores daughter-level extension diagnostics while assigning momentum.
2. `highlandPD/src/pdUtils/src/pdMomEstimation.hxx`
   - Added `EstimateMomentumWithExtensionDetailed(...)` declaration.
3. `highlandPD/src/pdUtils/src/pdMomEstimation.cxx`
   - Added `EstimateMomentumWithExtensionDetailed(...)` implementation and shared template-loading helper.
4. `highlandPD/src/pdEventModel/src/pdDataClasses.hxx`
   - Added new diagnostic members to `AnaAnnihilationVertexPD`.
5. `highlandPD/src/pdEventModel/src/pdDataClasses.cxx`
   - Initializes new diagnostic members.
6. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
   - Added new microtree enum entries.
7. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
   - Added branch definitions and fill logic for all diagnostics.

### Important notes

- These branches are diagnostic-only and intended for short-term validation.
- `extensionchi2ndf` is set to `-999` when `ndf <= 0` or extension was not attempted.

### How to test quickly

1. Run:
   - `neutralKaon.exe -o test.root minitree_pure_signal.root`
2. In ROOT:
   - `root -l test.root`
   - `default->Draw("k0vtxdaughter1momentummethod")`
   - `default->Draw("k0vtxdau1extensionvalid:k0vtxdau1extensionattempted","","colz")`
   - `default->Draw("k0vtxdau1extensionchi2ndf","k0vtxdau1extensionattempted==1")`

### Rollback plan (full removal)

1. Remove diagnostic members from:
   - `highlandPD/src/pdEventModel/src/pdDataClasses.hxx`
   - `highlandPD/src/pdEventModel/src/pdDataClasses.cxx`
2. Remove diagnostic capture logic from:
   - `highlandPD/src/pdUtils/src/pdAnnihilationUtils.cxx`
3. Remove detailed extension API if no longer needed:
   - `highlandPD/src/pdUtils/src/pdMomEstimation.hxx`
   - `highlandPD/src/pdUtils/src/pdMomEstimation.cxx`
4. Remove diagnostic microtree vars from:
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
5. Rebuild/install:
   - `cmake --build . -j8 && cmake --install .` (from `highlandPD/build`)

---

## ID: NK-XBIAS-001

- **Title**: Global reconstructed X shift using `VertexXBiasCorrectionCm`
- **Date**: 2026-03-27
- **Status**: Active (tentative)
- **Intent**: Correct observed positive X bias in reconstructed neutral-kaon vertex position.

### What it does

Applies a configurable X shift to reconstructed particle spatial quantities:

- `x_corrected = x_reco - VertexXBiasCorrectionCm`
- Parameter:
  - `neutralKaonAnalysis.VertexXBiasCorrectionCm`
  - defined in `neutralKaonAnalysis.parameters.dat`

### Where implemented

1. `highlandPD/src/neutralKaonAnalysis/parameters/neutralKaonAnalysis.parameters.dat`
   - Added:
     - `neutralKaonAnalysis.VertexXBiasCorrectionCm`

2. `highlandPD/src/pdCorrections/src/ParticlePositionSCECorrection.cxx`
   - Applies X shift to reconstructed particle objects:
     - `PositionStart[0]`, `PositionEnd[0]`
     - `PositionStartSCE[0]`, `PositionEndSCE[0]`
     - `Hits[plane][i].Position.X()`
     - `Hits[plane][i].Position_NoSCE.X()`
     - `TrjPoints[i].Position.X()`
     - `TrjPoints[i].Position_NoSCE.X()`

### Important notes

- This correction is applied in the correction stage, after SCE correction flow.
- Previous local X-shift hooks in `pdAnalysisUtils` were removed to avoid double-shifting.
- Default parameter value is currently `0.0` (no shift unless explicitly enabled).

### How to test quickly

1. Set:
   - `neutralKaonAnalysis.VertexXBiasCorrectionCm = 0.92`
2. Re-run:
   - `neutralKaon.exe -o test.root minitree_pure_signal.root`
3. Check means:
   - `k0vtxfitresidualx`
   - `k0vtxpandoraresidualx`

### Rollback plan (full removal)

1. Remove parameter from:
   - `highlandPD/src/neutralKaonAnalysis/parameters/neutralKaonAnalysis.parameters.dat`
2. Remove X-bias block from:
   - `highlandPD/src/pdCorrections/src/ParticlePositionSCECorrection.cxx`
3. Rebuild/install:
   - `cmake --build . -j8 && cmake --install .` (from `highlandPD/build`)

---

## ID: NK-HITPROF-002

- **Title**: Signal-only hit-to-true-line profiles vs traveled distance
- **Date**: 2026-03-27
- **Status**: Active (tentative)
- **Intent**: Characterize reconstruction quality along daughter tracks with finer near-vertex granularity.

### What it does

Creates external histogram/profile objects for each annihilation-vertex daughter:

- Distance from each collection-plane hit to the daughter true line:
  - true line defined by `true startPos + true Direction`
- Filled versus traveled distance from reco start to hit position.

Objects:

- `h_k0dau1_hitDistToTrueLine_vs_travel_2d` (`TH2F`)
- `h_k0dau2_hitDistToTrueLine_vs_travel_2d` (`TH2F`)
- `p_k0dau1_hitDistToTrueLine_vs_travel` (`TProfile`)
- `p_k0dau2_hitDistToTrueLine_vs_travel` (`TProfile`)

### Where implemented

1. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
   - Added helper logic to build/fill per-daughter hit-to-true-line hist/profile objects.
   - Fill is gated to signal candidates only via:
     - `neutralKaonAnaUtils::IsSignalCandidate(candidate, event)`
2. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
   - Added declaration for profile write-out helper.
3. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.hxx`
   - Added `Finalize()` override declaration.
4. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`
   - `Finalize()` writes these external objects to the output ROOT file.

### Important notes

- Histograms are external ROOT objects (not microtree branches).
- **Related**: signal-only dE/dx vs residual range and per-track graphs are documented in **`NK-DEDX-RR-004`** (same fill gate and same `Finalize` write path).
- Current traveled-distance binning is variable:
  - `0-30 cm`: `2 cm` bins
  - `30-500 cm`: `5 cm` bins
- Only signal candidates contribute (background candidates are skipped).

### How to test quickly

1. Run:
   - `neutralKaon.exe -o test.root minitree_pure_signal.root`
2. In ROOT:
   - `root -l test.root`
   - `.ls`
   - Check objects exist:
     - `h_k0dau1_hitDistToTrueLine_vs_travel_2d`
     - `h_k0dau2_hitDistToTrueLine_vs_travel_2d`
     - `p_k0dau1_hitDistToTrueLine_vs_travel`
     - `p_k0dau2_hitDistToTrueLine_vs_travel`
3. Draw:
   - `h_k0dau1_hitDistToTrueLine_vs_travel_2d->Draw("COLZ")`
   - `p_k0dau1_hitDistToTrueLine_vs_travel->Draw()`

### Rollback plan (full removal)

1. Remove hit-profile helper/static objects and fill calls from:
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
2. Remove write helper declaration from:
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
3. Remove `Finalize()` hook from:
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.hxx`
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`
4. Rebuild/install:
   - `cmake --build . -j8 && cmake --install .` (from `highlandPD/build`)

---

## ID: NK-DEDX-RR-004

- **Title**: Signal-only dE/dx vs residual range diagnostics and temporary Bragg scoring prototype
- **Date**: 2026-03-31
- **Status**: Removed
- **Intent**: Study stopping vs interacting pions at the annihilation vertex using collection-plane hit dE/dx vs residual range and temporary Bragg-style discriminators.

### What it does

- Historical note only: this item introduced and iterated several Bragg-score definitions and temporary graph outputs.
- Current code keeps only signal-gated dE/dx-vs-RR `TH2F` diagnostics by true daughter end process (2, 9, 10, 11); Bragg-score branches and related logic were removed.

### Where implemented

1. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.hxx` / `.cxx`
   - `GetSignalTrueParent(...)`: refactored signal parent access; `IsSignalCandidate` calls it.
2. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.cxx`
   - Active: signal-gated dE/dx-vs-RR `TH2F` by true daughter end process.
   - Removed: Bragg-score helpers/branches and per-track `TGraph` outputs.
3. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonTree.hxx`
   - Active enum set after Bragg removal.
4. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysis.cxx`
   - `Finalize()` calls `neutralKaonTree::WriteHitDistanceProfiles(output())`.

### Important notes

- Marked `Removed` to avoid confusion with current tree content.

### How to test quickly

1. Confirm Bragg branches are absent in new files and only dE/dx-vs-RR process histograms remain.

### Rollback plan (full removal)

N/A (already removed).

---

## ID: NK-SIGNAL-SPLIT-005

- **Title**: Signal category split by stopping behavior of true pion daughters
- **Date**: 2026-04-01
- **Status**: Active (tentative)
- **Intent**: Keep baseline signal definition while exposing a stopping/interacting signal subtype split for quick studies and plotting.

### What it does

Within the existing `signal` object category, signal candidates are split into:

- `two_stopping` (code `1`): both daughters have `ProcessEnd` in `{2, 11}`
- `one_stopping` (code `5`): exactly one daughter has `ProcessEnd` in `{2, 11}`
- `interacting` (code `6`): neither daughter has `ProcessEnd` in `{2, 11}`

Non-signal categories remain:
- `background` (`2`)
- `legit_vertex_2body` (`3`)
- `legit_vertex_multibody` (`4`)

### Where implemented

1. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.cxx`
   - Added stopping-process helper and subtype code assignment in `FillSignalCandidateCategory(...)`.
   - Updated `AddSignalCandidateCategory()` labels/codes/colors.
2. `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.hxx`
   - No API expansion needed; behavior changed in implementation only.

### Important notes

- Base signal definition (`IsSignalCandidate`) is unchanged.
- Split uses true daughter end-process info and is intended for classification/validation studies.

### How to test quickly

1. Run:
   - `neutralKaon.exe -o test.root minitree_pure_signal.root`
2. In ROOT:
   - `draw.Draw(ana, "signal", 10, 0, 10, "signal", "")`
   - Check bins corresponding to codes `1,5,6` are populated as expected.

### Rollback plan

1. Revert `AddSignalCandidateCategory()` arrays and subtype assignment in:
   - `highlandPD/src/neutralKaonAnalysis/src/neutralKaonAnalysisUtils.cxx`
2. Rebuild/install:
   - `cmake --build . -j8 && cmake --install .` (from `highlandPD/build`)

---

## Template for new entries

- **ID**: `<short-unique-id>`
- **Title**:
- **Date**:
- **Status**: `Active (tentative)` | `Validated` | `Removed`
- **Intent**:
- **What it does**:
- **Where implemented**:
- **Important notes**:
- **How to test quickly**:
- **Rollback plan**:
