# Unsure Implementations Log

This file tracks code changes that are experimental or not fully validated yet.
If needed, you can ask to remove any item by its `ID`.

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
