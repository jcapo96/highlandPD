# \(K^0_S\) Analysis: Strategy, Physics Motivation, and Current Implementation

This note compiles the **analysis strategy** and **physics motivation** (aligned with *K0 Analysis Strategy*) with the **current highlandPD implementation** as of the repository state. It is split into **geometry** (annihilation vertex, neutral particle, creation vertex) and **momentum** (TLE-style free-range dE/dx, MCS, and kinematic / joint fit).

**Code entry point (neutral kaon candidate building):** `highlandPD/src/neutralKaonAnalysis/src/neutralKaonSelection.cxx` — `FindNeutralCandidatesAction::Apply`.

---

## 1. Physics motivation

### 1.1 What is being reconstructed

- **\(K^0_S\)** is a neutral short-lived state; it does **not** leave a continuous ionization track in the LArTPC.
- **Identification** relies on inferring its **trajectory** from a **production (creation) vertex** and a **decay vertex**, and on characterizing the **two charged daughters** from the dominant hadronic decay  
  \[
  K^0_S \to \pi^+ + \pi^-
  \]
  which produces a characteristic **V-shaped** signature (two outgoing tracks from a common point).

### 1.2 Production context (beam, charge exchange)

- \(K^0_S\) can be produced in beam–nucleus interactions, e.g. channels involving a **\(K^+\)** and a **\(K^0\)** (including charge-exchange–like topologies where a **proton** in the final state is prominent).
- If the **incoming \(K^+\)** momentum is constrained and the **proton** momentum is well measured, the **invisible neutral** momentum can be strongly constrained (exploiting known masses and kinematics). The analysis strategy document highlights this as motivation for **creation-vertex refinement** that can associate a **proton-like** secondary with the parent track.

### 1.3 Relevance to DUNE

- \(K^0_S\) appears in final states relevant to **nucleon decay** and related searches (e.g. modes involving \(K^0\) / \(K^0_S\) alongside charged hadrons or leptons).
- Correct **\(K^0_S\)** reconstruction can improve **final-state reconstruction** (e.g. disentangling backgrounds or recovering momentum balance) when a \(K^0_S\) is present in **FSI** or signal-like topologies.
- Typical momenta cited in the strategy context are of order **hundreds of MeV/c** for the \(K^0_S\) and **~200 MeV/c** for decay pions in the rest frame; the implementation should therefore support **sub-GeV** pion kinematics and robust vertexing in the TPC.

---

## 2. Analysis strategy (conceptual)

1. **Decay vertex (“annihilation” in code):** find **pairs of reconstructed tracks** consistent with two pions originating from the same point, assign a **3D decay vertex** and resolve **track sharing** when multiple pairs compete.
2. **Creation vertex:** infer where the **parent** (beam-related) track ends relative to the decay vertex; optionally **refine** the creation point using a **cylinder** along the parent and a **secondary track** (motivated by \(K^+ \to K^0 + p\)-like topologies with a visible **proton**).
3. **Neutral particle object:** package **one decay vertex + one creation hypothesis** into an **`AnaNeutralParticlePD`**, with a **poly-line** from creation to decay for length and axis-based quality variables.
4. **Momentum:** assign **pion** momenta to the two decay daughters (and optionally **proton** momentum to the creation-side secondary), using **calorimetric / range–dE/dx** information and, when enabled, **multiple scattering** and a **joint two-body kinematic** hypothesis at the decay vertex.

The following sections describe **how this is implemented today** in `highlandPD`.

---

## 3. Geometric reconstruction (current implementation)

Parameters are primarily under the **`neutralKaonAnalysis.*`** namespace (`neutralKaonAnalysis.parameters.dat`).

### 3.1 Annihilation vertex (`AnaAnnihilationVertexPD`)

**Role:** reconstructed **decay vertex** of the neutral candidate: two daughter tracks and two geometric estimates of the vertex (Pandora-direction vs hit-based extrapolation).

**API:** `pdAnnihilationUtils::CreateVertices` → `CreateVerticesCommon`  
**Headers / types:** `pdAnnihilationUtils.hxx`, `pdNeutralDataClasses.hxx` (`AnaAnnihilationVertexPD`).

**Procedure (summary):**

1. **Pair enumeration** over `event.Particles` (unordered pairs).
2. **Per-track guards:** valid `PositionStart`; minimum collection-plane hits `NHitsPerPlane[2] > AnnihilationVertexMinCollectionHits`; PID veto **`IsProtonLikeAndNotPionLike`** using **`pdAnaUtils::Chi2PIDChi2PerHit`** (proton 2212 and pion 211 χ²/hit thresholds from parameters).
3. **Parent consistency:** `daughter1->ParentID == daughter2->ParentID`.
4. **Endpoint distance:** if `AnnihilationVertexUseEndpointCombinatorics != 0`, **`BestAnnihilationPairEndpoint`** minimizes 3D distance over SS/SE/ES/EE Pandora endpoints; else start–start only. Pairs must satisfy distance \(\le\) **`AnnihilationVertexRadius`** (`maxDaughterDistance` argument).
5. **Vertex position — pre overlap-filter:**
   - **`FillPositionPandora`:** infinite lines from each daughter’s **`PositionStart`** + **`DirectionStart`**; **skew-line closest approach**; vertex = **midpoint** of closest points; stores `MinimumDistancePandora`.
   - **`FillPositionFit`:** **`pdAnaUtils::ExtrapolateTrack`** over `[TrackFitDistanceFromStart, TrackFitDistanceFromStart + TrackFitLength]` from the **track start**; **`pdAnaUtils::FindClosestPointsBetweenLines`**; vertex = **midpoint**; stores `MinimumDistanceFit`.
6. **Overlap resolution:** **`FilterVerticesByMinimumDistanceFit`** sorts by **`MinimumDistanceFit`** (ascending) and keeps each **`AnaParticlePD`** in **at most one** vertex; conflicting vertices are **deleted**.
7. **Endpoint reversal:** **`ApplyAnnihilationPairingReversal`** for non-SS combos may **`ReverseParticleIfStartZGreaterThanEndZ`** (swap start/end, directions, SCE, momentum, reverse `TrjPoints` / `Hits`, optional **`EstimateHitsDirection`** / **`ComputeResidualRange`**), then **recompute** Pandora and fit vertex positions.
8. **Degeneracy:** **`ComputeAnnihilationVertexDegeneracy`** counts other tracks whose **raw or projected** activity is consistent with lying “on” the vertex within configured radii and line–vertex distances (**`AnnihilationVertexDegeneracy*`** parameters). When a creation position is supplied later in other APIs, additional **creation vs annihilation** proximity logic can veto ambiguous tracks (see annihilation utils).

**Outputs on `AnaAnnihilationVertexPD`:** `PositionPandora`, `PositionFit`, closest points, `MinimumDistance*`, `Degeneracy`, daughter list, pairing metadata, plus momentum-related fields filled in the momentum stage (Section 4).

---

### 3.2 Creation vertex (`AnaCreationVertexPD`)

**Role:** reconstructed **production** side of the neutral: primarily the **parent track end** (with optional tail projection), optionally refined to a **beam–proton** closest-approach point inside a **cylinder** between creation anchor and decay.

**Important:** there is **no** separate `CreateCreationVertices` pass in the running code. Creation vertices are built **inside** `pdNeutralUtils::CreateNeutralsFromAnnihilationVertices`.

**API:** `pdNeutralUtils::CreateNeutralsFromAnnihilationVertices`  
**Helpers:** `pdCreationUtils::ProjectBeamTailOntoStartDirection`, `pdCreationUtils::ProjectPointOntoLine` (`pdCreationUtils.hxx` / `.cxx`); local helpers in `pdNeutralUtils.cxx`: **`FindCommonRecoParentForAnnihilationVertex`**, **`BuildRecoFitLine`**, **`IsInsideCreationCylinder`**, **`HasValidRecoPoint3`**.

**Procedure (summary):**

1. **Reference decay position `annihilationPos`:** prefer **`PositionFit`** on the annihilation vertex, else **`PositionPandora`**.
2. **Parent list `parentParticles`:**
   - **`CreationVertexBeamParticleMode == 1`:** single candidate with **`UniqueID == daughters’ common ParentID`**, valid end, and **`PositionEnd[2] < annihilationPos.Z()`** if `annihilationPos` is valid.
   - **Mode 0 (default):** all particles except the two daughters with valid **`PositionEnd`** and **`PositionEnd.Z() < annihilationPos.Z()`** (requires valid `annihilationPos`).
3. **For each parent:** `new AnaCreationVertexPD`; **`BeamParticle = parent`**.
4. **Beam end anchor:** optional **`ProjectBeamTailOntoStartDirection`** using **`TrackFitDistanceCreationFromEnd`** (caller) and internal **`TrackFitCreationLength`** / extrapolation from the **track end** (`pdAnaUtils::ExtrapolateTrack` with `fromStart=false`); else raw **`PositionEnd`**.
5. **Z veto:** if annihilation **`PositionPandora`** is valid and **`parentEndZ >= PositionPandora[2]`**, skip this parent.
6. **Default creation position:** set **`Position`** and **`PositionPandora`** to the beam end; **`Particles = {parent}`**; **`UniqueID = parent->UniqueID`**.
7. **Refinement (optional):** if `CreationVertexRadius > 0` and fit parameters allow, build **beam fit line** from **`BuildRecoFitLine`** (`ExtrapolateTrack` from **start**). Scan other tracks:
   - valid start; **proton-like PID** cut via **`Chi2PIDChi2PerHit(..., 2212)`** vs **`CreationVertexSecondParticleMaxProtonChi2Ndf`**;
   - start inside **`IsInsideCreationCylinder`** around the beam segment between **beam end** and **`annihilationPos`**;
   - best **`FindClosestPointsBetweenLines`** between beam line and second track line → update **`Position` / `PositionPandora`** to **midpoint** of closest points; store **`SecondParticle`**, scores, **`FittedLineParams`**, extend **`Particles`**.
   - Optionally **`pdAnnihilationUtils::AssignProtonMomentumFromResidualRange(second)`** if annihilation daughters received momentum methods.
8. **Degeneracy:** **`pdAnnihilationUtils::ComputeCreationVertexDegeneracy(event, creationVtx, annihilationVtx, -1)`** (reuses **annihilation degeneracy** geometry parameters; excludes annihilation daughters from the count).

---

### 3.3 Neutral particle (`AnaNeutralParticlePD`)

**Role:** **analysis-level candidate** tying **one annihilation vertex** and **one creation vertex** (and a chosen **parent**), plus derived **axis**, **lengths**, and **alignment** diagnostics.

**API:** same as creation — produced only in **`CreateNeutralsFromAnnihilationVertices`**.

**Procedure (summary):**

1. **`new AnaNeutralParticlePD`** per surviving parent path:
   - **`AnnihilationVertex`**, **`CreationVertex`**, **`Parent = creationVtx->BeamParticle`**.
2. **Endpoints:**
   - **`PositionStart`** from **`creationVtx->Position`** (if valid).
   - **`PositionEnd`** from **`annihilationVtx->PositionPandora`** (note: **Pandora** decay position is used for the neutral’s end, not necessarily `PositionFit`).
3. **Lengths / direction:** **`Length` / `LengthPandora`** from start–Pandora end; **`LengthFit`** from start–**`PositionFit`** when both valid.
4. **Alignment:** **`pdAnnihilationUtils::FillNeutralParticleAlignment`** — angle between **creation→decay axis** and **vector sum of daughter momenta** using Pandora vs fit directions (`DaughterPandoraAndFitDirs`, **`ExtrapolateTrack`**).
5. **Truth link (MC):** **`GetIntermediateTrueParticle`** may set **`TrueObject`** when the true hierarchy is self-consistent (diagnostics / performance studies).
6. **Winner selection:** if **`SelectSingleNeutralPerAnnihilationVertex`** is enabled, **`SelectWinnerNeutralCandidateIndex`** scores candidates (length, beam alignment, parent hit count); non-winners are **`delete`**d with their **`CreationVertex`**.

---

## 4. Momentum reconstruction (current implementation)

Momentum for the **two decay pions** is assigned in **`pdAnnihilationUtils::CreateVerticesCommon`** after geometry (post filter + reversal). Creation-side **proton** momentum uses **`AssignProtonMomentumFromResidualRange`** when the refinement path applies.

The implementation exposes **three complementary ideas** in code and parameters:

| Concept in this document | Primary code paths |
|--------------------------|-------------------|
| **TLE-style / free-range dE/dx** | `pdAnaUtils::GetdEdxLikelihoodFreeRangeFit` (per-track scan used by **`AssignMomentumFromResidualRange`**); **`pdMomReconstruction::EstimatePionMomentumTLEFit`** (used inside the joint fit and pion momentum analysis) |
| **MCS** | **`pdMomReconstruction::EstimatePionMomentumMCSReco`** (pion momentum analysis; **optional** contribution inside joint \(K^0_S\) grid via config) |
| **Kinematic / joint fit** | **`pdMomReconstruction::FitJointK0sTwoPionMomentaOnGrid`** + **`pdJointK0sPionMomentum::FitJointMomentaOnGrid`** with \(K^0_S\) mass prior |

Configuration bridges:

- **`pdMomReconstruction::FillJointK0sTwoPionGridFitConfig_FromNeutralKaonParams`** (`pdMomReconstructionFromParams.hxx`) loads neutral-kaon parameters into **`JointK0sTwoPionGridFitConfig`**.

---

### 4.1 TLE-style method (track-length / interior dE/dx likelihood)

**Physics idea:** along a stopping or contained track, **ionization vs residual range** (or path from an assumed **vertex offset**) follows a **pion Bethe–Bloch–like** curve. A **likelihood scan** over a hypothetical **extra range** / momentum (free-range / TLE spirit) picks the momentum that best explains the **interior collection-plane dE/dx** pattern.

**Implementation A — annihilation daughters (default fallback path):**

- **`AssignPionMomentumFromResidualRange`** → **`AssignMomentumFromResidualRange(..., PDG 211, ...)`** in `pdAnnihilationUtils.cxx`.
- Calls **`pdAnaUtils::GetdEdxLikelihoodFreeRangeFit`** (`pdUtilsDEdx.hxx` / `pdUtilsDEdx.cxx`) with scan limits and hit masks from **`neutralKaonAnalysis.FreeRange*`** parameters (e.g. `FreeRangeScanLmaxCm`, `FreeRangeDedxMinInteriorHits`, `FreeRangeDedxSkipHitsFirst/Last`, dE/dx windows for pions, optional low-momentum refinement).
- On success, sets **`particle->Momentum`** and returns method code **`kMomMethodFreeRangeML`**.

**Implementation B — explicit TLE API (shared with joint fit and pion momentum analysis):**

- **`pdMomReconstruction::EstimatePionMomentumTLEFit`** (`pdMomReconstructionTLEFit.hxx` / `.cxx`) with **`PionTLEFitConfig`** (mirrors pion-momentum-analysis style: interior hits, dE/dx bounds, scan step, fine scan, PDF path slice, etc.).
- Used by **`FitJointK0sTwoPionMomentaOnGrid`** to build **per-daughter log-likelihood vs momentum curves** on a grid.

**Outputs / diagnostics:** vertex stores **`Daughter{1,2}MomentumMethod`**, extension diagnostics, and joint-TLE fields when the joint path runs (`Daughter*MomentumTLE`, auxiliary curves in `JointK0sTwoPionFitAuxiliary`).

---

### 4.2 MCS method (multiple scattering)

**Physics idea:** **angular deviations** between successive **track segments** depend on **momentum** (roughly \(\theta \sim 1/p\)). A **likelihood** (or scan) over momentum can be built from **reco `TrjPoints`** (and optionally compared to truth SCE-ordered points in validation).

**Implementation:**

- **`pdMomReconstruction::EstimatePionMomentumMCSReco`** and **`EstimatePionMomentumMCSTrue`** (`pdMomReconstructionMCS.hxx` / `.cxx`) with **`PionMCSConfig`** (wraps **`MCSLikelihoodConfig`**, segment drops, momentum scan ceiling); the MCS likelihood surface (`BuildPionMcsScatteringSamples`, `ComputePionMcsNLLFromScatteringSamples`, `BuildPionMCSLogLikelihoodVsMomentumCurve`, …) lives in the same translation unit.
- **Primary standalone use** in **`pionMomentumAnalysis`** (`pionMomentumAnalysis.cxx`) for single-track MCS momentum studies and micro-trees.
- **Neutral kaon joint use:** **`JointK0sTwoPionGridFitConfig`** includes **`PionMCSConfig mcs`** and flag **`useMCS`** / weight **`mcsWeight`**. When enabled, the joint pipeline can incorporate MCS information alongside TLE curves (see joint implementation in `pdMomReconstructionJointK0s.cxx`).

**Outputs / diagnostics:** `PionMCSRecoBuffers` per segment; joint path can expose marginal MCS argmax momenta on vertex (`Daughter*MomentumMCS` fields in `AnaAnnihilationVertexPD` when filled in `CreateVerticesCommon`).

---

### 4.3 Kinematic fit method (joint two-pion grid with \(K^0_S\) mass)

**Physics idea:** the two pions from \(K^0_S \to \pi^+\pi^-\) are **kinematically coupled**: for given directions, their momenta are not independent — the **invariant mass** of the pair should peak near **\(m_{K^0_S}\)**. A **joint** estimator therefore combines **per-track likelihoods** (TLE and optionally MCS) with a **mass penalty** / constraint around **`mK0sMassGeV`** (default ~0.4976 GeV) with width **`sigmaMassGeV`** (and optional **event-level \(\sigma_m\)** refinement when `useEventSigmaM` is enabled).

**Implementation:**

- **`pdMomReconstruction::FitJointK0sTwoPionMomentaOnGrid`** (`pdMomReconstructionJointK0s.hxx` / `.cxx`):
  1. Builds **TLE log L vs momentum** for each daughter via **`EstimatePionMomentumTLEFit`** (optional diagnostic axes returned).
  2. Optionally incorporates **MCS** according to **`JointK0sTwoPionGridFitConfig`**.
  3. Invokes **`pdJointK0sPionMomentum::FitJointMomentaOnGrid`** (see `pdJointK0sPionMomentum.hxx`) on a **2D momentum grid** (or refined grid) with **mass-aware scoring** (`sigmaMass*`, `massPenaltyScale`, `refineFactor`, bounds `pMinGeV` / `pMaxGeV` / `pStepGeV`).
- **Gating in annihilation creation:** controlled by **`neutralKaonAnalysis.JointK0sMomentumEnable`** and optional **`EnsureMomentumSignalOnly`** (signal category from **`neutralKaonAnaUtils::GetSignalCategoryCodeForAnnihilationVertex`**).
- **On success:** updates **`daughter1->Momentum`**, **`daughter2->Momentum`**, sets **`Daughter*MomentumMethod = kMomMethodJointK0sGrid`**, fills **`JointK0s*`** diagnostic members on **`AnaAnnihilationVertexPD`**.
- **On failure / disable:** falls back to **independent** **`AssignPionMomentumFromResidualRange`** per daughter.

**Relation to “geometry”:** joint fit uses **fit directions** from **`DaughterPandoraAndFitDirs`** (from `ExtrapolateTrack` near the vertex), so it is **tied to the same track extrapolation** as **`FillPositionFit`**.

---

## 5. End-to-end call sequence (neutral kaon box)

1. **`pdAnnihilationUtils::CreateVertices`** — build **`AnaAnnihilationVertexPD`** list (geometry + degeneracy + daughter momenta).
2. **`pdNeutralUtils::CreateNeutralsFromAnnihilationVertices`** — for each annihilation vertex, build **`AnaCreationVertexPD`** + **`AnaNeutralParticlePD`** candidates (geometry + creation degeneracy + alignment + optional proton momentum on second track).

Downstream cuts and trees read **`ToyBoxNeutralKaon::neutralParticleCandidates`** and follow pointers into **`AnnihilationVertex`** / **`CreationVertex`**.

---

## 6. Map to original strategy PDF

The PDF’s **Sections 2–6** (“Identification of \(K^0_S\) trajectory”, decay vertex radius, closest-approach sorting, creation vertex from parent end + cylinder refinement, `AnaAnnihilationVertexPD` / `CreateVertices`, `AnaNeutralParticlePD` assembly) match the **current code structure** above, with these **naming** correspondences:

| PDF / physics language | Code identifier |
|--------------------------|-----------------|
| Decay vertex | **`AnaAnnihilationVertexPD`** (namespace “annihilation”) |
| Creation vertex | **`AnaCreationVertexPD`** |
| Neutral trajectory object | **`AnaNeutralParticlePD`** |

The **momentum** sections of the PDF are superseded in detail by **Section 4** of this note, which reflects **`pdAnnihilationUtils`**, **`pdMomReconstruction`**, **`pdJointK0sPionMomentum`**, and **`pdAnaUtils`** as wired for **`neutralKaonAnalysis`**.

---

## 7. Suggested use

- **For documentation:** merge Section **1–2** with your narrative from *K0 Analysis Strategy.pdf*; use **Sections 3–5** as the **implementation appendix**.
- **For maintenance:** when parameters or function names change, update **Section 3–4** in lockstep with `neutralKaonAnalysis.parameters.dat` and the files cited above.

---

*Generated to align the strategy PDF with the current `highlandPD` implementation. Adjust wording to match your thesis/AN style (e.g. replace code identifiers with footnotes) as needed.*
