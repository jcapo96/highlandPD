#ifndef neutralKaonTruthTree_h
#define neutralKaonTruthTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"

namespace neutralKaonTruthTree{

  // Minimal truth-tree content for true K0s -> pi+ pi-.
  void AddNeutralKaonTruthVariables(OutputManager& output, UInt_t nmax);
  void FillNeutralKaonTruthVariables(OutputManager& output,
                                     bool parentHasRecoObject,
                                     bool daughter1HasRecoObject,
                                     bool daughter2HasRecoObject,
                                     Float_t daughterRecoStartDistance,
                                     Int_t daughterRecoParentIdDifference,
                                     Int_t daughter1RecoParentVsTrueParentRecoIdDiff,
                                     Int_t daughter2RecoParentVsTrueParentRecoIdDiff,
                                     Float_t parentRecoTrueStartDistance,
                                     Float_t parentRecoTrueEndDistance,
                                     Float_t daughter1RecoTrueStartDistance,
                                     Float_t daughter1RecoTrueEndDistance,
                                     Float_t daughter2RecoTrueStartDistance,
                                     Float_t daughter2RecoTrueEndDistance,
                                     Float_t daughterPairFitMinDistance,
                                     Float_t daughter1OtherPairFitMinDistance,
                                     Float_t daughter2OtherPairFitMinDistance,
                                     Float_t otherPairsFitMinDistanceGlobal,
                                     Float_t daughter1RecoProtonChi2Ndf,
                                     Float_t daughter1RecoKaonChi2Ndf,
                                     Float_t daughter1RecoPionChi2Ndf,
                                     Float_t daughter2RecoProtonChi2Ndf,
                                     Float_t daughter2RecoKaonChi2Ndf,
                                     Float_t daughter2RecoPionChi2Ndf,
                                     Float_t k0TrueMomentum,
                                     Float_t parentTrueStartMomentum,
                                     Float_t parentTrueEndMomentum,
                                     Float_t daughter1TrueStartMomentum,
                                     Float_t daughter1TrueEndMomentum,
                                     Float_t daughter2TrueStartMomentum,
                                     Float_t daughter2TrueEndMomentum,
                                     Float_t k0TrueEnergy,
                                     Float_t daughter1TrueEnergy,
                                     Float_t daughter2TrueEnergy,
                                     Float_t daughter1TrueSubRecoEnergySum,
                                     Float_t daughter2TrueSubRecoEnergySum,
                                     Int_t k0TruePdg,
                                     Int_t parentTruePdg,
                                     Int_t daughter1TruePdg,
                                     Int_t daughter2TruePdg,
                                     Float_t daughter1TrueLength,
                                     Float_t daughter1RecoLength,
                                     Float_t daughter2TrueLength,
                                     Float_t daughter2RecoLength,
                                     Float_t k0TrueLength,
                                     Int_t daughter1RecoNHits,
                                     Int_t daughter2RecoNHits,
                                     Int_t daughter1TrueSubRecoQualCount,
                                     Int_t daughter1TrueSubRecoQualHitsTot,
                                     Int_t daughter2TrueSubRecoQualCount,
                                     Int_t daughter2TrueSubRecoQualHitsTot,
                                     Int_t parentRecoNHits,
                                     Float_t parentTrueLength,
                                     Float_t parentRecoLength,
                                     Int_t daughter1OtherMinSepTruePdg,
                                     Int_t vertexMinIsK0DaughterPair,
                                     Int_t parentRecoZStartGreaterThanEnd,
                                     Int_t daughter1RecoZStartGreaterThanEnd,
                                     Int_t daughter2RecoZStartGreaterThanEnd,
                                     Int_t parentIsPandoraBeam,
                                     Int_t parentExactlyOneRecoProtonNearK0Creation);

  // enum with unique indexes for output tree variables (true K0s -> pi+ pi- only; see FillTruthTree)
  enum enumNeutralKaonTruthMicroTrees{
    ntruek0 = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1, // Counter: truth-tree rows for this channel (incremented per filled K0)
    k0parhasrecoobject, // 1 if true K0 parent has a reco track with valid start/end positions and hits; else 0
    k0dau1hasrecoobject, // Same for true daughter 1 (Daughters[0] order); else 0
    k0dau2hasrecoobject, // Same for true daughter 2 (Daughters[1] order); else 0
    k0daurecostartdistance, // 3D distance between reco PositionStart of the two daughters (-999 if unavailable)
    k0daurecoparentiddiff, // daughter1 reco ParentID minus daughter2 reco ParentID (-999 if either reco missing)
    k0dau1recoparentvstrueparentrecoiddiff, // Reco parent of dau1 UniqueID minus UniqueID of reco matched to true parent of K0 (-999 if missing)
    k0dau2recoparentvstrueparentrecoiddiff, // same for daughter 2 (-999 if missing)
    k0parrecotruestartdist, // 3D |reco start - true start| for true K0 parent (-999 if invalid)
    k0parrecotrueenddist, // 3D |reco end - true end| for true K0 parent (-999 if invalid)
    k0dau1recotruestartdist, // 3D |reco start - true start| for daughter 1 (-999 if invalid)
    k0dau1recotrueenddist, // 3D |reco end - true end| for daughter 1 (-999 if invalid)
    k0dau2recotruestartdist, // 3D |reco start - true start| for daughter 2 (-999 if invalid)
    k0dau2recotrueenddist, // 3D |reco end - true end| for daughter 2 (-999 if invalid)
    k0daupairfitmindist, // Closest approach of annihilation-style line fits for the two daughter reco tracks (-999 if not computed; only when starts within AnnihilationVertexRadius)
    k0dau1otherpairfitmindist, // Minimum of those line-separation distances for daughter1 vs any other reco (not the other K0 daughter), within vertex radius (-999 if none)
    k0dau2otherpairfitmindist, // Same for daughter2 vs other reco tracks (-999 if none)
    k0otherpairsfitmindistglobal, // Minimum of k0dau1otherpairfitmindist and k0dau2otherpairfitmindist when both valid (-999 if neither valid)
    k0dau1recoprotonchi2ndf, // Chi2PID(2212)/npts collection plane (-999 if invalid)
    k0dau1recokaonchi2ndf, // Chi2PID(321)/npts collection plane (-999 if invalid)
    k0dau1recopionchi2ndf, // Chi2PID(211)/npts collection plane, same as annihilation vertex logic (-999 if invalid)
    k0dau2recoprotonchi2ndf, // Chi2PID(2212)/npts collection plane (-999 if invalid)
    k0dau2recokaonchi2ndf, // Chi2PID(321)/npts collection plane (-999 if invalid)
    k0dau2recopionchi2ndf, // Chi2PID(211)/npts collection plane, same as annihilation vertex logic (-999 if invalid)
    k0truemomentum, // True |p| at K0 trajectory start (MC Momentum; -999 if unavailable)
    k0trueenergy, // True total energy sqrt(p^2+m^2) at K0 trajectory start (GeV; -999 if unavailable)
    k0partruestartmom, // True |p| at start for true K0 parent (-999 if no parent)
    k0partrueendmom, // True |p| at end for true K0 parent (-999 if no parent)
    k0dau1truestartmom, // True |p| at start for true daughter 1 (Daughters[0])
    k0dau1trueendmom, // True |p| at end for true daughter 1
    k0dau2truestartmom, // True |p| at start for true daughter 2 (Daughters[1])
    k0dau2trueendmom, // True |p| at end for true daughter 2
    k0dau1trueenergy, // True total energy at start for daughter 1 (GeV; -999 if unavailable)
    k0dau2trueenergy, // True total energy at start for daughter 2 (GeV; -999 if unavailable)
    k0truepdg, // True PDG of the K0 (this row)
    k0partruepdg, // True PDG of K0 parent (-999 if parent not found in spill)
    k0dau1truepdg, // True PDG of daughter 1 (Daughters[0])
    k0dau2truepdg, // True PDG of daughter 2 (Daughters[1])
    k0dau1truelength, // |PositionEnd-Position| for daughter 1 (cm; same as ana tree; -999 if invalid)
    k0dau1recolength, // Reco Length for daughter 1 (cm; -999 if no reco)
    k0dau2truelength, // |PositionEnd-Position| for daughter 2 (cm; -999 if invalid)
    k0dau2recolength, // Reco Length for daughter 2 (cm; -999 if no reco)
    k0truelength, // |PositionEnd-Position| for K0 (cm; -999 if invalid)
    k0dau1recohits, // Reco NHits for daughter 1 (-999 if no reco)
    k0dau2recohits, // Reco NHits for daughter 2 (-999 if no reco)
    k0dau1truensubreco, // Number of true sub-daughters of daughter 1 with reco valid start/end + hits
    k0dau1truensubrecohitstot, // Sum of reco NHits over those sub-daughters (0 if none)
    k0dau2truensubreco, // Same for true daughter 2
    k0dau2truensubrecohitstot,
    k0dau1truesubrecoenergy, // Sum true E at start (GeV) over direct true daughters of daughter 1 with valid reco (0 if none)
    k0dau2truesubrecoenergy, // Same for daughter 2
    k0parrecohits, // Reco NHits for true K0 parent (-999 if no reco)
    k0partruelength, // |PositionEnd-Position| for K0 parent (cm; -999 if no parent / invalid verts)
    k0parrecolength, // Reco Length for K0 parent (cm; -999 if no reco)
    k0dau1otherminpairtruepdg, // True PDG of reco particle attaining k0dau1otherpairfitmindist (-999 if none)
    k0vertexminisdaupair, // 1 if min of {dau-dau, dau1-other, dau2-other} fit-line distances equals dau-dau (and dau-dau valid); else 0
    k0parrecozstartgtzend, // 1 if reco PositionStart[2] > PositionEnd[2] for true-K0 parent track (0 if no/invalid reco Z)
    k0dau1recozstartgtzend, // same for daughter1 reco
    k0dau2recozstartgtzend, // same for daughter2 reco
    k0parispandorabeam, // 1 if reco matched to true K0 parent has isPandora (Pandora beam); else 0
    k0pardonerecoprotonatcreation, // 1 iff exactly one true daughter of K0's true parent is proton (2212) with valid reco and reco start within CreationVertexRadius of true K0 birth; else 0

    enumNeutralKaonTruthMicroTreesLast_neutralKaonTruthTree
  };
}

#endif