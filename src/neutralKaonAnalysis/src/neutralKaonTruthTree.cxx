#include "neutralKaonTruthTree.hxx"

//********************************************************************
void neutralKaonTruthTree::AddNeutralKaonTruthVariables(OutputManager& output, UInt_t nmax){
    AddVarMaxSizeVI(output, k0parhasrecoobject, "K0 true parent has valid reco object", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1hasrecoobject, "K0 true daughter1 has valid reco object", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2hasrecoobject, "K0 true daughter2 has valid reco object", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0daurecostartdistance,
                     "Reco distance between K0 daughter PositionStart (cm)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0daurecoparentiddiff,
                    "Reco ParentID difference between K0 daughters (dau1 - dau2)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1recoparentvstrueparentrecoiddiff,
                    "Daughter1 reco parent UniqueID minus UniqueID of reco matched to true parent of K0", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2recoparentvstrueparentrecoiddiff,
                    "Daughter2 reco parent UniqueID minus UniqueID of reco matched to true parent of K0", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0parrecotruestartdist,
                    "Distance reco vs true PositionStart for K0 true parent (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0parrecotrueenddist,
                    "Distance reco vs true PositionEnd for K0 true parent (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recotruestartdist,
                    "Distance reco vs true PositionStart for K0 daughter1 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recotrueenddist,
                    "Distance reco vs true PositionEnd for K0 daughter1 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recotruestartdist,
                    "Distance reco vs true PositionStart for K0 daughter2 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recotrueenddist,
                    "Distance reco vs true PositionEnd for K0 daughter2 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0daupairfitmindist,
                    "Minimum distance between fitted lines for K0 daughter pair (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1otherpairfitmindist,
                    "Minimum fit-line distance for daughter1 with other nearby particles (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2otherpairfitmindist,
                    "Minimum fit-line distance for daughter2 with other nearby particles (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0otherpairsfitmindistglobal,
                    "Global minimum over daughter-other pair fit distances (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recoprotonchi2ndf,
                    "Reco proton chi2/ndf for K0 daughter1", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recokaonchi2ndf,
                    "Reco kaon chi2/ndf for K0 daughter1", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recopionchi2ndf,
                    "Reco pion chi2/ndf for K0 daughter1", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recoprotonchi2ndf,
                    "Reco proton chi2/ndf for K0 daughter2", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recokaonchi2ndf,
                    "Reco kaon chi2/ndf for K0 daughter2", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recopionchi2ndf,
                    "Reco pion chi2/ndf for K0 daughter2", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0truemomentum, "True momentum magnitude for K0 (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partruestartmom,
                    "True start momentum magnitude for K0 parent (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partrueendmom,
                    "True end momentum magnitude for K0 parent (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1truestartmom,
                    "True start momentum magnitude for K0 daughter1 (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1trueendmom,
                    "True end momentum magnitude for K0 daughter1 (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2truestartmom,
                    "True start momentum magnitude for K0 daughter2 (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2trueendmom,
                    "True end momentum magnitude for K0 daughter2 (GeV/c)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0truepdg, "True PDG of K0", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0partruepdg, "True PDG of K0 parent", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1truepdg, "True PDG of K0 daughter1", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2truepdg, "True PDG of K0 daughter2", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1truelength, "True length K0 daughter1 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau1recolength, "Reco length K0 daughter1 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2truelength, "True length K0 daughter2 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0dau2recolength, "Reco length K0 daughter2 (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0truelength, "True length K0 (cm)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1recohits, "Reco NHits K0 daughter1", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2recohits, "Reco NHits K0 daughter2", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parrecohits, "Reco NHits true parent of K0", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0partruelength, "True length K0 parent (cm)", ntruek0, nmax);
    AddVarMaxSizeVF(output, k0parrecolength, "Reco length K0 parent (cm)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1otherminpairtruepdg,
                    "True PDG of partner for k0dau1otherpairfitmindist", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0vertexminisdaupair,
                    "1 if global min fit-line sep equals dau-dau pair", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parrecozstartgtzend,
                    "1 if K0 parent reco startZ > endZ (Pandora direction in Z)", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau1recozstartgtzend,
                    "1 if K0 daughter1 reco startZ > endZ", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0dau2recozstartgtzend,
                    "1 if K0 daughter2 reco startZ > endZ", ntruek0, nmax);
    AddVarMaxSizeVI(output, k0parispandorabeam,
                    "1 if K0 parent reco is Pandora beam (isPandora)", ntruek0, nmax);
}

//********************************************************************
void neutralKaonTruthTree::FillNeutralKaonTruthVariables(OutputManager& output,
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
                                                         Int_t parentRecoNHits,
                                                         Float_t parentTrueLength,
                                                         Float_t parentRecoLength,
                                                         Int_t daughter1OtherMinSepTruePdg,
                                                         Int_t vertexMinIsK0DaughterPair,
                                                         Int_t parentRecoZStartGreaterThanEnd,
                                                         Int_t daughter1RecoZStartGreaterThanEnd,
                                                         Int_t daughter2RecoZStartGreaterThanEnd,
                                                         Int_t parentIsPandoraBeam){
    output.FillVectorVar(k0parhasrecoobject, parentHasRecoObject ? 1 : 0);
    output.FillVectorVar(k0dau1hasrecoobject, daughter1HasRecoObject ? 1 : 0);
    output.FillVectorVar(k0dau2hasrecoobject, daughter2HasRecoObject ? 1 : 0);
    output.FillVectorVar(k0daurecostartdistance, daughterRecoStartDistance);
    output.FillVectorVar(k0daurecoparentiddiff, daughterRecoParentIdDifference);
    output.FillVectorVar(k0dau1recoparentvstrueparentrecoiddiff, daughter1RecoParentVsTrueParentRecoIdDiff);
    output.FillVectorVar(k0dau2recoparentvstrueparentrecoiddiff, daughter2RecoParentVsTrueParentRecoIdDiff);
    output.FillVectorVar(k0parrecotruestartdist, parentRecoTrueStartDistance);
    output.FillVectorVar(k0parrecotrueenddist, parentRecoTrueEndDistance);
    output.FillVectorVar(k0dau1recotruestartdist, daughter1RecoTrueStartDistance);
    output.FillVectorVar(k0dau1recotrueenddist, daughter1RecoTrueEndDistance);
    output.FillVectorVar(k0dau2recotruestartdist, daughter2RecoTrueStartDistance);
    output.FillVectorVar(k0dau2recotrueenddist, daughter2RecoTrueEndDistance);
    output.FillVectorVar(k0daupairfitmindist, daughterPairFitMinDistance);
    output.FillVectorVar(k0dau1otherpairfitmindist, daughter1OtherPairFitMinDistance);
    output.FillVectorVar(k0dau2otherpairfitmindist, daughter2OtherPairFitMinDistance);
    output.FillVectorVar(k0otherpairsfitmindistglobal, otherPairsFitMinDistanceGlobal);
    output.FillVectorVar(k0dau1recoprotonchi2ndf, daughter1RecoProtonChi2Ndf);
    output.FillVectorVar(k0dau1recokaonchi2ndf, daughter1RecoKaonChi2Ndf);
    output.FillVectorVar(k0dau1recopionchi2ndf, daughter1RecoPionChi2Ndf);
    output.FillVectorVar(k0dau2recoprotonchi2ndf, daughter2RecoProtonChi2Ndf);
    output.FillVectorVar(k0dau2recokaonchi2ndf, daughter2RecoKaonChi2Ndf);
    output.FillVectorVar(k0dau2recopionchi2ndf, daughter2RecoPionChi2Ndf);
    output.FillVectorVar(k0truemomentum, k0TrueMomentum);
    output.FillVectorVar(k0partruestartmom, parentTrueStartMomentum);
    output.FillVectorVar(k0partrueendmom, parentTrueEndMomentum);
    output.FillVectorVar(k0dau1truestartmom, daughter1TrueStartMomentum);
    output.FillVectorVar(k0dau1trueendmom, daughter1TrueEndMomentum);
    output.FillVectorVar(k0dau2truestartmom, daughter2TrueStartMomentum);
    output.FillVectorVar(k0dau2trueendmom, daughter2TrueEndMomentum);
    output.FillVectorVar(k0truepdg, k0TruePdg);
    output.FillVectorVar(k0partruepdg, parentTruePdg);
    output.FillVectorVar(k0dau1truepdg, daughter1TruePdg);
    output.FillVectorVar(k0dau2truepdg, daughter2TruePdg);
    output.FillVectorVar(k0dau1truelength, daughter1TrueLength);
    output.FillVectorVar(k0dau1recolength, daughter1RecoLength);
    output.FillVectorVar(k0dau2truelength, daughter2TrueLength);
    output.FillVectorVar(k0dau2recolength, daughter2RecoLength);
    output.FillVectorVar(k0truelength, k0TrueLength);
    output.FillVectorVar(k0dau1recohits, daughter1RecoNHits);
    output.FillVectorVar(k0dau2recohits, daughter2RecoNHits);
    output.FillVectorVar(k0parrecohits, parentRecoNHits);
    output.FillVectorVar(k0partruelength, parentTrueLength);
    output.FillVectorVar(k0parrecolength, parentRecoLength);
    output.FillVectorVar(k0dau1otherminpairtruepdg, daughter1OtherMinSepTruePdg);
    output.FillVectorVar(k0vertexminisdaupair, vertexMinIsK0DaughterPair);
    output.FillVectorVar(k0parrecozstartgtzend, parentRecoZStartGreaterThanEnd);
    output.FillVectorVar(k0dau1recozstartgtzend, daughter1RecoZStartGreaterThanEnd);
    output.FillVectorVar(k0dau2recozstartgtzend, daughter2RecoZStartGreaterThanEnd);
    output.FillVectorVar(k0parispandorabeam, parentIsPandoraBeam);
}