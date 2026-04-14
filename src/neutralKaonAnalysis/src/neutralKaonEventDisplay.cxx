#include "neutralKaonEventDisplay.hxx"
#include "ToyBoxNeutralKaon.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdCreationUtils.hxx"
#include "OutputManager.hxx"
#include <TEveGeoShape.h>
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveSelection.h>
#include <TGeoSphere.h>
#include <TGeoVolume.h>
#include "Parameters.hxx"
#include <TEveElement.h>
#include <TEveLine.h>
#include <TEveScene.h>
#include <TEveTrans.h>
#include <TEveText.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TTimer.h>
#include <TH1F.h>
#include <TVector3.h>
#include <iostream>
#include <cmath>
#include <array>
#include <set>
#include <string>
#include <unordered_map>

ClassImp(neutralKaonEventDisplay);

namespace {

std::string ProcessEnumToString(Int_t process){
    if (process < 0) return "";
    static AnaTrueParticlePD dummy;
    return dummy.ConvertProcess(static_cast<AnaTrueParticleB::ProcessEnum>(process));
}

TEveElementList* PrepareGroup(TEveScene* scene, const char* name) {
    TEveElementList* group = new TEveElementList(name);
    if (scene) scene->AddElement(group);
    return group;
}

void AppendElementToGroup(TEveScene* /*scene*/, TEveElementList* group, TEveElement* element) {
    if (group && element) group->AddElement(element);
}

TEveElement* FindElementByNameRecursive(TEveElement* parent, const char* name) {
    if (!parent || !name) return nullptr;

    const char* parentName = parent->GetElementName();
    if (parentName && std::string(parentName) == name) {
        return parent;
    }

    for (TEveElement::List_ci it = parent->BeginChildren(); it != parent->EndChildren(); ++it) {
        TEveElement* child = *it;
        if (!child) continue;
        if (TEveElement* found = FindElementByNameRecursive(child, name)) {
            return found;
        }
    }

    return nullptr;
}

class GroupCheckboxSyncTimer : public TTimer {
public:
    explicit GroupCheckboxSyncTimer(neutralKaonEventDisplay* owner)
        : TTimer(250, kTRUE), _owner(owner) {}

    Bool_t Notify() override {
        if (_owner) {
            _owner->SyncGroupCheckboxBehavior();
        }
        Reset();
        return kTRUE;
    }

private:
    neutralKaonEventDisplay* _owner;
};

void HideAllChildrenInGroup(TEveScene* scene, const char* groupName) {
    if (!scene || !groupName) return;

    TEveElement* found = FindElementByNameRecursive(scene, groupName);
    TEveElementList* group = dynamic_cast<TEveElementList*>(found);
    if (!group) return;

    for (TEveElement::List_i it = group->BeginChildren(); it != group->EndChildren(); ++it) {
        TEveElement* child = *it;
        if (!child) continue;
        child->SetRnrSelf(kFALSE);
        child->SetRnrChildren(kFALSE);
        child->ElementChanged(kTRUE, kTRUE);
    }

    if (gEve) {
        gEve->Redraw3D(kTRUE);
    }
}

void BuildParentTrajectoryHistogram(const AnaParticlePD* parent, Float_t* outHist) {
    if (!outHist) return;
    std::fill(outHist, outHist + (neutralKaonEventDisplay::kParentTrajHistBins * 3), 0.f);
    if (!parent) return;

    const auto& points = parent->TrjPoints;
    if (points.empty()) return;

    const Float_t binWidth = (neutralKaonEventDisplay::kParentTrajHistMax -
                              neutralKaonEventDisplay::kParentTrajHistMin) /
                             neutralKaonEventDisplay::kParentTrajHistBins;
    std::array<Int_t, neutralKaonEventDisplay::kParentTrajHistBins> countsX{};
    std::array<Int_t, neutralKaonEventDisplay::kParentTrajHistBins> countsY{};
    std::array<Int_t, neutralKaonEventDisplay::kParentTrajHistBins> countsZ{};

    auto fillBin = [&](Float_t value, std::array<Int_t, neutralKaonEventDisplay::kParentTrajHistBins>& counts) {
        if (!std::isfinite(value)) return;
        Float_t clamped = std::max(neutralKaonEventDisplay::kParentTrajHistMin,
                                   std::min(neutralKaonEventDisplay::kParentTrajHistMax - 1e-6f, value));
        Int_t bin = static_cast<Int_t>((clamped - neutralKaonEventDisplay::kParentTrajHistMin) / binWidth);
        bin = std::max(0, std::min(bin, neutralKaonEventDisplay::kParentTrajHistBins - 1));
        counts[bin]++;
    };

    for (const auto& point : points) {
        TVector3 dir = point.Direction;
        if (!std::isfinite(dir.X()) || !std::isfinite(dir.Y()) || !std::isfinite(dir.Z())) {
            continue;
        }
        if (dir.Mag() <= 0) continue;
        dir = dir.Unit();

        fillBin(dir.X(), countsX);
        fillBin(dir.Y(), countsY);
        fillBin(dir.Z(), countsZ);
    }

    for (Int_t b = 0; b < neutralKaonEventDisplay::kParentTrajHistBins; ++b) {
        outHist[b] = static_cast<Float_t>(countsX[b]);
        outHist[neutralKaonEventDisplay::kParentTrajHistBins + b] = static_cast<Float_t>(countsY[b]);
        outHist[2 * neutralKaonEventDisplay::kParentTrajHistBins + b] = static_cast<Float_t>(countsZ[b]);
    }
}

// True K0 decay vertex from truth daughters: prefer two charged pions, else first valid daughters.
bool TrueK0McDecayVertex(Float_t out[3], Int_t nDaughters, const Int_t* daughterPDG, const Float_t* daughterStarts) {
    const auto valid = [](Float_t v) { return v > -900.f; };
    if (!daughterStarts || nDaughters < 1) return false;
    if (daughterPDG) {
        Int_t piIdx[2];
        Int_t nPi = 0;
        constexpr Int_t kMaxDauSlots = 8; // matches neutralKaonEventDisplay::kMaxTrueDaughters
        for (Int_t d = 0; d < nDaughters && d < kMaxDauSlots && nPi < 2; ++d) {
            if (std::abs(daughterPDG[d]) != 211) continue;
            const Float_t* p = daughterStarts + 3 * d;
            if (valid(p[0]) && valid(p[1]) && valid(p[2])) piIdx[nPi++] = d;
        }
        if (nPi == 2) {
            const Float_t* p0 = daughterStarts + 3 * piIdx[0];
            const Float_t* p1 = daughterStarts + 3 * piIdx[1];
            out[0] = 0.5f * (p0[0] + p1[0]);
            out[1] = 0.5f * (p0[1] + p1[1]);
            out[2] = 0.5f * (p0[2] + p1[2]);
            return true;
        }
    }
    if (nDaughters >= 2) {
        const Float_t* p0 = daughterStarts;
        const Float_t* p1 = daughterStarts + 3;
        if (valid(p0[0]) && valid(p0[1]) && valid(p0[2]) && valid(p1[0]) && valid(p1[1]) && valid(p1[2])) {
            out[0] = 0.5f * (p0[0] + p1[0]);
            out[1] = 0.5f * (p0[1] + p1[1]);
            out[2] = 0.5f * (p0[2] + p1[2]);
            return true;
        }
    }
    const Float_t* p0 = daughterStarts;
    if (valid(p0[0]) && valid(p0[1]) && valid(p0[2])) {
        out[0] = p0[0];
        out[1] = p0[1];
        out[2] = p0[2];
        return true;
    }
    return false;
}

}

//********************************************************************
neutralKaonEventDisplay::neutralKaonEventDisplay() : pdEventDisplay() {
//********************************************************************
    // Set class name for auto-detection
    _eventDisplayClassName = "neutralKaonEventDisplay";

    _nK0Candidates = 0;
    _nTrueK0 = 0;
    _nAllTrueParticles = 0;

    // Initialize arrays
    for (Int_t i = 0; i < kMaxK0; i++) {
        _k0_daughter1ID[i] = -1;
        _k0_daughter2ID[i] = -1;
        _k0_parentID[i] = -1;
        _k0_secondParticleID[i] = -1;
        _k0_creationVtxRadius[i] = 0;
        _k0_annihilationVtxRadius[i] = 0;
        _k0_creationVtxDeg[i] = -999;
        _k0_annihilationVtxDeg[i] = -999;
        _k0_fitLineLength[i] = 0;
        _k0_hasTrueObject[i] = 0;
        _k0_truePDG[i] = -999;
        _k0_trueProcessEnd[i] = -1;
        _k0_trueParentPDG[i] = 0;
        _k0_trueNDaughters[i] = 0;
        _k0_trueNSiblings[i] = 0;
        _k0_trueParentNDaughters[i] = 0;
        for (Int_t j = 0; j < 5; j++) {
            _k0_creationVtxDegDist[i][j] = -999;
            _k0_annihilationVtxDegDist[i][j] = -999;
        }
        for (Int_t j = 0; j < 3; j++) {
            _k0_creationVtxPos[i][j] = -999;
            _k0_annihilationVtxPos[i][j] = -999;
            _k0_annihilationVtxFitPos[i][j] = -999;
            _k0_startPos[i][j] = -999;
            _k0_endPos[i][j] = -999;
            _k0_annVtx_fitLine1Start[i][j] = -999;
            _k0_annVtx_fitLine1Dir[i][j] = -999;
            _k0_annVtx_fitLine2Start[i][j] = -999;
            _k0_annVtx_fitLine2Dir[i][j] = -999;
            _k0_annVtx_pandoraLine1Dir[i][j] = -999;
            _k0_annVtx_pandoraLine2Dir[i][j] = -999;
            _k0_annVtx_closestPt1[i][j] = -999;
            _k0_annVtx_closestPt2[i][j] = -999;
            _k0_annVtx_pandoraClosestPt1[i][j] = -999;
            _k0_annVtx_pandoraClosestPt2[i][j] = -999;
            _k0_annVtx_momentumDirFit[i][j] = -999;
            _k0_trueK0Dir[i][j] = -999;
            _k0_trueDecayVtxFromRecoDaughters[i][j] = -999;
            _k0_creationVtx_fitLineBeamStart[i][j] = -999;
            _k0_creationVtx_fitLineBeamDir[i][j] = -999;
            _k0_creationVtx_fitLineSecondStart[i][j] = -999;
            _k0_creationVtx_fitLineSecondDir[i][j] = -999;
            _k0_creationVtx_closestPtBeam[i][j] = -999;
            _k0_creationVtx_closestPtSecond[i][j] = -999;
            _k0_trueStartPos[i][j] = -999;
            _k0_trueEndPos[i][j] = -999;
            _k0_trueParentStartPos[i][j] = -999;
            _k0_trueParentEndPos[i][j] = -999;
            _k0_trueParentDaughterStartPos[i][j] = -999;
            _k0_trueParentDaughterEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueDaughters; ++j) {
            _k0_trueParentDaughterPDG[i][j] = 0;
        }
        for (Int_t j = 0; j < kMaxTrueDaughters*3; ++j) {
            _k0_trueDaughterStartPos[i][j] = -999;
            _k0_trueDaughterEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueDaughters; ++j) {
            _k0_trueDaughterPDG[i][j] = 0;
        }
        for (Int_t j = 0; j < kMaxTrueSiblings*3; ++j) {
            _k0_trueSiblingStartPos[i][j] = -999;
            _k0_trueSiblingEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueSiblings; ++j) {
            _k0_trueSiblingPDG[i][j] = 0;
        }
        for (Int_t j = 0; j < 3; ++j) {
            _k0_parentStartPos[i][j] = -999;
            _k0_parentEndPos[i][j] = -999;
            _k0_parentEndPosCorrected[i][j] = -999;
            _k0_parentTailFitAnchor[i][j] = -999;
            _k0_parentTailFitDir[i][j] = -999;
            _k0_parentTrajDir[i][j] = -999;
            _k0_secondTrajDir[i][j] = -999;
            _k0_dau1TrajDir[i][j] = -999;
            _k0_dau2TrajDir[i][j] = -999;
        }
        _k0_parentTailFitLength[i] = -999.f;
        std::fill_n(_k0_parentTailHitsRaw[i], kMaxBeamTailHits * 3, -999.f);
        std::fill_n(_k0_parentTailHitsProjected[i], kMaxBeamTailHits * 3, -999.f);
        _k0_parentTailHitsRawN[i] = 0;
        _k0_parentTailHitsProjectedN[i] = 0;
        _k0_parentLength[i] = -999.f;
        std::fill_n(_k0_parentTrajDirHist[i], kParentTrajHistBins * 3, 0.f);
        _k0_parentTrajDirNPts[i] = 0;
        _k0_secondTrajDirNPts[i] = 0;
        _k0_dau1TrajDirNPts[i] = 0;
        _k0_dau2TrajDirNPts[i] = 0;
    }

    for (Int_t i = 0; i < kMaxTrueK0; ++i) {
        _trueK0_PDG[i] = 0;
        _trueK0_processEnd[i] = -1;
        _trueK0_parentPDG[i] = 0;
        _trueK0_nDaughters[i] = 0;
        _trueK0_nSiblings[i] = 0;
        for (Int_t j = 0; j < 3; ++j) {
            _trueK0_startPos[i][j] = -999;
            _trueK0_endPos[i][j] = -999;
            _trueK0_parentStartPos[i][j] = -999;
            _trueK0_parentEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueDaughters*3; ++j) {
            _trueK0_daughterStartPos[i][j] = -999;
            _trueK0_daughterEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueDaughters; ++j) {
            _trueK0_daughterPDG[i][j] = 0;
        }
        for (Int_t j = 0; j < kMaxTrueSiblings*3; ++j) {
            _trueK0_siblingStartPos[i][j] = -999;
            _trueK0_siblingEndPos[i][j] = -999;
        }
        for (Int_t j = 0; j < kMaxTrueSiblings; ++j) {
            _trueK0_siblingPDG[i][j] = 0;
        }
    }

    for (Int_t i = 0; i < kMaxAllTrueParticles; ++i) {
        _allTrueParticle_PDG[i] = 0;
        _allTrueParticle_processEnd[i] = -1;
        for (Int_t j = 0; j < 3; ++j) {
            _allTrueParticle_startPos[i][j] = -999;
            _allTrueParticle_endPos[i][j] = -999;
        }
    }
}

void neutralKaonEventDisplay::EnsureSelectionHooks() {
    if (_parentDirSelectionHooked || !gEve) {
        return;
    }

    TEveSelection* selection = gEve->GetSelection();
    if (!selection) {
        return;
    }

    selection->Connect("SelectionAdded(TEveElement*)", "neutralKaonEventDisplay", this,
                       "OnSelectionAdded(TEveElement*)");
    selection->Connect("SelectionRepeated(TEveElement*)", "neutralKaonEventDisplay", this,
                       "OnSelectionRepeated(TEveElement*)");
    _parentDirSelectionHooked = kTRUE;
}

void neutralKaonEventDisplay::EnsureGroupCheckboxSync() {
    if (_groupCheckboxSyncTimer) {
        return;
    }
    _groupCheckboxSyncTimer = new GroupCheckboxSyncTimer(this);
    _groupCheckboxSyncTimer->TurnOn();
}

Bool_t neutralKaonEventDisplay::SyncGroupCheckboxBehavior() {
    if (!_scene3D) {
        return kTRUE;
    }

    auto syncOneGroup = [&](const char* groupName, Bool_t& lastChecked) {
        TEveElement* found = FindElementByNameRecursive(_scene3D, groupName);
        TEveElementList* group = dynamic_cast<TEveElementList*>(found);
        if (!group) return;

        if (!group->GetRnrChildren()) {
            group->SetRnrChildren(kTRUE);
            group->ElementChanged(kTRUE, kTRUE);
        }

        const Bool_t nowChecked = group->GetRnrSelf();
        if (_groupCheckboxStateInitialized && lastChecked && !nowChecked) {
            HideAllChildrenInGroup(_scene3D, groupName);
        }
        lastChecked = nowChecked;
    };

    syncOneGroup("Neutral Particles", _lastNeutralGroupChecked);
    syncOneGroup("True Particles", _lastTrueGroupChecked);
    _groupCheckboxStateInitialized = kTRUE;

    return kTRUE;
}

Bool_t neutralKaonEventDisplay::IsGroupVisible(const char* groupName) const {
    if (!_scene3D || !groupName || groupName[0] == '\0') {
        return kTRUE;
    }

    TEveElement* found = FindElementByNameRecursive(_scene3D, groupName);
    if (!found) {
        return kTRUE;
    }

    return found->GetRnrSelf();
}

void neutralKaonEventDisplay::OnSelectionAdded(TEveElement* element) {
    if (!element) return;

    auto it = _parentDirElementToIndex.find(element);
    if (it == _parentDirElementToIndex.end()) {
        return;
    }
    ShowParentTrajectoryHistograms(it->second);
}

void neutralKaonEventDisplay::OnSelectionRepeated(TEveElement* element) {
    OnSelectionAdded(element);
}

void neutralKaonEventDisplay::ClearParentDirectionVisuals() {
    for (Int_t comp = 0; comp < 3; ++comp) {
        if (_parentDirTrueLines[comp]) {
            delete _parentDirTrueLines[comp];
            _parentDirTrueLines[comp] = nullptr;
        }
        if (_parentDirLegends[comp]) {
            delete _parentDirLegends[comp];
            _parentDirLegends[comp] = nullptr;
        }
    }
}

void neutralKaonEventDisplay::ShowParentTrajectoryHistograms(Int_t index) {
    if (index < 0 || index >= _nK0Candidates) {
        return;
    }
    if (_k0_parentTrajDirNPts[index] <= 0) {
        std::cout << "No trajectory direction samples for parent of candidate " << index << std::endl;
        return;
    }

    Bool_t hasTrueDirVec = kFALSE;
    TVector3 trueDirVec;
    const Float_t startX = _k0_trueParentStartPos[index][0];
    const Float_t startY = _k0_trueParentStartPos[index][1];
    const Float_t startZ = _k0_trueParentStartPos[index][2];
    const Float_t endX = _k0_trueParentEndPos[index][0];
    const Float_t endY = _k0_trueParentEndPos[index][1];
    const Float_t endZ = _k0_trueParentEndPos[index][2];
    if (startX > -900.f && startY > -900.f && startZ > -900.f &&
        endX > -900.f && endY > -900.f && endZ > -900.f) {
        TVector3 delta(endX - startX, endY - startY, endZ - startZ);
        if (delta.Mag() > 0) {
            trueDirVec = delta.Unit();
            hasTrueDirVec = kTRUE;
        }
    }

    if (!_parentDirCanvas) {
        _parentDirCanvas = new TCanvas("ParentTrajectoryDirection", "Parent Trajectory Direction", 1200, 400);
        _parentDirCanvas->Divide(3, 1);
    }

    ClearParentDirectionVisuals();

    static const char axisLabels[3] = {'X', 'Y', 'Z'};
    for (Int_t comp = 0; comp < 3; ++comp) {
        _parentDirCanvas->cd(comp + 1);
        if (!_parentDirHists[comp]) {
            _parentDirHists[comp] = new TH1F(Form("ParentTrajDir_%c", axisLabels[comp]),
                                             Form("Parent Direction %c Component", axisLabels[comp]),
                                             kParentTrajHistBins,
                                             kParentTrajHistMin,
                                             kParentTrajHistMax);
            _parentDirHists[comp]->SetStats(kFALSE);
            _parentDirHists[comp]->SetLineWidth(2);
        }

        TH1F* hist = _parentDirHists[comp];
        hist->Reset("ICES");
        for (Int_t bin = 0; bin < kParentTrajHistBins; ++bin) {
            hist->SetBinContent(bin + 1, _k0_parentTrajDirHist[index][comp * kParentTrajHistBins + bin]);
        }

        Int_t color = (comp == 0) ? kAzure + 2 : (comp == 1 ? kGreen + 2 : kOrange + 7);
        hist->SetLineColor(color);
        hist->SetTitle(Form("Parent Direction %c Component", axisLabels[comp]));
        hist->GetXaxis()->SetTitle("Direction component value");
        hist->GetYaxis()->SetTitle("Entries");
        hist->Draw("hist");

        Double_t mpv = _k0_parentTrajDir[index][comp];
        Double_t trueVal = hasTrueDirVec ? trueDirVec[comp] : -999.f;
        Double_t maxVal = hist->GetMaximum();
        if (maxVal <= 0) maxVal = 1.0;

        if (_parentDirTrueLines[comp]) {
            delete _parentDirTrueLines[comp];
            _parentDirTrueLines[comp] = nullptr;
        }
        if (trueVal > -900.f) {
            _parentDirTrueLines[comp] = new TLine(trueVal, 0, trueVal, maxVal * 1.05);
            _parentDirTrueLines[comp]->SetLineColor(kRed + 1);
            _parentDirTrueLines[comp]->SetLineStyle(2);
            _parentDirTrueLines[comp]->SetLineWidth(2);
            _parentDirTrueLines[comp]->Draw();
        }

        if (_parentDirLegends[comp]) {
            delete _parentDirLegends[comp];
            _parentDirLegends[comp] = nullptr;
        }
        _parentDirLegends[comp] = new TLegend(0.15, 0.75, 0.85, 0.9);
        _parentDirLegends[comp]->AddEntry(hist, Form("MPV = %.3f", mpv), "l");
        if (_parentDirTrueLines[comp]) {
            _parentDirLegends[comp]->AddEntry(_parentDirTrueLines[comp], Form("True = %.3f", trueVal), "l");
        }
        _parentDirLegends[comp]->SetBorderSize(0);
        _parentDirLegends[comp]->SetFillStyle(0);
        _parentDirLegends[comp]->Draw();
    }

    _parentDirCanvas->Update();
}

//********************************************************************
neutralKaonEventDisplay::~neutralKaonEventDisplay() {
//********************************************************************
    ClearParentDirectionVisuals();
    for (auto& hist : _parentDirHists) {
        if (hist) {
            delete hist;
            hist = nullptr;
        }
    }
    if (_parentDirCanvas) {
        delete _parentDirCanvas;
        _parentDirCanvas = nullptr;
    }
    if (_groupCheckboxSyncTimer) {
        _groupCheckboxSyncTimer->TurnOff();
        delete _groupCheckboxSyncTimer;
        _groupCheckboxSyncTimer = nullptr;
    }
}

//********************************************************************
void neutralKaonEventDisplay::AddAnalysisVariables(OutputManager& output, Int_t tree_index) {
//********************************************************************
    // Add K0 candidate variables
    output.AddVectorVar(tree_index, edk0_daughter1ID, "ED_k0_daughter1ID", "I", "K0 daughter 1 unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_daughter2ID, "ED_k0_daughter2ID", "I", "K0 daughter 2 unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_parentID, "ED_k0_parentID", "I", "K0 parent unique ID", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_secondParticleID, "ED_k0_secondParticleID", "I", "K0 second particle unique ID (creation vtx)", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_creationVtxRadius, "ED_k0_creationVtxRadius", "F", "K0 creation vertex radius", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_annihilationVtxRadius, "ED_k0_annihilationVtxRadius", "F", "K0 annihilation vertex radius", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_creationVtxDeg, "ED_k0_creationVtxDeg", "I", "K0 creation vertex degeneracy", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_annihilationVtxDeg, "ED_k0_annihilationVtxDeg", "I", "K0 annihilation vertex degeneracy", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_creationVtxPos, "ED_k0_creationVtxPos", "F", "K0 creation vertex positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annihilationVtxPos, "ED_k0_annihilationVtxPos", "F", "K0 annihilation vertex positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annihilationVtxFitPos, "ED_k0_annihilationVtxFitPos", "F", "K0 annihilation vertex fit positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_startPos, "ED_k0_startPos", "F", "K0 start positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_endPos, "ED_k0_endPos", "F", "K0 end positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);

    // Fitted lines for vertex finding visualization
    output.AddMatrixVar(tree_index, edk0_annVtx_fitLine1Start, "ED_k0_annVtx_fitLine1Start", "F", "Annihilation vtx daughter1 fit line start", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_fitLine1Dir, "ED_k0_annVtx_fitLine1Dir", "F", "Annihilation vtx daughter1 fit line direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_fitLine2Start, "ED_k0_annVtx_fitLine2Start", "F", "Annihilation vtx daughter2 fit line start", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_fitLine2Dir, "ED_k0_annVtx_fitLine2Dir", "F", "Annihilation vtx daughter2 fit line direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_pandoraLine1Dir, "ED_k0_annVtx_pandoraLine1Dir", "F", "Annihilation vtx daughter1 Pandora direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_pandoraLine2Dir, "ED_k0_annVtx_pandoraLine2Dir", "F", "Annihilation vtx daughter2 Pandora direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_closestPt1, "ED_k0_annVtx_closestPt1", "F", "Closest point on daughter1 fit line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_closestPt2, "ED_k0_annVtx_closestPt2", "F", "Closest point on daughter2 fit line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_pandoraClosestPt1, "ED_k0_annVtx_pandoraClosestPt1", "F", "Closest point on daughter1 Pandora line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_pandoraClosestPt2, "ED_k0_annVtx_pandoraClosestPt2", "F", "Closest point on daughter2 Pandora line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_annVtx_momentumDirFit, "ED_k0_annVtx_momentumDirFit", "F", "Annihilation-vertex momentum direction from fit daughters", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_trueK0Dir, "ED_k0_trueK0Dir", "F", "Associated true K0 direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_trueDecayVtxFromRecoDaughters, "ED_k0_trueDecayVtxFromRecoDaughters", "F",
                        "True decay vertex: midpoint of reco daughters' true start positions", ednK0Candidates,
                        "ED_nK0Candidates", -kMaxK0, 3);

    // Fitted lines for creation vertex
    output.AddMatrixVar(tree_index, edk0_creationVtx_fitLineBeamStart, "ED_k0_creationVtx_fitLineBeamStart", "F", "Creation vtx beam fit line start", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_creationVtx_fitLineBeamDir, "ED_k0_creationVtx_fitLineBeamDir", "F", "Creation vtx beam fit line direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_creationVtx_fitLineSecondStart, "ED_k0_creationVtx_fitLineSecondStart", "F", "Creation vtx second fit line start", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_creationVtx_fitLineSecondDir, "ED_k0_creationVtx_fitLineSecondDir", "F", "Creation vtx second fit line direction", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_creationVtx_closestPtBeam, "ED_k0_creationVtx_closestPtBeam", "F", "Closest point on beam fit line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_creationVtx_closestPtSecond, "ED_k0_creationVtx_closestPtSecond", "F", "Closest point on second fit line", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);

    output.AddVectorVar(tree_index, edk0_fitLineLength, "ED_k0_fitLineLength", "F", "Fit line length used for vertex finding", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);

    // Degeneracy distances
    output.AddMatrixVar(tree_index, edk0_creationVtxDegDist, "ED_k0_creationVtxDegDist", "F", "Creation vertex degeneracy distances", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 5);
    output.AddMatrixVar(tree_index, edk0_annihilationVtxDegDist, "ED_k0_annihilationVtxDegDist", "F", "Annihilation vertex degeneracy distances", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 5);

    // True K0 trajectory
    output.AddVectorVar(tree_index, edk0_hasTrueObject, "ED_k0_hasTrueObject", "I", "K0 has true object", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_trueStartPos, "ED_k0_trueStartPos", "F", "True K0 start position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_trueEndPos, "ED_k0_trueEndPos", "F", "True K0 end position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_truePDG, "ED_k0_truePDG", "I", "True K0 PDG", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_trueProcessEnd, "ED_k0_trueProcessEnd", "I", "True K0 end process enum", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_trueParentPDG, "ED_k0_trueParentPDG", "I", "True K0 parent PDG", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_trueParentStartPos, "ED_k0_trueParentStartPos", "F", "True K0 parent start position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_trueParentEndPos, "ED_k0_trueParentEndPos", "F", "True K0 parent end position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_trueParentNDaughters, "ED_k0_trueParentNDaughters", "I", "Number of true daughters of K0 parent", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_trueParentDaughterStartPos, "ED_k0_trueParentDaughterStartPos", "F", "True K0 parent daughters start positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edk0_trueParentDaughterEndPos, "ED_k0_trueParentDaughterEndPos", "F", "True K0 parent daughters end positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edk0_trueParentDaughterPDG, "ED_k0_trueParentDaughterPDG", "I", "True K0 parent daughters PDG codes", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters);
    output.AddMatrixVar(tree_index, edk0_parentStartPos, "ED_k0_parentStartPos", "F", "Reco K0 parent start position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_parentEndPos, "ED_k0_parentEndPos", "F", "Reconstructed K0 parent end position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_parentEndPosCorrected, "ED_k0_parentEndPosCorrected", "F", "Projected parent end position", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_parentTailFitAnchor, "ED_k0_parentTailFitAnchor", "F", "Beam backward-fit anchor used for tail projection", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_parentTailFitDir, "ED_k0_parentTailFitDir", "F", "Beam backward-fit direction used for tail projection", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_parentTailFitLength, "ED_k0_parentTailFitLength", "F", "Beam backward-fit length used for tail projection", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_parentTailHitsRaw, "ED_k0_parentTailHitsRaw", "F", "Raw parent tail trajectory points", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxBeamTailHits * 3);
    output.AddMatrixVar(tree_index, edk0_parentTailHitsProjected, "ED_k0_parentTailHitsProjected", "F", "Projected parent tail trajectory points", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxBeamTailHits * 3);
    output.AddVectorVar(tree_index, edk0_parentTailHitsRawN, "ED_k0_parentTailHitsRawN", "I", "Raw parent tail hit count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_parentTailHitsProjectedN, "ED_k0_parentTailHitsProjectedN", "I", "Projected parent tail hit count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_parentLength, "ED_k0_parentLength", "F", "Reconstructed K0 parent length", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_parentTrajDir, "ED_k0_parentTrajDir", "F", "K0 parent trajectory direction MPV", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddMatrixVar(tree_index, edk0_parentTrajDirHist, "ED_k0_parentTrajDirHist", "F", "K0 parent trajectory direction histograms (XYZ, 60 bins each)", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kParentTrajHistBins * 3);
    output.AddVectorVar(tree_index, edk0_parentTrajDirNPts, "ED_k0_parentTrajDirNPts", "I", "K0 parent trajectory direction sample count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_secondTrajDir, "ED_k0_secondTrajDir", "F", "Creation partner trajectory direction MPV", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_secondTrajDirNPts, "ED_k0_secondTrajDirNPts", "I", "Creation partner trajectory direction sample count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_dau1TrajDir, "ED_k0_dau1TrajDir", "F", "Daughter1 trajectory direction MPV", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_dau1TrajDirNPts, "ED_k0_dau1TrajDirNPts", "I", "Daughter1 trajectory direction sample count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_dau2TrajDir, "ED_k0_dau2TrajDir", "F", "Daughter2 trajectory direction MPV", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, 3);
    output.AddVectorVar(tree_index, edk0_dau2TrajDirNPts, "ED_k0_dau2TrajDirNPts", "I", "Daughter2 trajectory direction sample count", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddVectorVar(tree_index, edk0_trueNDaughters, "ED_k0_trueNDaughters", "I", "Number of true daughters for K0", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_trueDaughterStartPos, "ED_k0_trueDaughterStartPos", "F", "True K0 daughters start positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edk0_trueDaughterEndPos, "ED_k0_trueDaughterEndPos", "F", "True K0 daughters end positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edk0_trueDaughterPDG, "ED_k0_trueDaughterPDG", "I", "True K0 daughters PDG codes", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueDaughters);
    output.AddVectorVar(tree_index, edk0_trueNSiblings, "ED_k0_trueNSiblings", "I", "Number of true siblings for K0", ednK0Candidates, "ED_nK0Candidates", -kMaxK0);
    output.AddMatrixVar(tree_index, edk0_trueSiblingStartPos, "ED_k0_trueSiblingStartPos", "F", "True K0 siblings start positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueSiblings*3);
    output.AddMatrixVar(tree_index, edk0_trueSiblingEndPos, "ED_k0_trueSiblingEndPos", "F", "True K0 siblings end positions", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueSiblings*3);
    output.AddMatrixVar(tree_index, edk0_trueSiblingPDG, "ED_k0_trueSiblingPDG", "I", "True K0 siblings PDG codes", ednK0Candidates, "ED_nK0Candidates", -kMaxK0, kMaxTrueSiblings);

    // Standalone true K0 information
    output.AddVar(tree_index, ednTrueK0, "ED_nTrueK0", "I", "Number of standalone true K0s");
    output.AddMatrixVar(tree_index, edtrueK0_startPos, "ED_trueK0_startPos", "F", "Standalone true K0 start positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, 3);
    output.AddMatrixVar(tree_index, edtrueK0_endPos, "ED_trueK0_endPos", "F", "Standalone true K0 end positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, 3);
    output.AddVectorVar(tree_index, edtrueK0_PDG, "ED_trueK0_PDG", "I", "Standalone true K0 PDG codes", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0);
    output.AddVectorVar(tree_index, edtrueK0_processEnd, "ED_trueK0_processEnd", "I", "Standalone true K0 end process enum", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0);
    output.AddVectorVar(tree_index, edtrueK0_parentPDG, "ED_trueK0_parentPDG", "I", "Standalone true K0 parent PDG", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0);
    output.AddMatrixVar(tree_index, edtrueK0_parentStartPos, "ED_trueK0_parentStartPos", "F", "Standalone true K0 parent start position", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, 3);
    output.AddMatrixVar(tree_index, edtrueK0_parentEndPos, "ED_trueK0_parentEndPos", "F", "Standalone true K0 parent end position", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, 3);
    output.AddVectorVar(tree_index, edtrueK0_nDaughters, "ED_trueK0_nDaughters", "I", "Number of standalone true K0 daughters", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0);
    output.AddMatrixVar(tree_index, edtrueK0_daughterStartPos, "ED_trueK0_daughterStartPos", "F", "Standalone true K0 daughters start positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edtrueK0_daughterEndPos, "ED_trueK0_daughterEndPos", "F", "Standalone true K0 daughters end positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueDaughters*3);
    output.AddMatrixVar(tree_index, edtrueK0_daughterPDG, "ED_trueK0_daughterPDG", "I", "Standalone true K0 daughters PDG codes", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueDaughters);
    output.AddVectorVar(tree_index, edtrueK0_nSiblings, "ED_trueK0_nSiblings", "I", "Number of standalone true K0 siblings", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0);
    output.AddMatrixVar(tree_index, edtrueK0_siblingStartPos, "ED_trueK0_siblingStartPos", "F", "Standalone true K0 siblings start positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueSiblings*3);
    output.AddMatrixVar(tree_index, edtrueK0_siblingEndPos, "ED_trueK0_siblingEndPos", "F", "Standalone true K0 siblings end positions", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueSiblings*3);
    output.AddMatrixVar(tree_index, edtrueK0_siblingPDG, "ED_trueK0_siblingPDG", "I", "Standalone true K0 siblings PDG codes", ednTrueK0, "ED_nTrueK0", -kMaxTrueK0, kMaxTrueSiblings);

    output.AddVar(tree_index, ednAllTrueParticles, "ED_nAllTrueParticles", "I", "Number of all true particles");
    output.AddMatrixVar(tree_index, edallTrueParticle_startPos, "ED_allTrueParticle_startPos", "F", "All true particle start positions", ednAllTrueParticles, "ED_nAllTrueParticles", -kMaxAllTrueParticles, 3);
    output.AddMatrixVar(tree_index, edallTrueParticle_endPos, "ED_allTrueParticle_endPos", "F", "All true particle end positions", ednAllTrueParticles, "ED_nAllTrueParticles", -kMaxAllTrueParticles, 3);
    output.AddVectorVar(tree_index, edallTrueParticle_PDG, "ED_allTrueParticle_PDG", "I", "All true particle PDG codes", ednAllTrueParticles, "ED_nAllTrueParticles", -kMaxAllTrueParticles);
    output.AddVectorVar(tree_index, edallTrueParticle_processEnd, "ED_allTrueParticle_processEnd", "I", "All true particle end process enums", ednAllTrueParticles, "ED_nAllTrueParticles", -kMaxAllTrueParticles);

    std::cout << "neutralKaonEventDisplay::AddAnalysisVariables() - K0 variables added" << std::endl;
}

//********************************************************************
void neutralKaonEventDisplay::FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) {
//********************************************************************
    // Cast box to ToyBoxNeutralKaon
    const ToyBoxNeutralKaon* k0Box = static_cast<const ToyBoxNeutralKaon*>(box);

    std::set<const AnaTrueParticlePD*> matchedTrueKaons;
    std::unordered_map<Int_t, AnaTrueParticlePD*> trueParticleByID;
    for (int i = 0; i < event.nTrueParticles; ++i) {
        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(event.TrueParticles[i]);
        if (!truePart) continue;
        trueParticleByID[truePart->ID] = truePart;
    }

    Int_t allTrueParticlesStored = 0;
    for (Int_t i = 0; i < event.nTrueParticles && allTrueParticlesStored < kMaxAllTrueParticles; ++i) {
        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(event.TrueParticles[i]);
        if (!truePart) continue;

        Float_t startPos[3] = {
            static_cast<Float_t>(truePart->Position[0]),
            static_cast<Float_t>(truePart->Position[1]),
            static_cast<Float_t>(truePart->Position[2])
        };
        Float_t endPos[3] = {
            static_cast<Float_t>(truePart->PositionEnd[0]),
            static_cast<Float_t>(truePart->PositionEnd[1]),
            static_cast<Float_t>(truePart->PositionEnd[2])
        };

        output.FillMatrixVarFromArray(edallTrueParticle_startPos, startPos, 3);
        output.FillMatrixVarFromArray(edallTrueParticle_endPos, endPos, 3);
        output.FillVectorVar(edallTrueParticle_PDG, truePart->PDG);
        output.FillVectorVar(edallTrueParticle_processEnd, static_cast<Int_t>(truePart->ProcessEnd));
        output.IncrementCounter(ednAllTrueParticles);

        _allTrueParticle_PDG[allTrueParticlesStored] = truePart->PDG;
        _allTrueParticle_processEnd[allTrueParticlesStored] = static_cast<Int_t>(truePart->ProcessEnd);
        for (Int_t c = 0; c < 3; ++c) {
            _allTrueParticle_startPos[allTrueParticlesStored][c] = startPos[c];
            _allTrueParticle_endPos[allTrueParticlesStored][c] = endPos[c];
        }
        ++allTrueParticlesStored;
    }

    auto fillTrueRelations = [&](const AnaTrueParticlePD* trueK0,
                                 Int_t& processEnd,
                                 Int_t& parentPDG,
                                 Float_t parentStart[3],
                                 Float_t parentEnd[3],
                                 Int_t& nDaughters,
                                 Float_t daughterStart[kMaxTrueDaughters*3],
                                 Float_t daughterEnd[kMaxTrueDaughters*3],
                                 Int_t daughterPDG[kMaxTrueDaughters],
                                 Int_t& nSiblings,
                                 Float_t siblingStart[kMaxTrueSiblings*3],
                                 Float_t siblingEnd[kMaxTrueSiblings*3],
                                 Int_t siblingPDG[kMaxTrueSiblings]) {
        processEnd = -1;
        parentPDG = 0;
        for (Int_t i = 0; i < 3; ++i) {
            parentStart[i] = -999.f;
            parentEnd[i] = -999.f;
        }
        nDaughters = 0;
        std::fill_n(daughterStart, kMaxTrueDaughters*3, -999.f);
        std::fill_n(daughterEnd, kMaxTrueDaughters*3, -999.f);
        std::fill_n(daughterPDG, kMaxTrueDaughters, 0);
        nSiblings = 0;
        std::fill_n(siblingStart, kMaxTrueSiblings*3, -999.f);
        std::fill_n(siblingEnd, kMaxTrueSiblings*3, -999.f);
        std::fill_n(siblingPDG, kMaxTrueSiblings, 0);

        if (!trueK0) return;

        processEnd = static_cast<Int_t>(trueK0->ProcessEnd);

        // Parent information
        if (trueK0->ParentID != -1) {
            auto pit = trueParticleByID.find(trueK0->ParentID);
            if (pit != trueParticleByID.end() && pit->second) {
                const AnaTrueParticlePD* parent = pit->second;
                parentPDG = parent->PDG;
                for (Int_t i = 0; i < 3; ++i) {
                    parentStart[i] = parent->Position[i];
                    parentEnd[i] = parent->PositionEnd[i];
                }

                for (size_t idx = 0; idx < parent->Daughters.size(); ++idx) {
                    Int_t siblingID = parent->Daughters[idx];
                    if (siblingID == trueK0->ID) continue;
                    auto sit = trueParticleByID.find(siblingID);
                    if (sit == trueParticleByID.end() || !sit->second) continue;
                    if (nSiblings >= kMaxTrueSiblings) {
                        ++nSiblings;
                        continue;
                    }
                    const AnaTrueParticlePD* sibling = sit->second;
                    const Int_t slot = nSiblings;
                    siblingPDG[slot] = sibling->PDG;
                    siblingStart[slot*3 + 0] = sibling->Position[0];
                    siblingStart[slot*3 + 1] = sibling->Position[1];
                    siblingStart[slot*3 + 2] = sibling->Position[2];
                    siblingEnd[slot*3 + 0] = sibling->PositionEnd[0];
                    siblingEnd[slot*3 + 1] = sibling->PositionEnd[1];
                    siblingEnd[slot*3 + 2] = sibling->PositionEnd[2];
                    ++nSiblings;
                }
            }
        }

        // Daughters information
        for (size_t idx = 0; idx < trueK0->Daughters.size(); ++idx) {
            if (nDaughters >= kMaxTrueDaughters) break; // Cap at kMaxTrueDaughters for storage
            Int_t daughterID = trueK0->Daughters[idx];
            auto dit = trueParticleByID.find(daughterID);
            if (dit == trueParticleByID.end() || !dit->second) continue;
            const AnaTrueParticlePD* daughter = dit->second;
            const Int_t slot = nDaughters;
            daughterPDG[slot] = daughter->PDG;
            daughterStart[slot*3 + 0] = daughter->Position[0];
            daughterStart[slot*3 + 1] = daughter->Position[1];
            daughterStart[slot*3 + 2] = daughter->Position[2];
            daughterEnd[slot*3 + 0] = daughter->PositionEnd[0];
            daughterEnd[slot*3 + 1] = daughter->PositionEnd[1];
            daughterEnd[slot*3 + 2] = daughter->PositionEnd[2];
            ++nDaughters;
        }
    };

    auto resolveAssociatedTrueNeutral = [&](AnaNeutralParticlePD* neutralParticle,
                                            AnaParticlePD* daughter1,
                                            AnaParticlePD* daughter2) -> const AnaTrueParticlePD* {
        if (!neutralParticle || !daughter1 || !daughter2 ||
            !daughter1->TrueObject || !daughter2->TrueObject) {
            return NULL;
        }

        const AnaTrueParticlePD* trueDaughter1 =
            static_cast<const AnaTrueParticlePD*>(daughter1->TrueObject);
        const AnaTrueParticlePD* trueDaughter2 =
            static_cast<const AnaTrueParticlePD*>(daughter2->TrueObject);
        if (!trueDaughter1 || !trueDaughter2) return NULL;

        if (trueDaughter1->ParentID <= 0 || trueDaughter1->ParentID != trueDaughter2->ParentID) {
            return NULL;
        }

        auto trueNeutralIt = trueParticleByID.find(trueDaughter1->ParentID);
        if (trueNeutralIt == trueParticleByID.end() || !trueNeutralIt->second) {
            return NULL;
        }
        const AnaTrueParticlePD* trueNeutral = trueNeutralIt->second;

        const Int_t absPDG = std::abs(trueNeutral->PDG);
        if (!(absPDG == 311 || absPDG == 310 || absPDG == 130)) {
            return NULL;
        }

        if (!neutralParticle->Parent || !neutralParticle->Parent->TrueObject) {
            return NULL;
        }
        const AnaTrueParticlePD* trueRecoParent =
            static_cast<const AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
        if (!trueRecoParent) return NULL;

        bool isTrueDaughterOfRecoParent = false;
        for (size_t idx = 0; idx < trueRecoParent->Daughters.size(); ++idx) {
            if (trueRecoParent->Daughters[idx] == trueNeutral->ID) {
                isTrueDaughterOfRecoParent = true;
                break;
            }
        }
        if (!isTrueDaughterOfRecoParent) return NULL;

        return trueNeutral;
    };

    // Fill K0 candidate data from neutral particle candidates
    if (k0Box && k0Box->neutralParticleCandidates.size() > 0) {
        for (size_t i = 0; i < k0Box->neutralParticleCandidates.size() && i < (size_t)kMaxK0; i++) {
            AnaNeutralParticlePD* neutralParticle = k0Box->neutralParticleCandidates[i];
            if (!neutralParticle) continue;

            // Get daughter particles from annihilation vertex
            Int_t daughter1ID = -1;
            Int_t daughter2ID = -1;
            AnaParticlePD* daughter1 = nullptr;
            AnaParticlePD* daughter2 = nullptr;
            if (neutralParticle->AnnihilationVertex && neutralParticle->AnnihilationVertex->Particles.size() >= 2) {
                daughter1 = neutralParticle->AnnihilationVertex->Particles[0];
                daughter2 = neutralParticle->AnnihilationVertex->Particles[1];
                if (daughter1) {
                    daughter1ID = daughter1->UniqueID;
                }
                if (daughter2) {
                    daughter2ID = daughter2->UniqueID;
                }
            }

            // Get parent ID
            Int_t parentID = neutralParticle->Parent ? neutralParticle->Parent->UniqueID : -1;

            // Get second particle ID from creation vertex
            Int_t secondParticleID = -1;
            if (neutralParticle->CreationVertex && neutralParticle->CreationVertex->SecondParticle) {
                secondParticleID = neutralParticle->CreationVertex->SecondParticle->UniqueID;
            }

            // Fill K0 data
            output.FillVectorVar(edk0_daughter1ID, daughter1ID);
            output.FillVectorVar(edk0_daughter2ID, daughter2ID);
            output.FillVectorVar(edk0_parentID, parentID);
            output.FillVectorVar(edk0_secondParticleID, secondParticleID);

            Float_t parentStartPos[3] = {-999.f, -999.f, -999.f};
            Float_t parentEndPos[3] = {-999.f, -999.f, -999.f};
            Float_t parentEndPosCorrected[3] = {-999.f, -999.f, -999.f};
            Float_t parentTailFitAnchor[3] = {-999.f, -999.f, -999.f};
            Float_t parentTailFitDir[3] = {-999.f, -999.f, -999.f};
            Float_t parentTailFitLength = -999.f;
            Float_t parentTailHitsRaw[kMaxBeamTailHits * 3];
            Float_t parentTailHitsProjected[kMaxBeamTailHits * 3];
            std::fill_n(parentTailHitsRaw, kMaxBeamTailHits * 3, -999.f);
            std::fill_n(parentTailHitsProjected, kMaxBeamTailHits * 3, -999.f);
            Int_t parentTailHitsRawN = 0;
            Int_t parentTailHitsProjectedN = 0;
            Float_t parentTrajDir[3] = {-999.f, -999.f, -999.f};
            Float_t secondTrajDir[3] = {-999.f, -999.f, -999.f};
            Float_t dau1TrajDir[3] = {-999.f, -999.f, -999.f};
            Float_t dau2TrajDir[3] = {-999.f, -999.f, -999.f};
            Int_t parentTrajNpts = 0;
            Int_t secondTrajNpts = 0;
            Int_t dau1TrajNpts = 0;
            Int_t dau2TrajNpts = 0;
            Float_t parentLength = -999.f;
            if (neutralParticle->Parent) {
                parentStartPos[0] = neutralParticle->Parent->PositionStart[0];
                parentStartPos[1] = neutralParticle->Parent->PositionStart[1];
                parentStartPos[2] = neutralParticle->Parent->PositionStart[2];
                parentEndPos[0] = neutralParticle->Parent->PositionEnd[0];
                parentEndPos[1] = neutralParticle->Parent->PositionEnd[1];
                parentEndPos[2] = neutralParticle->Parent->PositionEnd[2];
                parentLength = neutralParticle->Parent->Length;

                if (neutralParticle->Parent->TrajectoryDirectionNPoints > 0) {
                    parentTrajDir[0] = neutralParticle->Parent->TrajectoryDirection.X();
                    parentTrajDir[1] = neutralParticle->Parent->TrajectoryDirection.Y();
                    parentTrajDir[2] = neutralParticle->Parent->TrajectoryDirection.Z();
                    parentTrajNpts = neutralParticle->Parent->TrajectoryDirectionNPoints;
                }

                double fitDistanceFromEndCm = 10.0;
                if (ND::params().HasParameter("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd")) {
                    fitDistanceFromEndCm = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceCreationFromEnd");
                }

                double fitCreationLengthCm = 25.0;
                if (ND::params().HasParameter("neutralKaonAnalysis.TrackFitCreationLength")) {
                    fitCreationLengthCm = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitCreationLength");
                }

                std::vector<double> beamTailFitParams;
                pdAnaUtils::ExtrapolateTrack(neutralParticle->Parent,
                                             beamTailFitParams,
                                             fitCreationLengthCm,
                                             false,
                                             fitDistanceFromEndCm);
                if (beamTailFitParams.size() >= 6 &&
                    beamTailFitParams[0] > -900.0 &&
                    beamTailFitParams[1] > -900.0 &&
                    beamTailFitParams[2] > -900.0) {
                    TVector3 fitDir(beamTailFitParams[3], beamTailFitParams[4], beamTailFitParams[5]);
                    if (fitDir.Mag2() > 1e-10) {
                        fitDir = fitDir.Unit();
                        parentTailFitAnchor[0] = static_cast<Float_t>(beamTailFitParams[0]);
                        parentTailFitAnchor[1] = static_cast<Float_t>(beamTailFitParams[1]);
                        parentTailFitAnchor[2] = static_cast<Float_t>(beamTailFitParams[2]);
                        parentTailFitDir[0] = static_cast<Float_t>(fitDir.X());
                        parentTailFitDir[1] = static_cast<Float_t>(fitDir.Y());
                        parentTailFitDir[2] = static_cast<Float_t>(fitDir.Z());
                        parentTailFitLength = static_cast<Float_t>(fitCreationLengthCm);
                    }
                }

                TVector3 correctedEnd(-999.0, -999.0, -999.0);
                std::vector<TVector3> rawTail;
                std::vector<TVector3> projectedTail;
                const bool projectedOk = pdCreationUtils::ProjectBeamTailOntoStartDirection(
                    neutralParticle->Parent,
                    fitDistanceFromEndCm,
                    correctedEnd,
                    &rawTail,
                    &projectedTail);

                if (projectedOk && correctedEnd.X() > -900 && correctedEnd.Y() > -900 && correctedEnd.Z() > -900) {
                    parentEndPosCorrected[0] = correctedEnd.X();
                    parentEndPosCorrected[1] = correctedEnd.Y();
                    parentEndPosCorrected[2] = correctedEnd.Z();
                }

                parentTailHitsRawN = static_cast<Int_t>(std::min(rawTail.size(), static_cast<size_t>(kMaxBeamTailHits)));
                for (Int_t t = 0; t < parentTailHitsRawN; ++t) {
                    parentTailHitsRaw[t * 3 + 0] = rawTail[t].X();
                    parentTailHitsRaw[t * 3 + 1] = rawTail[t].Y();
                    parentTailHitsRaw[t * 3 + 2] = rawTail[t].Z();
                }

                parentTailHitsProjectedN =
                    static_cast<Int_t>(std::min(projectedTail.size(), static_cast<size_t>(kMaxBeamTailHits)));
                for (Int_t t = 0; t < parentTailHitsProjectedN; ++t) {
                    parentTailHitsProjected[t * 3 + 0] = projectedTail[t].X();
                    parentTailHitsProjected[t * 3 + 1] = projectedTail[t].Y();
                    parentTailHitsProjected[t * 3 + 2] = projectedTail[t].Z();
                }

            }
            output.FillMatrixVarFromArray(edk0_parentStartPos, parentStartPos, 3);
            output.FillMatrixVarFromArray(edk0_parentEndPos, parentEndPos, 3);
            output.FillMatrixVarFromArray(edk0_parentEndPosCorrected, parentEndPosCorrected, 3);
            output.FillMatrixVarFromArray(edk0_parentTailFitAnchor, parentTailFitAnchor, 3);
            output.FillMatrixVarFromArray(edk0_parentTailFitDir, parentTailFitDir, 3);
            output.FillVectorVar(edk0_parentTailFitLength, parentTailFitLength);
            output.FillMatrixVarFromArray(edk0_parentTailHitsRaw, parentTailHitsRaw, kMaxBeamTailHits * 3);
            output.FillMatrixVarFromArray(edk0_parentTailHitsProjected, parentTailHitsProjected, kMaxBeamTailHits * 3);
            output.FillVectorVar(edk0_parentTailHitsRawN, parentTailHitsRawN);
            output.FillVectorVar(edk0_parentTailHitsProjectedN, parentTailHitsProjectedN);
            output.FillVectorVar(edk0_parentLength, parentLength);
            output.FillMatrixVarFromArray(edk0_parentTrajDir, parentTrajDir, 3);
            Float_t parentTrajHist[kParentTrajHistBins * 3];
            BuildParentTrajectoryHistogram(neutralParticle->Parent, parentTrajHist);
            output.FillMatrixVarFromArray(edk0_parentTrajDirHist, parentTrajHist, kParentTrajHistBins * 3);
            output.FillVectorVar(edk0_parentTrajDirNPts, parentTrajNpts);

            AnaParticlePD* secondParticle = (neutralParticle->CreationVertex) ? neutralParticle->CreationVertex->SecondParticle : nullptr;
            if (secondParticle && secondParticle->TrajectoryDirectionNPoints > 0) {
                secondTrajDir[0] = secondParticle->TrajectoryDirection.X();
                secondTrajDir[1] = secondParticle->TrajectoryDirection.Y();
                secondTrajDir[2] = secondParticle->TrajectoryDirection.Z();
                secondTrajNpts = secondParticle->TrajectoryDirectionNPoints;
            }
            if (daughter1 && daughter1->TrajectoryDirectionNPoints > 0) {
                dau1TrajDir[0] = daughter1->TrajectoryDirection.X();
                dau1TrajDir[1] = daughter1->TrajectoryDirection.Y();
                dau1TrajDir[2] = daughter1->TrajectoryDirection.Z();
                dau1TrajNpts = daughter1->TrajectoryDirectionNPoints;
            }
            if (daughter2 && daughter2->TrajectoryDirectionNPoints > 0) {
                dau2TrajDir[0] = daughter2->TrajectoryDirection.X();
                dau2TrajDir[1] = daughter2->TrajectoryDirection.Y();
                dau2TrajDir[2] = daughter2->TrajectoryDirection.Z();
                dau2TrajNpts = daughter2->TrajectoryDirectionNPoints;
            }
            output.FillMatrixVarFromArray(edk0_secondTrajDir, secondTrajDir, 3);
            output.FillVectorVar(edk0_secondTrajDirNPts, secondTrajNpts);
            output.FillMatrixVarFromArray(edk0_dau1TrajDir, dau1TrajDir, 3);
            output.FillVectorVar(edk0_dau1TrajDirNPts, dau1TrajNpts);
            output.FillMatrixVarFromArray(edk0_dau2TrajDir, dau2TrajDir, 3);
            output.FillVectorVar(edk0_dau2TrajDirNPts, dau2TrajNpts);

            // Vertex radii (read from parameters file)
            Float_t creationRadius = 0.0f;
            Float_t annihilationRadius = ND::params().GetParameterD("neutralKaonAnalysis.AnnihilationVertexRadius");

            output.FillVectorVar(edk0_creationVtxRadius, creationRadius);
            output.FillVectorVar(edk0_annihilationVtxRadius, annihilationRadius);

            // Creation vertex position and degeneracy
            Float_t creationPos[3] = {-999, -999, -999};
            Int_t creationDeg = -999;
            if (neutralParticle->CreationVertex) {
                creationPos[0] = (Float_t)neutralParticle->CreationVertex->Position[0];
                creationPos[1] = (Float_t)neutralParticle->CreationVertex->Position[1];
                creationPos[2] = (Float_t)neutralParticle->CreationVertex->Position[2];
                creationDeg = neutralParticle->CreationVertex->Degeneracy;
            } else if (neutralParticle->Parent) {
                creationPos[0] = (Float_t)neutralParticle->Parent->PositionEnd[0];
                creationPos[1] = (Float_t)neutralParticle->Parent->PositionEnd[1];
                creationPos[2] = (Float_t)neutralParticle->Parent->PositionEnd[2];
            }

            output.FillVectorVar(edk0_creationVtxDeg, creationDeg);

            // Creation vertex degeneracy distances
            Float_t creationDegDist[5] = {-999, -999, -999, -999, -999};
            if (neutralParticle->CreationVertex) {
                for (Int_t d = 0; d < 5; d++) {
                    creationDegDist[d] = neutralParticle->CreationVertex->DegDist[d];
                }
            }
            output.FillMatrixVarFromArray(edk0_creationVtxDegDist, creationDegDist, 5);

            // Annihilation vertex position and degeneracy
            Float_t annihilationPos[3] = {-999, -999, -999};
            Float_t annihilationFitPos[3] = {-999, -999, -999};
            Int_t annihilationDeg = -999;
            if (neutralParticle->AnnihilationVertex) {
                annihilationPos[0] = (Float_t)neutralParticle->AnnihilationVertex->PositionPandora[0];
                annihilationPos[1] = (Float_t)neutralParticle->AnnihilationVertex->PositionPandora[1];
                annihilationPos[2] = (Float_t)neutralParticle->AnnihilationVertex->PositionPandora[2];
                annihilationFitPos[0] = (Float_t)neutralParticle->AnnihilationVertex->PositionFit[0];
                annihilationFitPos[1] = (Float_t)neutralParticle->AnnihilationVertex->PositionFit[1];
                annihilationFitPos[2] = (Float_t)neutralParticle->AnnihilationVertex->PositionFit[2];
                annihilationDeg = neutralParticle->AnnihilationVertex->Degeneracy;
            }
            output.FillVectorVar(edk0_annihilationVtxDeg, annihilationDeg);

            // Annihilation vertex degeneracy distances
            Float_t annihilationDegDist[5] = {-999, -999, -999, -999, -999};
            output.FillMatrixVarFromArray(edk0_annihilationVtxDegDist, annihilationDegDist, 5);

            // K0 trajectory (from creation to annihilation)
            Float_t startPos[3] = {creationPos[0], creationPos[1], creationPos[2]};
            Float_t endPos[3] = {annihilationPos[0], annihilationPos[1], annihilationPos[2]};

            output.FillMatrixVarFromArray(edk0_creationVtxPos, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annihilationVtxPos, annihilationPos, 3);
            output.FillMatrixVarFromArray(edk0_annihilationVtxFitPos, annihilationFitPos, 3);
            output.FillMatrixVarFromArray(edk0_startPos, startPos, 3);
            output.FillMatrixVarFromArray(edk0_endPos, endPos, 3);

            // Extract fitted line parameters from annihilation vertex
            Float_t fitLine1Start[3] = {-999, -999, -999};
            Float_t fitLine1Dir[3] = {-999, -999, -999};
            Float_t fitLine2Start[3] = {-999, -999, -999};
            Float_t fitLine2Dir[3] = {-999, -999, -999};
            Float_t pandoraLine1Dir[3] = {-999, -999, -999};
            Float_t pandoraLine2Dir[3] = {-999, -999, -999};
            Float_t closestPt1[3] = {-999, -999, -999};
            Float_t closestPt2[3] = {-999, -999, -999};
            Float_t pandoraClosestPt1[3] = {-999, -999, -999};
            Float_t pandoraClosestPt2[3] = {-999, -999, -999};
            Float_t vtxMomentumDirFit[3] = {-999, -999, -999};
            Float_t trueK0Dir[3] = {-999, -999, -999};
            Float_t fitLineLength = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitLength");
            Float_t fitLineDistanceFromStart = ND::params().GetParameterD("neutralKaonAnalysis.TrackFitDistanceFromStart");

            // For annihilation vertex, draw both Pandora and fit-based line definitions.
            if (neutralParticle->AnnihilationVertex &&
                neutralParticle->AnnihilationVertex->Particles.size() >= 2) {
                AnaParticlePD* d1 = neutralParticle->AnnihilationVertex->Particles[0];
                AnaParticlePD* d2 = neutralParticle->AnnihilationVertex->Particles[1];

                if (d1 && d2) {
                    // Start points are the physical anchors nearest TrackFitDistanceFromStart.
                    fitLine1Start[0] = d1->PositionStart[0];
                    fitLine1Start[1] = d1->PositionStart[1];
                    fitLine1Start[2] = d1->PositionStart[2];
                    fitLine2Start[0] = d2->PositionStart[0];
                    fitLine2Start[1] = d2->PositionStart[1];
                    fitLine2Start[2] = d2->PositionStart[2];

                    // Pandora directions
                    pandoraLine1Dir[0] = d1->DirectionStart[0];
                    pandoraLine1Dir[1] = d1->DirectionStart[1];
                    pandoraLine1Dir[2] = d1->DirectionStart[2];
                    pandoraLine2Dir[0] = d2->DirectionStart[0];
                    pandoraLine2Dir[1] = d2->DirectionStart[1];
                    pandoraLine2Dir[2] = d2->DirectionStart[2];

                    // Extrapolated-fit directions (first trackFitLength cm from start)
                    std::vector<double> line1Params;
                    std::vector<double> line2Params;
                    pdAnaUtils::ExtrapolateTrack(d1, line1Params, fitLineLength, true, fitLineDistanceFromStart);
                    pdAnaUtils::ExtrapolateTrack(d2, line2Params, fitLineLength, true, fitLineDistanceFromStart);
                    const bool line1Valid = (line1Params.size() >= 6 &&
                                            line1Params[0] > -900.0 &&
                                            line1Params[1] > -900.0 &&
                                            line1Params[2] > -900.0);
                    const bool line2Valid = (line2Params.size() >= 6 &&
                                            line2Params[0] > -900.0 &&
                                            line2Params[1] > -900.0 &&
                                            line2Params[2] > -900.0);
                    if (line1Valid) {
                        fitLine1Start[0] = line1Params[0];
                        fitLine1Start[1] = line1Params[1];
                        fitLine1Start[2] = line1Params[2];
                        fitLine1Dir[0] = line1Params[3];
                        fitLine1Dir[1] = line1Params[4];
                        fitLine1Dir[2] = line1Params[5];
                    } else {
                        fitLine1Dir[0] = pandoraLine1Dir[0];
                        fitLine1Dir[1] = pandoraLine1Dir[1];
                        fitLine1Dir[2] = pandoraLine1Dir[2];
                    }
                    if (line2Valid) {
                        fitLine2Start[0] = line2Params[0];
                        fitLine2Start[1] = line2Params[1];
                        fitLine2Start[2] = line2Params[2];
                        fitLine2Dir[0] = line2Params[3];
                        fitLine2Dir[1] = line2Params[4];
                        fitLine2Dir[2] = line2Params[5];
                    } else {
                        fitLine2Dir[0] = pandoraLine2Dir[0];
                        fitLine2Dir[1] = pandoraLine2Dir[1];
                        fitLine2Dir[2] = pandoraLine2Dir[2];
                    }

                    // Force fitted directions to follow particle travel direction.
                    TVector3 p1(pandoraLine1Dir[0], pandoraLine1Dir[1], pandoraLine1Dir[2]);
                    TVector3 f1(fitLine1Dir[0], fitLine1Dir[1], fitLine1Dir[2]);
                    if (p1.Mag2() > 1e-10 && f1.Mag2() > 1e-10 && f1.Dot(p1) < 0) {
                        fitLine1Dir[0] *= -1.f;
                        fitLine1Dir[1] *= -1.f;
                        fitLine1Dir[2] *= -1.f;
                    }

                    TVector3 p2(pandoraLine2Dir[0], pandoraLine2Dir[1], pandoraLine2Dir[2]);
                    TVector3 f2(fitLine2Dir[0], fitLine2Dir[1], fitLine2Dir[2]);
                    if (p2.Mag2() > 1e-10 && f2.Mag2() > 1e-10 && f2.Dot(p2) < 0) {
                        fitLine2Dir[0] *= -1.f;
                        fitLine2Dir[1] *= -1.f;
                        fitLine2Dir[2] *= -1.f;
                    }

                    // Reconstructed vertex momentum direction from daughter momenta and fit directions.
                    if (d1->Momentum > 0.f && d2->Momentum > 0.f) {
                        TVector3 dir1(fitLine1Dir[0], fitLine1Dir[1], fitLine1Dir[2]);
                        TVector3 dir2(fitLine2Dir[0], fitLine2Dir[1], fitLine2Dir[2]);
                        if (dir1.Mag2() > 1e-10 && dir2.Mag2() > 1e-10) {
                            const TVector3 pVec = d1->Momentum * dir1.Unit() + d2->Momentum * dir2.Unit();
                            if (pVec.Mag2() > 1e-10) {
                                const TVector3 u = pVec.Unit();
                                vtxMomentumDirFit[0] = u.X();
                                vtxMomentumDirFit[1] = u.Y();
                                vtxMomentumDirFit[2] = u.Z();
                            }
                        }
                    }
                }

                closestPt1[0] = neutralParticle->AnnihilationVertex->ClosestPointFit1[0];
                closestPt1[1] = neutralParticle->AnnihilationVertex->ClosestPointFit1[1];
                closestPt1[2] = neutralParticle->AnnihilationVertex->ClosestPointFit1[2];
                closestPt2[0] = neutralParticle->AnnihilationVertex->ClosestPointFit2[0];
                closestPt2[1] = neutralParticle->AnnihilationVertex->ClosestPointFit2[1];
                closestPt2[2] = neutralParticle->AnnihilationVertex->ClosestPointFit2[2];
                pandoraClosestPt1[0] = neutralParticle->AnnihilationVertex->ClosestPointPandora1[0];
                pandoraClosestPt1[1] = neutralParticle->AnnihilationVertex->ClosestPointPandora1[1];
                pandoraClosestPt1[2] = neutralParticle->AnnihilationVertex->ClosestPointPandora1[2];
                pandoraClosestPt2[0] = neutralParticle->AnnihilationVertex->ClosestPointPandora2[0];
                pandoraClosestPt2[1] = neutralParticle->AnnihilationVertex->ClosestPointPandora2[1];
                pandoraClosestPt2[2] = neutralParticle->AnnihilationVertex->ClosestPointPandora2[2];
            }

            output.FillMatrixVarFromArray(edk0_annVtx_fitLine1Start, fitLine1Start, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine1Dir, fitLine1Dir, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine2Start, fitLine2Start, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine2Dir, fitLine2Dir, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraLine1Dir, pandoraLine1Dir, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraLine2Dir, pandoraLine2Dir, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_closestPt1, closestPt1, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_closestPt2, closestPt2, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraClosestPt1, pandoraClosestPt1, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraClosestPt2, pandoraClosestPt2, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_momentumDirFit, vtxMomentumDirFit, 3);

            // Extract fitted line parameters from creation vertex
            Float_t creationFitLineBeamStart[3] = {-999, -999, -999};
            Float_t creationFitLineBeamDir[3] = {-999, -999, -999};
            Float_t creationFitLineSecondStart[3] = {-999, -999, -999};
            Float_t creationFitLineSecondDir[3] = {-999, -999, -999};
            Float_t creationClosestPtBeam[3] = {-999, -999, -999};
            Float_t creationClosestPtSecond[3] = {-999, -999, -999};

            if (neutralParticle->CreationVertex &&
                neutralParticle->CreationVertex->FittedLineParams.size() >= 2) {
                // Get fitted line params (similar to annihilation vertex)
                std::vector<double> beamLine = neutralParticle->CreationVertex->FittedLineParams[0];
                std::vector<double> secondLine = neutralParticle->CreationVertex->FittedLineParams[1];

                if (beamLine.size() >= 6 && secondLine.size() >= 6) {
                    creationFitLineBeamStart[0] = beamLine[0];
                    creationFitLineBeamStart[1] = beamLine[1];
                    creationFitLineBeamStart[2] = beamLine[2];
                    creationFitLineBeamDir[0] = beamLine[3];
                    creationFitLineBeamDir[1] = beamLine[4];
                    creationFitLineBeamDir[2] = beamLine[5];

                    creationFitLineSecondStart[0] = secondLine[0];
                    creationFitLineSecondStart[1] = secondLine[1];
                    creationFitLineSecondStart[2] = secondLine[2];
                    creationFitLineSecondDir[0] = secondLine[3];
                    creationFitLineSecondDir[1] = secondLine[4];
                    creationFitLineSecondDir[2] = secondLine[5];
                }
            }

            if (neutralParticle->CreationVertex) {
                creationClosestPtBeam[0] = neutralParticle->CreationVertex->ClosestPointBeam[0];
                creationClosestPtBeam[1] = neutralParticle->CreationVertex->ClosestPointBeam[1];
                creationClosestPtBeam[2] = neutralParticle->CreationVertex->ClosestPointBeam[2];
                creationClosestPtSecond[0] = neutralParticle->CreationVertex->ClosestPointSecond[0];
                creationClosestPtSecond[1] = neutralParticle->CreationVertex->ClosestPointSecond[1];
                creationClosestPtSecond[2] = neutralParticle->CreationVertex->ClosestPointSecond[2];
            }

            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineBeamStart, creationFitLineBeamStart, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineBeamDir, creationFitLineBeamDir, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineSecondStart, creationFitLineSecondStart, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineSecondDir, creationFitLineSecondDir, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_closestPtBeam, creationClosestPtBeam, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_closestPtSecond, creationClosestPtSecond, 3);

            output.FillVectorVar(edk0_fitLineLength, fitLineLength);

            // Fill true K0 trajectory from the neutral-particle associated TrueObject.
            // This TrueObject is assigned during neutral creation.
            Int_t hasTrueObject = 0;
            Float_t trueStartPos[3] = {-999, -999, -999};
            Float_t trueEndPos[3] = {-999, -999, -999};
            Int_t truePDG = -999;
            const AnaTrueParticlePD* associatedTrueK0 = NULL;

            associatedTrueK0 = resolveAssociatedTrueNeutral(neutralParticle, daughter1, daughter2);
            if (associatedTrueK0) {
                hasTrueObject = 1;
                trueStartPos[0] = associatedTrueK0->Position[0];
                trueStartPos[1] = associatedTrueK0->Position[1];
                trueStartPos[2] = associatedTrueK0->Position[2];
                trueEndPos[0] = associatedTrueK0->PositionEnd[0];
                trueEndPos[1] = associatedTrueK0->PositionEnd[1];
                trueEndPos[2] = associatedTrueK0->PositionEnd[2];
                truePDG = associatedTrueK0->PDG;
                TVector3 trueDirVec(associatedTrueK0->PositionEnd[0] - associatedTrueK0->Position[0],
                                    associatedTrueK0->PositionEnd[1] - associatedTrueK0->Position[1],
                                    associatedTrueK0->PositionEnd[2] - associatedTrueK0->Position[2]);
                if (trueDirVec.Mag2() > 1e-10) {
                    trueDirVec = trueDirVec.Unit();
                    trueK0Dir[0] = trueDirVec.X();
                    trueK0Dir[1] = trueDirVec.Y();
                    trueK0Dir[2] = trueDirVec.Z();
                }
            }

            Float_t trueDecayVtxFromRecoDaughters[3] = {-999.f, -999.f, -999.f};
            if (daughter1 && daughter2 && daughter1->TrueObject && daughter2->TrueObject) {
                const AnaTrueParticlePD* tr1 = static_cast<const AnaTrueParticlePD*>(daughter1->TrueObject);
                const AnaTrueParticlePD* tr2 = static_cast<const AnaTrueParticlePD*>(daughter2->TrueObject);
                if (tr1 && tr2 && tr1->Position[0] > -900.f && tr1->Position[1] > -900.f &&
                    tr1->Position[2] > -900.f && tr2->Position[0] > -900.f && tr2->Position[1] > -900.f &&
                    tr2->Position[2] > -900.f) {
                    trueDecayVtxFromRecoDaughters[0] =
                        0.5f * static_cast<Float_t>(tr1->Position[0] + tr2->Position[0]);
                    trueDecayVtxFromRecoDaughters[1] =
                        0.5f * static_cast<Float_t>(tr1->Position[1] + tr2->Position[1]);
                    trueDecayVtxFromRecoDaughters[2] =
                        0.5f * static_cast<Float_t>(tr1->Position[2] + tr2->Position[2]);
                }
            }

            Int_t processEndCode;
            Int_t parentPDG;
            Float_t parentStart[3];
            Float_t parentEnd[3];
            Int_t nDaughters;
            Float_t daughterStartFlat[kMaxTrueDaughters*3];
            Float_t daughterEndFlat[kMaxTrueDaughters*3];
            Int_t daughterPDGArr[kMaxTrueDaughters];
            Int_t nParentDaughters;
            Float_t parentDaughterStartFlat[kMaxTrueDaughters*3];
            Float_t parentDaughterEndFlat[kMaxTrueDaughters*3];
            Int_t parentDaughterPDGArr[kMaxTrueDaughters];
            Int_t nSiblings;
            Float_t siblingStartFlat[kMaxTrueSiblings*3];
            Float_t siblingEndFlat[kMaxTrueSiblings*3];
            Int_t siblingPDGArr[kMaxTrueSiblings];
            fillTrueRelations(associatedTrueK0, processEndCode, parentPDG, parentStart, parentEnd, nDaughters, daughterStartFlat, daughterEndFlat, daughterPDGArr, nSiblings, siblingStartFlat, siblingEndFlat, siblingPDGArr);

            auto fillParentDaughters = [&](const AnaTrueParticlePD* trueParent,
                                           Int_t& outN,
                                           Float_t startFlat[kMaxTrueDaughters*3],
                                           Float_t endFlat[kMaxTrueDaughters*3],
                                           Int_t pdgArr[kMaxTrueDaughters]) {
                outN = 0;
                std::fill_n(startFlat, kMaxTrueDaughters * 3, -999.f);
                std::fill_n(endFlat, kMaxTrueDaughters * 3, -999.f);
                std::fill_n(pdgArr, kMaxTrueDaughters, 0);
                if (!trueParent) return;

                std::set<Int_t> daughterIDs;
                for (size_t d = 0; d < trueParent->Daughters.size(); ++d) {
                    daughterIDs.insert(trueParent->Daughters[d]);
                }

                for (const auto& entry : trueParticleByID) {
                    const AnaTrueParticlePD* candidate = entry.second;
                    if (!candidate) continue;
                    if (candidate->ParentID == trueParent->ID) {
                        daughterIDs.insert(candidate->ID);
                    }
                }

                for (std::set<Int_t>::const_iterator it = daughterIDs.begin();
                     it != daughterIDs.end() && outN < kMaxTrueDaughters; ++it) {
                    const Int_t daughterID = *it;
                    auto dit = trueParticleByID.find(daughterID);
                    if (dit == trueParticleByID.end() || !dit->second) continue;
                    const AnaTrueParticlePD* trueDaughter = dit->second;
                    const Int_t slot = outN;
                    pdgArr[slot] = trueDaughter->PDG;
                    startFlat[slot * 3 + 0] = trueDaughter->Position[0];
                    startFlat[slot * 3 + 1] = trueDaughter->Position[1];
                    startFlat[slot * 3 + 2] = trueDaughter->Position[2];
                    endFlat[slot * 3 + 0] = trueDaughter->PositionEnd[0];
                    endFlat[slot * 3 + 1] = trueDaughter->PositionEnd[1];
                    endFlat[slot * 3 + 2] = trueDaughter->PositionEnd[2];
                    ++outN;
                }
            };

            const AnaTrueParticlePD* trueParentForDisplay = nullptr;
            if (associatedTrueK0 && associatedTrueK0->ParentID > 0) {
                auto pit = trueParticleByID.find(associatedTrueK0->ParentID);
                if (pit != trueParticleByID.end()) trueParentForDisplay = pit->second;
            }
            if (!trueParentForDisplay && neutralParticle->Parent && neutralParticle->Parent->TrueObject) {
                trueParentForDisplay = static_cast<const AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
            }
            fillParentDaughters(trueParentForDisplay, nParentDaughters, parentDaughterStartFlat, parentDaughterEndFlat, parentDaughterPDGArr);

            if (!associatedTrueK0) {
                if (neutralParticle->Parent && neutralParticle->Parent->TrueObject) {
                    const AnaTrueParticlePD* trueRecoParent =
                        static_cast<const AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
                    if (trueRecoParent) {
                        parentPDG = trueRecoParent->PDG;
                        processEndCode = static_cast<Int_t>(trueRecoParent->ProcessEnd);
                        for (Int_t c = 0; c < 3; ++c) {
                            parentStart[c] = trueRecoParent->Position[c];
                            parentEnd[c] = trueRecoParent->PositionEnd[c];
                        }
                        if (parentEnd[0] > -900.f && parentEnd[1] > -900.f && parentEnd[2] > -900.f) {
                            trueStartPos[0] = parentEnd[0];
                            trueStartPos[1] = parentEnd[1];
                            trueStartPos[2] = parentEnd[2];
                        }
                    }
                }

                nDaughters = 0;
                std::fill_n(daughterStartFlat, kMaxTrueDaughters * 3, -999.f);
                std::fill_n(daughterEndFlat, kMaxTrueDaughters * 3, -999.f);
                std::fill_n(daughterPDGArr, kMaxTrueDaughters, 0);

                if (neutralParticle->Parent && neutralParticle->Parent->TrueObject) {
                    const AnaTrueParticlePD* trueRecoParent =
                        static_cast<const AnaTrueParticlePD*>(neutralParticle->Parent->TrueObject);
                    if (trueRecoParent) {
                        const Int_t maxDaughters = static_cast<Int_t>(std::min(trueRecoParent->Daughters.size(), static_cast<size_t>(kMaxTrueDaughters)));
                        for (Int_t d = 0; d < maxDaughters; ++d) {
                            const Int_t daughterID = trueRecoParent->Daughters[d];
                            auto dit = trueParticleByID.find(daughterID);
                            if (dit == trueParticleByID.end() || !dit->second) continue;
                            const AnaTrueParticlePD* trueDaughter = dit->second;
                            const Int_t slot = nDaughters;
                            daughterPDGArr[slot] = trueDaughter->PDG;
                            daughterStartFlat[slot * 3 + 0] = trueDaughter->Position[0];
                            daughterStartFlat[slot * 3 + 1] = trueDaughter->Position[1];
                            daughterStartFlat[slot * 3 + 2] = trueDaughter->Position[2];
                            daughterEndFlat[slot * 3 + 0] = trueDaughter->PositionEnd[0];
                            daughterEndFlat[slot * 3 + 1] = trueDaughter->PositionEnd[1];
                            daughterEndFlat[slot * 3 + 2] = trueDaughter->PositionEnd[2];
                            ++nDaughters;
                        }
                    }
                }

                if (trueDecayVtxFromRecoDaughters[0] > -900.f && trueDecayVtxFromRecoDaughters[1] > -900.f &&
                    trueDecayVtxFromRecoDaughters[2] > -900.f) {
                    trueEndPos[0] = trueDecayVtxFromRecoDaughters[0];
                    trueEndPos[1] = trueDecayVtxFromRecoDaughters[1];
                    trueEndPos[2] = trueDecayVtxFromRecoDaughters[2];
                }

                TVector3 trueFallbackDir(trueEndPos[0] - trueStartPos[0],
                                         trueEndPos[1] - trueStartPos[1],
                                         trueEndPos[2] - trueStartPos[2]);
                if (trueFallbackDir.Mag2() > 1e-10 &&
                    trueStartPos[0] > -900.f && trueStartPos[1] > -900.f && trueStartPos[2] > -900.f &&
                    trueEndPos[0] > -900.f && trueEndPos[1] > -900.f && trueEndPos[2] > -900.f) {
                    trueFallbackDir = trueFallbackDir.Unit();
                    trueK0Dir[0] = trueFallbackDir.X();
                    trueK0Dir[1] = trueFallbackDir.Y();
                    trueK0Dir[2] = trueFallbackDir.Z();
                }
            }

            output.FillVectorVar(edk0_hasTrueObject, hasTrueObject);
            output.FillMatrixVarFromArray(edk0_trueStartPos, trueStartPos, 3);
            output.FillMatrixVarFromArray(edk0_trueEndPos, trueEndPos, 3);
            output.FillMatrixVarFromArray(edk0_trueK0Dir, trueK0Dir, 3);
            output.FillMatrixVarFromArray(edk0_trueDecayVtxFromRecoDaughters, trueDecayVtxFromRecoDaughters, 3);
            output.FillVectorVar(edk0_truePDG, truePDG);
            output.FillVectorVar(edk0_trueProcessEnd, processEndCode);
            output.FillVectorVar(edk0_trueParentPDG, parentPDG);
            output.FillMatrixVarFromArray(edk0_trueParentStartPos, parentStart, 3);
            output.FillMatrixVarFromArray(edk0_trueParentEndPos, parentEnd, 3);
            output.FillVectorVar(edk0_trueParentNDaughters, nParentDaughters);
            output.FillMatrixVarFromArray(edk0_trueParentDaughterStartPos, parentDaughterStartFlat, kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueParentDaughterEndPos, parentDaughterEndFlat, kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueParentDaughterPDG, parentDaughterPDGArr, kMaxTrueDaughters);
            output.FillVectorVar(edk0_trueNDaughters, nDaughters);
            output.FillMatrixVarFromArray(edk0_trueDaughterStartPos, daughterStartFlat, kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueDaughterEndPos, daughterEndFlat, kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueDaughterPDG, daughterPDGArr, kMaxTrueDaughters);
            output.FillVectorVar(edk0_trueNSiblings, nSiblings);
            output.FillMatrixVarFromArray(edk0_trueSiblingStartPos, siblingStartFlat, kMaxTrueSiblings*3);
            output.FillMatrixVarFromArray(edk0_trueSiblingEndPos, siblingEndFlat, kMaxTrueSiblings*3);
            output.FillMatrixVarFromArray(edk0_trueSiblingPDG, siblingPDGArr, kMaxTrueSiblings);

            if (associatedTrueK0) {
                matchedTrueKaons.insert(associatedTrueK0);
            }
            if (i < kMaxK0) {
                _k0_truePDG[i] = truePDG;
                _k0_hasTrueObject[i] = hasTrueObject;
                _k0_trueProcessEnd[i] = processEndCode;
                _k0_trueParentPDG[i] = parentPDG;
                for (Int_t c = 0; c < 3; ++c) {
                    _k0_trueStartPos[i][c] = trueStartPos[c];
                    _k0_trueEndPos[i][c] = trueEndPos[c];
                    _k0_trueDecayVtxFromRecoDaughters[i][c] = trueDecayVtxFromRecoDaughters[c];
                    _k0_trueParentStartPos[i][c] = parentStart[c];
                    _k0_trueParentEndPos[i][c] = parentEnd[c];
                }
                _k0_trueParentNDaughters[i] = nParentDaughters;
                for (Int_t c = 0; c < kMaxTrueDaughters*3; ++c) {
                    _k0_trueParentDaughterStartPos[i][c] = parentDaughterStartFlat[c];
                    _k0_trueParentDaughterEndPos[i][c] = parentDaughterEndFlat[c];
                }
                for (Int_t c = 0; c < kMaxTrueDaughters; ++c) {
                    _k0_trueParentDaughterPDG[i][c] = parentDaughterPDGArr[c];
                }
                _k0_trueNDaughters[i] = nDaughters;
                for (Int_t c = 0; c < kMaxTrueDaughters*3; ++c) {
                    _k0_trueDaughterStartPos[i][c] = daughterStartFlat[c];
                    _k0_trueDaughterEndPos[i][c] = daughterEndFlat[c];
                }
                for (Int_t c = 0; c < kMaxTrueDaughters; ++c) {
                    _k0_trueDaughterPDG[i][c] = daughterPDGArr[c];
                }
                for (Int_t c = 0; c < kMaxTrueSiblings*3; ++c) {
                    _k0_trueSiblingStartPos[i][c] = siblingStartFlat[c];
                    _k0_trueSiblingEndPos[i][c] = siblingEndFlat[c];
                }
                for (Int_t c = 0; c < kMaxTrueSiblings; ++c) {
                    _k0_trueSiblingPDG[i][c] = siblingPDGArr[c];
                }
            }

            // Increment K0 counter
            output.IncrementCounter(ednK0Candidates);
        }
    }

    // Store standalone true K0 trajectories (not associated to neutral particles)
    Int_t standaloneStored = 0;
    for (int i = 0; i < event.nTrueParticles && standaloneStored < kMaxTrueK0; ++i) {
        AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(event.TrueParticles[i]);
        if (!truePart) continue;

        if (std::abs(truePart->PDG) != 310) continue;

        if (matchedTrueKaons.find(truePart) != matchedTrueKaons.end()) continue;

        Float_t trueStartPos[3] = {
            (Float_t)truePart->Position[0],
            (Float_t)truePart->Position[1],
            (Float_t)truePart->Position[2]
        };
        Float_t trueEndPos[3] = {
            (Float_t)truePart->PositionEnd[0],
            (Float_t)truePart->PositionEnd[1],
            (Float_t)truePart->PositionEnd[2]
        };

        output.FillMatrixVarFromArray(edtrueK0_startPos, trueStartPos, 3);
        output.FillMatrixVarFromArray(edtrueK0_endPos, trueEndPos, 3);
        output.FillVectorVar(edtrueK0_PDG, truePart->PDG);
        Int_t processEndCode;
        Int_t parentPDG;
        Float_t parentStart[3];
        Float_t parentEnd[3];
        Int_t nDaughters;
        Float_t daughterStartFlat[kMaxTrueDaughters*3];
        Float_t daughterEndFlat[kMaxTrueDaughters*3];
        Int_t daughterPDGArr[kMaxTrueDaughters];
        Int_t nSiblings;
        Float_t siblingStartFlat[kMaxTrueSiblings*3];
        Float_t siblingEndFlat[kMaxTrueSiblings*3];
        Int_t siblingPDGArr[kMaxTrueSiblings];
        fillTrueRelations(truePart, processEndCode, parentPDG, parentStart, parentEnd, nDaughters, daughterStartFlat, daughterEndFlat, daughterPDGArr, nSiblings, siblingStartFlat, siblingEndFlat, siblingPDGArr);

        output.FillVectorVar(edtrueK0_processEnd, processEndCode);
        output.FillVectorVar(edtrueK0_parentPDG, parentPDG);
        output.FillMatrixVarFromArray(edtrueK0_parentStartPos, parentStart, 3);
        output.FillMatrixVarFromArray(edtrueK0_parentEndPos, parentEnd, 3);
        output.FillVectorVar(edtrueK0_nDaughters, nDaughters);
        output.FillMatrixVarFromArray(edtrueK0_daughterStartPos, daughterStartFlat, kMaxTrueDaughters*3);
        output.FillMatrixVarFromArray(edtrueK0_daughterEndPos, daughterEndFlat, kMaxTrueDaughters*3);
        output.FillMatrixVarFromArray(edtrueK0_daughterPDG, daughterPDGArr, kMaxTrueDaughters);
        output.FillVectorVar(edtrueK0_nSiblings, nSiblings);
        output.FillMatrixVarFromArray(edtrueK0_siblingStartPos, siblingStartFlat, kMaxTrueSiblings*3);
        output.FillMatrixVarFromArray(edtrueK0_siblingEndPos, siblingEndFlat, kMaxTrueSiblings*3);
        output.FillMatrixVarFromArray(edtrueK0_siblingPDG, siblingPDGArr, kMaxTrueSiblings);

        if (standaloneStored < kMaxTrueK0) {
            _trueK0_PDG[standaloneStored] = truePart->PDG;
            _trueK0_processEnd[standaloneStored] = processEndCode;
            _trueK0_parentPDG[standaloneStored] = parentPDG;
            for (Int_t c = 0; c < 3; ++c) {
                _trueK0_startPos[standaloneStored][c] = trueStartPos[c];
                _trueK0_endPos[standaloneStored][c] = trueEndPos[c];
                _trueK0_parentStartPos[standaloneStored][c] = parentStart[c];
                _trueK0_parentEndPos[standaloneStored][c] = parentEnd[c];
            }
            _trueK0_nDaughters[standaloneStored] = nDaughters;
            for (Int_t c = 0; c < kMaxTrueDaughters*3; ++c) {
                _trueK0_daughterStartPos[standaloneStored][c] = daughterStartFlat[c];
                _trueK0_daughterEndPos[standaloneStored][c] = daughterEndFlat[c];
            }
            for (Int_t c = 0; c < kMaxTrueDaughters; ++c) {
                _trueK0_daughterPDG[standaloneStored][c] = daughterPDGArr[c];
            }
            for (Int_t c = 0; c < kMaxTrueSiblings*3; ++c) {
                _trueK0_siblingStartPos[standaloneStored][c] = siblingStartFlat[c];
                _trueK0_siblingEndPos[standaloneStored][c] = siblingEndFlat[c];
            }
            for (Int_t c = 0; c < kMaxTrueSiblings; ++c) {
                _trueK0_siblingPDG[standaloneStored][c] = siblingPDGArr[c];
            }
        }
        output.IncrementCounter(ednTrueK0);

        ++standaloneStored;
    }

    bool hasRecoCandidates = (k0Box && !k0Box->neutralParticleCandidates.empty());
    if (!hasRecoCandidates && standaloneStored > 0) {
        Int_t pseudoCount = std::min(standaloneStored, kMaxK0);
        Float_t zero5[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
        Float_t zero3[3] = {0.f, 0.f, 0.f};
        Float_t invalid3[3] = {-999.f, -999.f, -999.f};
        Float_t invalidTail[kMaxBeamTailHits * 3];
        std::fill_n(invalidTail, kMaxBeamTailHits * 3, -999.f);
        for (Int_t idx = 0; idx < pseudoCount; ++idx) {
            Float_t creationPos[3] = {
                _trueK0_startPos[idx][0],
                _trueK0_startPos[idx][1],
                _trueK0_startPos[idx][2]
            };
            Float_t annihilationPos[3] = {
                _trueK0_endPos[idx][0],
                _trueK0_endPos[idx][1],
                _trueK0_endPos[idx][2]
            };

            output.FillVectorVar(edk0_daughter1ID, -1);
            output.FillVectorVar(edk0_daughter2ID, -1);
            output.FillVectorVar(edk0_parentID, -1);
            output.FillVectorVar(edk0_secondParticleID, -1);
            output.FillVectorVar(edk0_creationVtxRadius, 0.f);
            output.FillVectorVar(edk0_annihilationVtxRadius, 0.f);
            output.FillVectorVar(edk0_creationVtxDeg, 0);
            output.FillVectorVar(edk0_annihilationVtxDeg, 0);
            output.FillMatrixVarFromArray(edk0_creationVtxDegDist, zero5, 5);
            output.FillMatrixVarFromArray(edk0_annihilationVtxDegDist, zero5, 5);
            output.FillMatrixVarFromArray(edk0_creationVtxPos, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annihilationVtxPos, annihilationPos, 3);
            output.FillMatrixVarFromArray(edk0_annihilationVtxFitPos, annihilationPos, 3);
            output.FillMatrixVarFromArray(edk0_startPos, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_endPos, annihilationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine1Start, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine1Dir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine2Start, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_fitLine2Dir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraLine1Dir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraLine2Dir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_closestPt1, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_closestPt2, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraClosestPt1, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_pandoraClosestPt2, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_annVtx_momentumDirFit, zero3, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineBeamStart, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineBeamDir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineSecondStart, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_fitLineSecondDir, zero3, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_closestPtBeam, creationPos, 3);
            output.FillMatrixVarFromArray(edk0_creationVtx_closestPtSecond, creationPos, 3);
            output.FillVectorVar(edk0_fitLineLength, 0.f);
            output.FillMatrixVarFromArray(edk0_parentEndPosCorrected, invalid3, 3);
            output.FillMatrixVarFromArray(edk0_parentTailFitAnchor, invalid3, 3);
            output.FillMatrixVarFromArray(edk0_parentTailFitDir, invalid3, 3);
            output.FillVectorVar(edk0_parentTailFitLength, -999.f);
            output.FillMatrixVarFromArray(edk0_parentTailHitsRaw, invalidTail, kMaxBeamTailHits * 3);
            output.FillMatrixVarFromArray(edk0_parentTailHitsProjected, invalidTail, kMaxBeamTailHits * 3);
            output.FillVectorVar(edk0_parentTailHitsRawN, 0);
            output.FillVectorVar(edk0_parentTailHitsProjectedN, 0);
            output.FillVectorVar(edk0_hasTrueObject, 1);
            output.FillMatrixVarFromArray(edk0_trueStartPos, _trueK0_startPos[idx], 3);
            output.FillMatrixVarFromArray(edk0_trueEndPos, _trueK0_endPos[idx], 3);
            TVector3 truePseudoDir(_trueK0_endPos[idx][0] - _trueK0_startPos[idx][0],
                                   _trueK0_endPos[idx][1] - _trueK0_startPos[idx][1],
                                   _trueK0_endPos[idx][2] - _trueK0_startPos[idx][2]);
            Float_t truePseudoDirArr[3] = {-999.f, -999.f, -999.f};
            if (truePseudoDir.Mag2() > 1e-10) {
                truePseudoDir = truePseudoDir.Unit();
                truePseudoDirArr[0] = truePseudoDir.X();
                truePseudoDirArr[1] = truePseudoDir.Y();
                truePseudoDirArr[2] = truePseudoDir.Z();
            }
            output.FillMatrixVarFromArray(edk0_trueK0Dir, truePseudoDirArr, 3);
            Float_t noRecoDauVtx[3] = {-999.f, -999.f, -999.f};
            output.FillMatrixVarFromArray(edk0_trueDecayVtxFromRecoDaughters, noRecoDauVtx, 3);
            output.FillVectorVar(edk0_truePDG, _trueK0_PDG[idx]);
            output.FillVectorVar(edk0_trueProcessEnd, _trueK0_processEnd[idx]);
            output.FillVectorVar(edk0_trueParentPDG, _trueK0_parentPDG[idx]);
            output.FillMatrixVarFromArray(edk0_trueParentStartPos, _trueK0_parentStartPos[idx], 3);
            output.FillMatrixVarFromArray(edk0_trueParentEndPos, _trueK0_parentEndPos[idx], 3);
            output.FillVectorVar(edk0_trueNDaughters, _trueK0_nDaughters[idx]);
            output.FillMatrixVarFromArray(edk0_trueDaughterStartPos, _trueK0_daughterStartPos[idx], kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueDaughterEndPos, _trueK0_daughterEndPos[idx], kMaxTrueDaughters*3);
            output.FillMatrixVarFromArray(edk0_trueDaughterPDG, _trueK0_daughterPDG[idx], kMaxTrueDaughters);
            output.FillMatrixVarFromArray(edk0_trueSiblingStartPos, _trueK0_siblingStartPos[idx], kMaxTrueSiblings*3);
            output.FillMatrixVarFromArray(edk0_trueSiblingEndPos, _trueK0_siblingEndPos[idx], kMaxTrueSiblings*3);
            output.FillMatrixVarFromArray(edk0_trueSiblingPDG, _trueK0_siblingPDG[idx], kMaxTrueSiblings);
            output.IncrementCounter(ednK0Candidates);

            _k0_hasTrueObject[idx] = 1;
            _k0_truePDG[idx] = _trueK0_PDG[idx];
            _k0_trueProcessEnd[idx] = _trueK0_processEnd[idx];
            _k0_trueParentPDG[idx] = _trueK0_parentPDG[idx];
            for (Int_t c = 0; c < 3; ++c) {
                _k0_creationVtxPos[idx][c] = creationPos[c];
                _k0_annihilationVtxPos[idx][c] = annihilationPos[c];
                _k0_startPos[idx][c] = creationPos[c];
                _k0_endPos[idx][c] = annihilationPos[c];
                _k0_trueStartPos[idx][c] = _trueK0_startPos[idx][c];
                _k0_trueEndPos[idx][c] = _trueK0_endPos[idx][c];
                _k0_trueDecayVtxFromRecoDaughters[idx][c] = -999.f;
                _k0_trueParentStartPos[idx][c] = _trueK0_parentStartPos[idx][c];
                _k0_trueParentEndPos[idx][c] = _trueK0_parentEndPos[idx][c];
            }
            _k0_trueNDaughters[idx] = _trueK0_nDaughters[idx];
            for (Int_t c = 0; c < kMaxTrueDaughters*3; ++c) {
                _k0_trueDaughterStartPos[idx][c] = _trueK0_daughterStartPos[idx][c];
                _k0_trueDaughterEndPos[idx][c] = _trueK0_daughterEndPos[idx][c];
            }
            for (Int_t c = 0; c < kMaxTrueDaughters; ++c) {
                _k0_trueDaughterPDG[idx][c] = _trueK0_daughterPDG[idx][c];
            }
            for (Int_t c = 0; c < kMaxTrueSiblings*3; ++c) {
                _k0_trueSiblingStartPos[idx][c] = _trueK0_siblingStartPos[idx][c];
                _k0_trueSiblingEndPos[idx][c] = _trueK0_siblingEndPos[idx][c];
            }
            for (Int_t c = 0; c < kMaxTrueSiblings; ++c) {
                _k0_trueSiblingPDG[idx][c] = _trueK0_siblingPDG[idx][c];
            }
        }
    }
}

//********************************************************************
bool neutralKaonEventDisplay::ReadAnalysisData(TTree* tree) {
//********************************************************************
    // Set branch addresses for K0 candidate data
    tree->SetBranchAddress("ED_nK0Candidates", &_nK0Candidates);
    tree->SetBranchAddress("ED_k0_daughter1ID", _k0_daughter1ID);
    tree->SetBranchAddress("ED_k0_daughter2ID", _k0_daughter2ID);
    tree->SetBranchAddress("ED_k0_parentID", _k0_parentID);
    tree->SetBranchAddress("ED_k0_secondParticleID", _k0_secondParticleID);
    tree->SetBranchAddress("ED_k0_creationVtxRadius", _k0_creationVtxRadius);
    tree->SetBranchAddress("ED_k0_annihilationVtxRadius", _k0_annihilationVtxRadius);
    tree->SetBranchAddress("ED_k0_creationVtxDeg", _k0_creationVtxDeg);
    tree->SetBranchAddress("ED_k0_annihilationVtxDeg", _k0_annihilationVtxDeg);
    tree->SetBranchAddress("ED_k0_creationVtxDegDist", _k0_creationVtxDegDist);
    tree->SetBranchAddress("ED_k0_annihilationVtxDegDist", _k0_annihilationVtxDegDist);
    tree->SetBranchAddress("ED_k0_creationVtxPos", _k0_creationVtxPos);
    tree->SetBranchAddress("ED_k0_annihilationVtxPos", _k0_annihilationVtxPos);
    tree->SetBranchAddress("ED_k0_annihilationVtxFitPos", _k0_annihilationVtxFitPos);
    tree->SetBranchAddress("ED_k0_startPos", _k0_startPos);
    tree->SetBranchAddress("ED_k0_endPos", _k0_endPos);
    tree->SetBranchAddress("ED_k0_annVtx_fitLine1Start", _k0_annVtx_fitLine1Start);
    tree->SetBranchAddress("ED_k0_annVtx_fitLine1Dir", _k0_annVtx_fitLine1Dir);
    tree->SetBranchAddress("ED_k0_annVtx_fitLine2Start", _k0_annVtx_fitLine2Start);
    tree->SetBranchAddress("ED_k0_annVtx_fitLine2Dir", _k0_annVtx_fitLine2Dir);
    tree->SetBranchAddress("ED_k0_annVtx_pandoraLine1Dir", _k0_annVtx_pandoraLine1Dir);
    tree->SetBranchAddress("ED_k0_annVtx_pandoraLine2Dir", _k0_annVtx_pandoraLine2Dir);
    tree->SetBranchAddress("ED_k0_annVtx_closestPt1", _k0_annVtx_closestPt1);
    tree->SetBranchAddress("ED_k0_annVtx_closestPt2", _k0_annVtx_closestPt2);
    tree->SetBranchAddress("ED_k0_annVtx_pandoraClosestPt1", _k0_annVtx_pandoraClosestPt1);
    tree->SetBranchAddress("ED_k0_annVtx_pandoraClosestPt2", _k0_annVtx_pandoraClosestPt2);
    tree->SetBranchAddress("ED_k0_annVtx_momentumDirFit", _k0_annVtx_momentumDirFit);
    tree->SetBranchAddress("ED_k0_trueK0Dir", _k0_trueK0Dir);
    if (tree->GetBranch("ED_k0_trueDecayVtxFromRecoDaughters")) {
        tree->SetBranchAddress("ED_k0_trueDecayVtxFromRecoDaughters", _k0_trueDecayVtxFromRecoDaughters);
    }
    tree->SetBranchAddress("ED_k0_creationVtx_fitLineBeamStart", _k0_creationVtx_fitLineBeamStart);
    tree->SetBranchAddress("ED_k0_creationVtx_fitLineBeamDir", _k0_creationVtx_fitLineBeamDir);
    tree->SetBranchAddress("ED_k0_creationVtx_fitLineSecondStart", _k0_creationVtx_fitLineSecondStart);
    tree->SetBranchAddress("ED_k0_creationVtx_fitLineSecondDir", _k0_creationVtx_fitLineSecondDir);
    tree->SetBranchAddress("ED_k0_creationVtx_closestPtBeam", _k0_creationVtx_closestPtBeam);
    tree->SetBranchAddress("ED_k0_creationVtx_closestPtSecond", _k0_creationVtx_closestPtSecond);
    tree->SetBranchAddress("ED_k0_fitLineLength", _k0_fitLineLength);
    tree->SetBranchAddress("ED_k0_hasTrueObject", _k0_hasTrueObject);
    tree->SetBranchAddress("ED_k0_trueStartPos", _k0_trueStartPos);
    tree->SetBranchAddress("ED_k0_trueEndPos", _k0_trueEndPos);
    tree->SetBranchAddress("ED_k0_truePDG", _k0_truePDG);
    tree->SetBranchAddress("ED_k0_trueProcessEnd", _k0_trueProcessEnd);
    tree->SetBranchAddress("ED_k0_trueParentPDG", _k0_trueParentPDG);
    tree->SetBranchAddress("ED_k0_trueParentStartPos", _k0_trueParentStartPos);
    tree->SetBranchAddress("ED_k0_trueParentEndPos", _k0_trueParentEndPos);
    tree->SetBranchAddress("ED_k0_parentStartPos", _k0_parentStartPos);
    tree->SetBranchAddress("ED_k0_parentEndPos", _k0_parentEndPos);
    tree->SetBranchAddress("ED_k0_parentEndPosCorrected", _k0_parentEndPosCorrected);
    if (tree->GetBranch("ED_k0_parentTailFitAnchor")) {
        tree->SetBranchAddress("ED_k0_parentTailFitAnchor", _k0_parentTailFitAnchor);
    }
    if (tree->GetBranch("ED_k0_parentTailFitDir")) {
        tree->SetBranchAddress("ED_k0_parentTailFitDir", _k0_parentTailFitDir);
    }
    if (tree->GetBranch("ED_k0_parentTailFitLength")) {
        tree->SetBranchAddress("ED_k0_parentTailFitLength", _k0_parentTailFitLength);
    }
    tree->SetBranchAddress("ED_k0_parentTailHitsRaw", _k0_parentTailHitsRaw);
    tree->SetBranchAddress("ED_k0_parentTailHitsProjected", _k0_parentTailHitsProjected);
    tree->SetBranchAddress("ED_k0_parentTailHitsRawN", _k0_parentTailHitsRawN);
    tree->SetBranchAddress("ED_k0_parentTailHitsProjectedN", _k0_parentTailHitsProjectedN);
    tree->SetBranchAddress("ED_k0_parentLength", _k0_parentLength);
    tree->SetBranchAddress("ED_k0_parentTrajDir", _k0_parentTrajDir);
    tree->SetBranchAddress("ED_k0_parentTrajDirHist", _k0_parentTrajDirHist);
    tree->SetBranchAddress("ED_k0_parentTrajDirNPts", _k0_parentTrajDirNPts);
    tree->SetBranchAddress("ED_k0_secondTrajDir", _k0_secondTrajDir);
    tree->SetBranchAddress("ED_k0_secondTrajDirNPts", _k0_secondTrajDirNPts);
    tree->SetBranchAddress("ED_k0_dau1TrajDir", _k0_dau1TrajDir);
    tree->SetBranchAddress("ED_k0_dau1TrajDirNPts", _k0_dau1TrajDirNPts);
    tree->SetBranchAddress("ED_k0_dau2TrajDir", _k0_dau2TrajDir);
    tree->SetBranchAddress("ED_k0_dau2TrajDirNPts", _k0_dau2TrajDirNPts);
    tree->SetBranchAddress("ED_k0_trueNDaughters", _k0_trueNDaughters);
    tree->SetBranchAddress("ED_k0_trueDaughterStartPos", _k0_trueDaughterStartPos);
    tree->SetBranchAddress("ED_k0_trueDaughterEndPos", _k0_trueDaughterEndPos);
    tree->SetBranchAddress("ED_k0_trueDaughterPDG", _k0_trueDaughterPDG);
    tree->SetBranchAddress("ED_k0_trueParentNDaughters", _k0_trueParentNDaughters);
    tree->SetBranchAddress("ED_k0_trueParentDaughterStartPos", _k0_trueParentDaughterStartPos);
    tree->SetBranchAddress("ED_k0_trueParentDaughterEndPos", _k0_trueParentDaughterEndPos);
    tree->SetBranchAddress("ED_k0_trueParentDaughterPDG", _k0_trueParentDaughterPDG);
    tree->SetBranchAddress("ED_k0_trueNSiblings", _k0_trueNSiblings);
    tree->SetBranchAddress("ED_k0_trueSiblingStartPos", _k0_trueSiblingStartPos);
    tree->SetBranchAddress("ED_k0_trueSiblingEndPos", _k0_trueSiblingEndPos);
    tree->SetBranchAddress("ED_k0_trueSiblingPDG", _k0_trueSiblingPDG);

    tree->SetBranchAddress("ED_nTrueK0", &_nTrueK0);
    tree->SetBranchAddress("ED_trueK0_startPos", _trueK0_startPos);
    tree->SetBranchAddress("ED_trueK0_endPos", _trueK0_endPos);
    tree->SetBranchAddress("ED_trueK0_PDG", _trueK0_PDG);
    tree->SetBranchAddress("ED_trueK0_processEnd", _trueK0_processEnd);
    tree->SetBranchAddress("ED_trueK0_parentPDG", _trueK0_parentPDG);
    tree->SetBranchAddress("ED_trueK0_parentStartPos", _trueK0_parentStartPos);
    tree->SetBranchAddress("ED_trueK0_parentEndPos", _trueK0_parentEndPos);
    tree->SetBranchAddress("ED_trueK0_nDaughters", _trueK0_nDaughters);
    tree->SetBranchAddress("ED_trueK0_daughterStartPos", _trueK0_daughterStartPos);
    tree->SetBranchAddress("ED_trueK0_daughterEndPos", _trueK0_daughterEndPos);
    tree->SetBranchAddress("ED_trueK0_daughterPDG", _trueK0_daughterPDG);
    tree->SetBranchAddress("ED_trueK0_nSiblings", _trueK0_nSiblings);
    tree->SetBranchAddress("ED_trueK0_siblingStartPos", _trueK0_siblingStartPos);
    tree->SetBranchAddress("ED_trueK0_siblingEndPos", _trueK0_siblingEndPos);
    tree->SetBranchAddress("ED_trueK0_siblingPDG", _trueK0_siblingPDG);

    tree->SetBranchAddress("ED_nAllTrueParticles", &_nAllTrueParticles);
    tree->SetBranchAddress("ED_allTrueParticle_startPos", _allTrueParticle_startPos);
    tree->SetBranchAddress("ED_allTrueParticle_endPos", _allTrueParticle_endPos);
    tree->SetBranchAddress("ED_allTrueParticle_PDG", _allTrueParticle_PDG);
    tree->SetBranchAddress("ED_allTrueParticle_processEnd", _allTrueParticle_processEnd);

    // Load only analysis branches for the already-selected entry.
    // This avoids a full second TTree::GetEntry() pass while still refreshing
    // K0 display payload after branch addresses are bound.
    const Long64_t entry = tree->GetReadEntry();
    if (entry >= 0) {
        const std::initializer_list<const char*> analysisBranches = {
            "ED_nK0Candidates",
            "ED_k0_daughter1ID", "ED_k0_daughter2ID", "ED_k0_parentID", "ED_k0_secondParticleID",
            "ED_k0_creationVtxRadius", "ED_k0_annihilationVtxRadius",
            "ED_k0_creationVtxDeg", "ED_k0_annihilationVtxDeg",
            "ED_k0_creationVtxDegDist", "ED_k0_annihilationVtxDegDist",
            "ED_k0_creationVtxPos", "ED_k0_annihilationVtxPos", "ED_k0_annihilationVtxFitPos",
            "ED_k0_startPos", "ED_k0_endPos",
            "ED_k0_annVtx_fitLine1Start", "ED_k0_annVtx_fitLine1Dir",
            "ED_k0_annVtx_fitLine2Start", "ED_k0_annVtx_fitLine2Dir",
            "ED_k0_annVtx_pandoraLine1Dir", "ED_k0_annVtx_pandoraLine2Dir",
            "ED_k0_annVtx_closestPt1", "ED_k0_annVtx_closestPt2",
            "ED_k0_annVtx_pandoraClosestPt1", "ED_k0_annVtx_pandoraClosestPt2",
            "ED_k0_annVtx_momentumDirFit", "ED_k0_trueK0Dir",
            "ED_k0_trueDecayVtxFromRecoDaughters",
            "ED_k0_creationVtx_fitLineBeamStart", "ED_k0_creationVtx_fitLineBeamDir",
            "ED_k0_creationVtx_fitLineSecondStart", "ED_k0_creationVtx_fitLineSecondDir",
            "ED_k0_creationVtx_closestPtBeam", "ED_k0_creationVtx_closestPtSecond",
            "ED_k0_fitLineLength", "ED_k0_hasTrueObject",
            "ED_k0_trueStartPos", "ED_k0_trueEndPos", "ED_k0_truePDG", "ED_k0_trueProcessEnd",
            "ED_k0_trueParentPDG", "ED_k0_trueParentStartPos", "ED_k0_trueParentEndPos",
            "ED_k0_trueParentNDaughters", "ED_k0_trueParentDaughterStartPos",
            "ED_k0_trueParentDaughterEndPos", "ED_k0_trueParentDaughterPDG",
            "ED_k0_parentStartPos", "ED_k0_parentEndPos", "ED_k0_parentLength",
            "ED_k0_parentEndPosCorrected", "ED_k0_parentTailFitAnchor", "ED_k0_parentTailFitDir",
            "ED_k0_parentTailFitLength", "ED_k0_parentTailHitsRaw", "ED_k0_parentTailHitsProjected",
            "ED_k0_parentTailHitsRawN", "ED_k0_parentTailHitsProjectedN",
            "ED_k0_parentTrajDir", "ED_k0_parentTrajDirHist", "ED_k0_parentTrajDirNPts",
            "ED_k0_secondTrajDir", "ED_k0_secondTrajDirNPts",
            "ED_k0_dau1TrajDir", "ED_k0_dau1TrajDirNPts",
            "ED_k0_dau2TrajDir", "ED_k0_dau2TrajDirNPts",
            "ED_k0_trueNDaughters", "ED_k0_trueDaughterStartPos",
            "ED_k0_trueDaughterEndPos", "ED_k0_trueDaughterPDG",
            "ED_k0_trueNSiblings", "ED_k0_trueSiblingStartPos",
            "ED_k0_trueSiblingEndPos", "ED_k0_trueSiblingPDG",
            "ED_nTrueK0", "ED_trueK0_startPos", "ED_trueK0_endPos",
            "ED_trueK0_PDG", "ED_trueK0_processEnd",
            "ED_trueK0_parentPDG", "ED_trueK0_parentStartPos", "ED_trueK0_parentEndPos",
            "ED_trueK0_nDaughters", "ED_trueK0_daughterStartPos",
            "ED_trueK0_daughterEndPos", "ED_trueK0_daughterPDG",
            "ED_trueK0_nSiblings", "ED_trueK0_siblingStartPos",
            "ED_trueK0_siblingEndPos", "ED_trueK0_siblingPDG",
            "ED_nAllTrueParticles", "ED_allTrueParticle_startPos", "ED_allTrueParticle_endPos",
            "ED_allTrueParticle_PDG", "ED_allTrueParticle_processEnd"
        };
        for (const char* bname : analysisBranches) {
            if (TBranch* br = tree->GetBranch(bname)) {
                br->GetEntry(entry);
            }
        }
    }

    std::cout << "Read K0 analysis data: " << _nK0Candidates << " K0 candidates" << std::endl;

    return true;
}

//********************************************************************
void neutralKaonEventDisplay::DrawAnalysisContent3D(TEveScene* scene) {
//********************************************************************
    EnsureSelectionHooks();
    EnsureGroupCheckboxSync();

    // Draw K0 candidates
    TEveElementList* neutralParticlesGroup = PrepareGroup(scene, "Neutral Particles");
    TEveElementList* trueParticlesGroup = PrepareGroup(scene, "True Particles");

    auto addElement = [&](TEveElementList* group, TEveElement* element) {
        AppendElementToGroup(scene, group, element);
    };

    auto anchorPoint = [&](TEveElement* element, Float_t x, Float_t y, Float_t z) {
        if (!element) return;
        if (x <= -900.f || y <= -900.f || z <= -900.f) return;
        this->RegisterMeasurementAnchor(element, x, y, z);
    };

    auto drawArrow3D = [&](TEveElementList* group,
                           const char* name,
                           const TVector3& start,
                           const TVector3& dirUnit,
                           double length,
                           Int_t color,
                           Int_t style) {
        if (!group || dirUnit.Mag2() <= 1e-10 || length <= 0.0) return;
        const TVector3 u = dirUnit.Unit();
        const TVector3 tip = start + length * u;

        TEveElementList* arrowRoot = new TEveElementList(name);
        addElement(group, arrowRoot);

        TEveLine* shaft = new TEveLine("shaft");
        shaft->SetPoint(0, start.X(), start.Y(), start.Z());
        shaft->SetPoint(1, tip.X(), tip.Y(), tip.Z());
        shaft->SetMainColor(color);
        shaft->SetLineWidth(2);
        shaft->SetLineStyle(style);
        addElement(arrowRoot, shaft);

        // Build a stable perpendicular basis for arrowhead wings.
        TVector3 ref(0.0, 0.0, 1.0);
        if (std::abs(u.Dot(ref)) > 0.9) ref.SetXYZ(0.0, 1.0, 0.0);
        TVector3 n1 = u.Cross(ref).Unit();
        TVector3 n2 = u.Cross(n1).Unit();
        const double headLen = std::max(3.0, 0.10 * length);
        const double wing = 0.32 * headLen;
        const TVector3 base = tip - headLen * u;
        const TVector3 wing1 = base + wing * n1;
        const TVector3 wing2 = base - wing * n1;
        const TVector3 wing3 = base + wing * n2;
        const TVector3 wing4 = base - wing * n2;

        auto addWing = [&](const TVector3& w, const char* wname) {
            TEveLine* wingLine = new TEveLine(wname);
            wingLine->SetPoint(0, tip.X(), tip.Y(), tip.Z());
            wingLine->SetPoint(1, w.X(), w.Y(), w.Z());
            wingLine->SetMainColor(color);
            wingLine->SetLineWidth(2);
            wingLine->SetLineStyle(style);
            addElement(arrowRoot, wingLine);
        };
        addWing(wing1, "head wing A");
        addWing(wing2, "head wing B");
        addWing(wing3, "head wing C");
        addWing(wing4, "head wing D");
    };

    _parentDirElementToIndex.clear();
    _parentDirectionLines.clear();

    for (Int_t i = 0; i < _nK0Candidates && i < kMaxK0; i++) {
        Int_t parentColor = kBlue;
        TEveElementList* neutralParticleGroup =
            new TEveElementList(Form("Neutral UID=%d TruePDG=%d", i, _k0_truePDG[i]));
        addElement(neutralParticlesGroup, neutralParticleGroup);

        TEveElementList* trueParticleGroup =
            new TEveElementList(Form("True UID=%d PDG=%d", i, _k0_truePDG[i]));
        addElement(trueParticlesGroup, trueParticleGroup);

        TEveElementList* creationVertexGroup = neutralParticleGroup;
        TEveElementList* annihilationVertexGroup = neutralParticleGroup;
        TEveElementList* momentumArrowGroup = neutralParticleGroup;
        TEveElementList* creationFitGroup = neutralParticleGroup;
        TEveElementList* annihilationFitGroup = neutralParticleGroup;

        TEveElementList* truthGroup = trueParticleGroup;
        TEveElementList* truthParentGroup = trueParticleGroup;
        TEveElementList* truthDaughterGroup = trueParticleGroup;
        TEveElementList* truthSiblingGroup = trueParticleGroup;
        TEveElementList* annotationGroup = trueParticleGroup;

        if (_k0_startPos[i][0] > -900 && _k0_endPos[i][0] > -900) {
            Int_t recoColor = kGray + 1;
            if (_k0_hasTrueObject[i] && _k0_truePDG[i] != -999) {
                recoColor = GetParticleColor(_k0_truePDG[i]);
                if (recoColor == kBlack) recoColor = kGray + 1;
            }

            TEveLine* recoK0 = new TEveLine(Form("K0 #%d Reco Trajectory", i));
            recoK0->SetPoint(0, _k0_startPos[i][0], _k0_startPos[i][1], _k0_startPos[i][2]);
            recoK0->SetPoint(1, _k0_endPos[i][0], _k0_endPos[i][1], _k0_endPos[i][2]);
            recoK0->SetMainColor(recoColor);
            recoK0->SetLineWidth(3);
            recoK0->SetLineStyle(2);
            addElement(neutralParticleGroup, recoK0);
            anchorPoint(recoK0, _k0_startPos[i][0], _k0_startPos[i][1], _k0_startPos[i][2]);
            anchorPoint(recoK0, _k0_endPos[i][0], _k0_endPos[i][1], _k0_endPos[i][2]);
        }

        // Creation vertex sphere
        Float_t creationX = _k0_creationVtxPos[i][0];
        Float_t creationY = _k0_creationVtxPos[i][1];
        Float_t creationZ = _k0_creationVtxPos[i][2];

        if (creationX > -900 && creationY > -900 && creationZ > -900) {
            Float_t radius = _k0_creationVtxRadius[i];
            TGeoSphere* creationSphere = new TGeoSphere(0, radius);
            TGeoVolume* creationVol = new TGeoVolume(Form("CreationSphere_%d", i), creationSphere);
            TEveGeoShape* creationShape = new TEveGeoShape(Form("K0 #%d Creation Radius (%.0f cm)", i, radius));
            creationShape->SetShape(creationVol->GetShape());
            creationShape->SetMainColor(parentColor);
            creationShape->SetMainTransparency(80);
            TEveTrans& creationTrans = creationShape->RefMainTrans();
            creationTrans.SetPos(creationX, creationY, creationZ);
            addElement(creationVertexGroup, creationShape);
            anchorPoint(creationShape, creationX, creationY, creationZ);

            // Creation vertex marker
            TEvePointSet* vtx = new TEvePointSet(Form("K0 #%d Creation Vertex", i));
            vtx->SetNextPoint(creationX, creationY, creationZ);
            vtx->SetMarkerStyle(29); // Star
            vtx->SetMarkerSize(3.0);
            vtx->SetMainColor(parentColor);
            addElement(creationVertexGroup, vtx);
            anchorPoint(vtx, creationX, creationY, creationZ);

            // Creation vertex degeneracy label
            Int_t creationDeg = _k0_creationVtxDeg[i];
            if (creationDeg > -1) {
                TEveText* degLabel = new TEveText(Form("deg=%d", creationDeg));
                degLabel->SetMainColor(parentColor);
                degLabel->SetFontSize(18);
                degLabel->RefMainTrans().SetPos(creationX + radius + 5, creationY + radius + 5, creationZ);
                addElement(creationVertexGroup, degLabel);
            }
        }

        Float_t parentEndX = _k0_parentEndPos[i][0];
        Float_t parentEndY = _k0_parentEndPos[i][1];
        Float_t parentEndZ = _k0_parentEndPos[i][2];
        Float_t parentDirX = _k0_parentTrajDir[i][0];
        Float_t parentDirY = _k0_parentTrajDir[i][1];
        Float_t parentDirZ = _k0_parentTrajDir[i][2];

        const Float_t parentStartX = _k0_parentStartPos[i][0];
        const Float_t parentStartY = _k0_parentStartPos[i][1];
        const Float_t parentStartZ = _k0_parentStartPos[i][2];
        const Bool_t parentStartValid = parentStartX > -900 && parentStartY > -900 && parentStartZ > -900;
        const Bool_t parentEndValid = parentEndX > -900 && parentEndY > -900 && parentEndZ > -900;
        TVector3 dir(parentDirX, parentDirY, parentDirZ);
        Bool_t dirValid = parentDirX > -900 && parentDirY > -900 && parentDirZ > -900 && dir.Mag() > 0;
        if (!dirValid && parentStartValid && parentEndValid) {
            dir = TVector3(parentEndX - parentStartX, parentEndY - parentStartY, parentEndZ - parentStartZ);
            if (dir.Mag() > 0) {
                dir = dir.Unit();
                dirValid = kTRUE;
            }
        } else if (dirValid) {
            dir = dir.Unit();
        }

        Double_t parentLength = _k0_parentLength[i];
        if (parentLength <= 0 && parentStartValid && parentEndValid) {
            parentLength = TVector3(parentEndX - parentStartX, parentEndY - parentStartY, parentEndZ - parentStartZ).Mag();
        }

        if (!dirValid || parentLength <= 0) {
            // Cannot draw reliable parent direction line
        } else {
            TVector3 startPoint;
            TVector3 endPoint;
            if (parentStartValid) {
                startPoint.SetXYZ(parentStartX, parentStartY, parentStartZ);
            } else if (parentEndValid) {
                startPoint.SetXYZ(parentEndX, parentEndY, parentEndZ);
                startPoint -= dir * parentLength;
            } else {
                continue;
            }

            if (parentEndValid) {
                endPoint.SetXYZ(parentEndX, parentEndY, parentEndZ);
            } else {
                endPoint = startPoint + dir * parentLength;
            }

            TVector3 extendedStart = startPoint - dir * (0.2 * parentLength);
            TVector3 extendedEnd = endPoint + dir * (0.2 * parentLength);

            TEveLine* parentLine = new TEveLine(Form("Parent Direction #%d", i));
            parentLine->SetLineWidth(3);
            Int_t parentPDG = _k0_trueParentPDG[i];
            parentLine->SetMainColor(GetParticleColor(parentPDG));
            parentLine->SetPickable(kTRUE);
            parentLine->SetElementTitle(Form("Parent Direction MPV [%0.3f, %0.3f, %0.3f] (%d points)",
                                             parentDirX, parentDirY, parentDirZ,
                                             _k0_parentTrajDirNPts[i]));
            parentLine->SetNextPoint(extendedStart.X(), extendedStart.Y(), extendedStart.Z());
            parentLine->SetNextPoint(extendedEnd.X(), extendedEnd.Y(), extendedEnd.Z());
            addElement(momentumArrowGroup, parentLine);
            anchorPoint(parentLine, extendedStart.X(), extendedStart.Y(), extendedStart.Z());
            anchorPoint(parentLine, extendedEnd.X(), extendedEnd.Y(), extendedEnd.Z());

            _parentDirElementToIndex[parentLine] = i;
            _parentDirectionLines.push_back(parentLine);
        }

        const Int_t rawTailN = std::max(0, std::min(_k0_parentTailHitsRawN[i], kMaxBeamTailHits));
        const Int_t projectedTailN = std::max(0, std::min(_k0_parentTailHitsProjectedN[i], kMaxBeamTailHits));

        if (rawTailN > 0) {
            TEvePointSet* rawTail = new TEvePointSet(Form("K0 #%d Parent Tail Raw", i));
            rawTail->SetMarkerStyle(20);
            rawTail->SetMarkerSize(1.4);
            rawTail->SetMainColor(kGray + 2);
            for (Int_t t = 0; t < rawTailN; ++t) {
                const Float_t x = _k0_parentTailHitsRaw[i][t * 3 + 0];
                const Float_t y = _k0_parentTailHitsRaw[i][t * 3 + 1];
                const Float_t z = _k0_parentTailHitsRaw[i][t * 3 + 2];
                if (x <= -900.f || y <= -900.f || z <= -900.f) continue;
                rawTail->SetNextPoint(x, y, z);
                anchorPoint(rawTail, x, y, z);
            }
            addElement(creationFitGroup, rawTail);
        }

        if (projectedTailN > 0) {
            TEvePointSet* projectedTail = new TEvePointSet(Form("K0 #%d Parent Tail Projected", i));
            projectedTail->SetMarkerStyle(24);
            projectedTail->SetMarkerSize(1.6);
            projectedTail->SetMainColor(kYellow + 1);
            for (Int_t t = 0; t < projectedTailN; ++t) {
                const Float_t x = _k0_parentTailHitsProjected[i][t * 3 + 0];
                const Float_t y = _k0_parentTailHitsProjected[i][t * 3 + 1];
                const Float_t z = _k0_parentTailHitsProjected[i][t * 3 + 2];
                if (x <= -900.f || y <= -900.f || z <= -900.f) continue;
                projectedTail->SetNextPoint(x, y, z);
                anchorPoint(projectedTail, x, y, z);
            }
            addElement(creationFitGroup, projectedTail);
        }

        if (parentEndValid) {
            TEvePointSet* rawEnd = new TEvePointSet(Form("K0 #%d Parent End (Raw)", i));
            rawEnd->SetNextPoint(parentEndX, parentEndY, parentEndZ);
            rawEnd->SetMarkerStyle(34);
            rawEnd->SetMarkerSize(2.2);
            rawEnd->SetMainColor(kMagenta + 1);
            addElement(creationFitGroup, rawEnd);
            anchorPoint(rawEnd, parentEndX, parentEndY, parentEndZ);
        }

        const Float_t correctedEndX = _k0_parentEndPosCorrected[i][0];
        const Float_t correctedEndY = _k0_parentEndPosCorrected[i][1];
        const Float_t correctedEndZ = _k0_parentEndPosCorrected[i][2];
        const Bool_t correctedEndValid = correctedEndX > -900.f && correctedEndY > -900.f && correctedEndZ > -900.f;
        if (correctedEndValid) {
            TEvePointSet* correctedEnd = new TEvePointSet(Form("K0 #%d Parent End (Projected)", i));
            correctedEnd->SetNextPoint(correctedEndX, correctedEndY, correctedEndZ);
            correctedEnd->SetMarkerStyle(29);
            correctedEnd->SetMarkerSize(2.6);
            correctedEnd->SetMainColor(kYellow + 1);
            addElement(creationFitGroup, correctedEnd);
            anchorPoint(correctedEnd, correctedEndX, correctedEndY, correctedEndZ);
        }

        const Float_t fitAnchorX = _k0_parentTailFitAnchor[i][0];
        const Float_t fitAnchorY = _k0_parentTailFitAnchor[i][1];
        const Float_t fitAnchorZ = _k0_parentTailFitAnchor[i][2];
        const Float_t fitDirX = _k0_parentTailFitDir[i][0];
        const Float_t fitDirY = _k0_parentTailFitDir[i][1];
        const Float_t fitDirZ = _k0_parentTailFitDir[i][2];
        const Float_t fitLen = _k0_parentTailFitLength[i];
        const Bool_t fitAnchorValid = fitAnchorX > -900.f && fitAnchorY > -900.f && fitAnchorZ > -900.f;
        TVector3 fitDir(fitDirX, fitDirY, fitDirZ);
        if (fitAnchorValid && fitDir.Mag2() > 1e-10 && fitLen > 0.f) {
            fitDir = fitDir.Unit();
            TVector3 anchor(fitAnchorX, fitAnchorY, fitAnchorZ);
            const double fitSpan = std::max<double>(std::max<double>(fitLen, parentLength), 20.0);
            const TVector3 lineStart = anchor - fitSpan * fitDir;
            const TVector3 lineEnd = anchor + fitSpan * fitDir;

            TEveLine* fitLine = new TEveLine(Form("K0 #%d Parent Tail Backward Fit Line", i));
            fitLine->SetPoint(0, lineStart.X(), lineStart.Y(), lineStart.Z());
            fitLine->SetPoint(1, lineEnd.X(), lineEnd.Y(), lineEnd.Z());
            fitLine->SetMainColor(kCyan + 1);
            fitLine->SetLineStyle(2);
            fitLine->SetLineWidth(3);
            addElement(creationFitGroup, fitLine);
            anchorPoint(fitLine, lineStart.X(), lineStart.Y(), lineStart.Z());
            anchorPoint(fitLine, lineEnd.X(), lineEnd.Y(), lineEnd.Z());

            TEvePointSet* fitAnchor = new TEvePointSet(Form("K0 #%d Parent Tail Backward-Fit Anchor", i));
            fitAnchor->SetNextPoint(anchor.X(), anchor.Y(), anchor.Z());
            fitAnchor->SetMarkerStyle(43);
            fitAnchor->SetMarkerSize(2.2);
            fitAnchor->SetMainColor(kCyan + 1);
            addElement(creationFitGroup, fitAnchor);
            anchorPoint(fitAnchor, anchor.X(), anchor.Y(), anchor.Z());
        }

        // Annihilation vertex sphere
        Float_t annihilationX = _k0_annihilationVtxPos[i][0];
        Float_t annihilationY = _k0_annihilationVtxPos[i][1];
        Float_t annihilationZ = _k0_annihilationVtxPos[i][2];

        if (annihilationX > -900 && annihilationY > -900 && annihilationZ > -900) {
            Float_t annihilationRadius = _k0_annihilationVtxRadius[i];
            TGeoSphere* annihilationSphere = new TGeoSphere(0, annihilationRadius);
            TGeoVolume* annihilationVol = new TGeoVolume(Form("AnnihilationSphere_%d", i), annihilationSphere);
            TEveGeoShape* annihilationShape = new TEveGeoShape(Form("K0 #%d Annihilation Radius (%.0f cm)", i, annihilationRadius));
            annihilationShape->SetShape(annihilationVol->GetShape());
            annihilationShape->SetMainColor(kRed);
            annihilationShape->SetMainTransparency(80);
            TEveTrans& annihilationTrans = annihilationShape->RefMainTrans();
            annihilationTrans.SetPos(annihilationX, annihilationY, annihilationZ);
            addElement(annihilationVertexGroup, annihilationShape);
            anchorPoint(annihilationShape, annihilationX, annihilationY, annihilationZ);

            // Annihilation vertex marker
            TEvePointSet* vtx = new TEvePointSet(Form("K0 #%d Decay Vertex", i));
            vtx->SetNextPoint(annihilationX, annihilationY, annihilationZ);
            vtx->SetMarkerStyle(29); // Star
            vtx->SetMarkerSize(3.0);
            vtx->SetMainColor(kRed);
            addElement(annihilationVertexGroup, vtx);
            anchorPoint(vtx, annihilationX, annihilationY, annihilationZ);

            // Fit-based annihilation vertex marker
            if (_k0_annihilationVtxFitPos[i][0] > -900 &&
                _k0_annihilationVtxFitPos[i][1] > -900 &&
                _k0_annihilationVtxFitPos[i][2] > -900) {
                TEvePointSet* fitVtx = new TEvePointSet(Form("K0 #%d Decay Vertex (Fit)", i));
                fitVtx->SetNextPoint(_k0_annihilationVtxFitPos[i][0],
                                     _k0_annihilationVtxFitPos[i][1],
                                     _k0_annihilationVtxFitPos[i][2]);
                fitVtx->SetMarkerStyle(20); // Full circle
                fitVtx->SetMarkerSize(2.2);
                fitVtx->SetMainColor(kOrange + 1);
                addElement(annihilationVertexGroup, fitVtx);
                anchorPoint(fitVtx, _k0_annihilationVtxFitPos[i][0], _k0_annihilationVtxFitPos[i][1], _k0_annihilationVtxFitPos[i][2]);
            }

            // Anchor-hit markers used by the fit-line method.
            if (_k0_annVtx_fitLine1Start[i][0] > -900 &&
                _k0_annVtx_fitLine1Start[i][1] > -900 &&
                _k0_annVtx_fitLine1Start[i][2] > -900) {
                TEvePointSet* dau1Start = new TEvePointSet(Form("K0 #%d Annihilation Daughter1 Fit Anchor", i));
                dau1Start->SetNextPoint(_k0_annVtx_fitLine1Start[i][0],
                                        _k0_annVtx_fitLine1Start[i][1],
                                        _k0_annVtx_fitLine1Start[i][2]);
                dau1Start->SetMarkerStyle(29);
                dau1Start->SetMarkerSize(2.6);
                dau1Start->SetMainColor(kYellow + 1);
                addElement(annihilationVertexGroup, dau1Start);
                anchorPoint(dau1Start, _k0_annVtx_fitLine1Start[i][0], _k0_annVtx_fitLine1Start[i][1], _k0_annVtx_fitLine1Start[i][2]);
            }
            if (_k0_annVtx_fitLine2Start[i][0] > -900 &&
                _k0_annVtx_fitLine2Start[i][1] > -900 &&
                _k0_annVtx_fitLine2Start[i][2] > -900) {
                TEvePointSet* dau2Start = new TEvePointSet(Form("K0 #%d Annihilation Daughter2 Fit Anchor", i));
                dau2Start->SetNextPoint(_k0_annVtx_fitLine2Start[i][0],
                                        _k0_annVtx_fitLine2Start[i][1],
                                        _k0_annVtx_fitLine2Start[i][2]);
                dau2Start->SetMarkerStyle(29);
                dau2Start->SetMarkerSize(2.6);
                dau2Start->SetMainColor(kYellow + 1);
                addElement(annihilationVertexGroup, dau2Start);
                anchorPoint(dau2Start, _k0_annVtx_fitLine2Start[i][0], _k0_annVtx_fitLine2Start[i][1], _k0_annVtx_fitLine2Start[i][2]);
            }

            // Annihilation vertex degeneracy label
            Int_t annihilationDeg = _k0_annihilationVtxDeg[i];
            if (annihilationDeg > 0) {
                TEveText* degLabel = new TEveText(Form("deg=%d", annihilationDeg));
                degLabel->SetMainColor(kRed);
                degLabel->SetFontSize(18);
                degLabel->RefMainTrans().SetPos(annihilationX + annihilationRadius + 5,
                                                annihilationY + annihilationRadius + 5,
                                                annihilationZ);
                addElement(annihilationVertexGroup, degLabel);
            }

            // Reco arrow: Pandora annihilation vertex. True arrow: MC decay from daughter start positions.
            TVector3 recoVtxPos(annihilationX, annihilationY, annihilationZ);
            const double vecLen = std::max(12.0, std::min(38.0, static_cast<double>(_k0_annihilationVtxRadius[i] * 0.5)));
            if (_k0_annVtx_momentumDirFit[i][0] > -900) {
                TVector3 recoDir(_k0_annVtx_momentumDirFit[i][0], _k0_annVtx_momentumDirFit[i][1], _k0_annVtx_momentumDirFit[i][2]);
                if (recoDir.Mag2() > 1e-10) {
                    drawArrow3D(momentumArrowGroup,
                                Form("K0 #%d Vertex Momentum Arrow (Reco)", i),
                                recoVtxPos, recoDir, vecLen, kMagenta + 1, 1);
                }
            }
            if (_k0_trueK0Dir[i][0] > -900 && _k0_hasTrueObject[i]) {
                Float_t trueMcVtx[3];
                bool haveTrueMcVtx = false;
                if (_k0_trueDecayVtxFromRecoDaughters[i][0] > -900.f && _k0_trueDecayVtxFromRecoDaughters[i][1] > -900.f &&
                    _k0_trueDecayVtxFromRecoDaughters[i][2] > -900.f) {
                    trueMcVtx[0] = _k0_trueDecayVtxFromRecoDaughters[i][0];
                    trueMcVtx[1] = _k0_trueDecayVtxFromRecoDaughters[i][1];
                    trueMcVtx[2] = _k0_trueDecayVtxFromRecoDaughters[i][2];
                    haveTrueMcVtx = true;
                }
                if (!haveTrueMcVtx) {
                    haveTrueMcVtx = TrueK0McDecayVertex(trueMcVtx, _k0_trueNDaughters[i], _k0_trueDaughterPDG[i],
                                                        _k0_trueDaughterStartPos[i]);
                }
                if (!haveTrueMcVtx && _k0_trueEndPos[i][0] > -900 && _k0_trueEndPos[i][1] > -900 &&
                    _k0_trueEndPos[i][2] > -900) {
                    trueMcVtx[0] = _k0_trueEndPos[i][0];
                    trueMcVtx[1] = _k0_trueEndPos[i][1];
                    trueMcVtx[2] = _k0_trueEndPos[i][2];
                    haveTrueMcVtx = true;
                }
                if (haveTrueMcVtx) {
                    TVector3 trueDir(_k0_trueK0Dir[i][0], _k0_trueK0Dir[i][1], _k0_trueK0Dir[i][2]);
                    if (trueDir.Mag2() > 1e-10) {
                        TVector3 trueDecayPos(trueMcVtx[0], trueMcVtx[1], trueMcVtx[2]);
                        drawArrow3D(momentumArrowGroup,
                                    Form("K0 #%d Vertex Momentum Arrow (True K0)", i),
                                    trueDecayPos, trueDir, vecLen, kBlue + 1, 2);
                    }
                }
            }
        }

        // K0 true trajectory (if TrueObject exists - solid line, K0 PDG color)
        if (_k0_hasTrueObject[i] && _k0_trueStartPos[i][0] > -900 && _k0_trueEndPos[i][0] > -900) {
            Int_t truePDG = _k0_truePDG[i];
            Int_t trueColor = GetParticleColor(truePDG);
            if (trueColor == kBlack) {
                if (abs(truePDG) == 311 || abs(truePDG) == 310 || abs(truePDG) == 130) {
                    trueColor = GetParticleColor(311);
                } else {
                    trueColor = kGray + 1; // neutral fallback for unknown PDG
                }
            }

            TEveLine* trueK0 = new TEveLine(Form("K0 #%d True Trajectory", i));
            trueK0->SetPoint(0, _k0_trueStartPos[i][0], _k0_trueStartPos[i][1], _k0_trueStartPos[i][2]);
            trueK0->SetPoint(1, _k0_trueEndPos[i][0], _k0_trueEndPos[i][1], _k0_trueEndPos[i][2]);
            trueK0->SetMainColor(trueColor);
            trueK0->SetLineWidth(3);
            trueK0->SetLineStyle(1); // Solid line for truth
            addElement(truthGroup, trueK0);
            anchorPoint(trueK0,
                        _k0_trueStartPos[i][0],
                        _k0_trueStartPos[i][1],
                        _k0_trueStartPos[i][2]);
            anchorPoint(trueK0,
                        _k0_trueEndPos[i][0],
                        _k0_trueEndPos[i][1],
                        _k0_trueEndPos[i][2]);

            // Add markers at true start and end positions
            TEvePointSet* trueStart = new TEvePointSet(Form("K0 #%d True Start", i));
            trueStart->SetNextPoint(_k0_trueStartPos[i][0], _k0_trueStartPos[i][1], _k0_trueStartPos[i][2]);
            trueStart->SetMarkerStyle(29); // Star
            trueStart->SetMarkerSize(2.5);
            trueStart->SetMainColor(trueColor);
            addElement(truthGroup, trueStart);
            anchorPoint(trueStart,
                        _k0_trueStartPos[i][0],
                        _k0_trueStartPos[i][1],
                        _k0_trueStartPos[i][2]);

            TEvePointSet* trueEnd = new TEvePointSet(Form("K0 #%d True End", i));
            trueEnd->SetNextPoint(_k0_trueEndPos[i][0], _k0_trueEndPos[i][1], _k0_trueEndPos[i][2]);
            trueEnd->SetMarkerStyle(29); // Star
            trueEnd->SetMarkerSize(2.5);
            trueEnd->SetMainColor(trueColor);
            addElement(truthGroup, trueEnd);
            anchorPoint(trueEnd,
                        _k0_trueEndPos[i][0],
                        _k0_trueEndPos[i][1],
                        _k0_trueEndPos[i][2]);

            std::string procLabel = ProcessEnumToString(_k0_trueProcessEnd[i]);
            Int_t ndRecoTrue = _k0_trueParentNDaughters[i];
            std::string labelText;
            if (!procLabel.empty()) {
                labelText = Form("%s, true daughters=%d", procLabel.c_str(), ndRecoTrue);
            } else {
                labelText = Form("true daughters=%d", ndRecoTrue);
            }
            if (annotationGroup && !labelText.empty()) {
                TEveText* procText = new TEveText(labelText.c_str());
                procText->SetMainColor(trueColor);
                procText->SetFontSize(18);
                TEveTrans& procTrans = procText->RefMainTrans();
                procTrans.SetPos(_k0_trueEndPos[i][0] + 5.f,
                                 _k0_trueEndPos[i][1] + 5.f,
                                 _k0_trueEndPos[i][2] + 5.f);
                addElement(annotationGroup, procText);
            }

            if (_k0_trueParentPDG[i] != 0 && _k0_trueParentStartPos[i][0] > -900 && _k0_trueParentEndPos[i][0] > -900) {
                Int_t parentColor = GetParticleColor(_k0_trueParentPDG[i]);
                if (parentColor == kBlack) parentColor = kGray + 1;
                TEveLine* parentLine = new TEveLine(Form("K0 #%d True Parent", i));
                parentLine->SetPoint(0, _k0_trueParentStartPos[i][0], _k0_trueParentStartPos[i][1], _k0_trueParentStartPos[i][2]);
                parentLine->SetPoint(1, _k0_trueParentEndPos[i][0], _k0_trueParentEndPos[i][1], _k0_trueParentEndPos[i][2]);
                parentLine->SetMainColor(parentColor);
                parentLine->SetLineWidth(3);
                parentLine->SetLineStyle(1);
                addElement(truthParentGroup, parentLine);
                anchorPoint(parentLine,
                            _k0_trueParentStartPos[i][0],
                            _k0_trueParentStartPos[i][1],
                            _k0_trueParentStartPos[i][2]);
                anchorPoint(parentLine,
                            _k0_trueParentEndPos[i][0],
                            _k0_trueParentEndPos[i][1],
                            _k0_trueParentEndPos[i][2]);
            }

            for (Int_t d = 0; d < _k0_trueParentNDaughters[i] && d < kMaxTrueDaughters; ++d) {
                Float_t sx = _k0_trueParentDaughterStartPos[i][d*3 + 0];
                Float_t sy = _k0_trueParentDaughterStartPos[i][d*3 + 1];
                Float_t sz = _k0_trueParentDaughterStartPos[i][d*3 + 2];
                Float_t ex = _k0_trueParentDaughterEndPos[i][d*3 + 0];
                Float_t ey = _k0_trueParentDaughterEndPos[i][d*3 + 1];
                Float_t ez = _k0_trueParentDaughterEndPos[i][d*3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Int_t daughterColor = GetParticleColor(_k0_trueParentDaughterPDG[i][d]);
                if (daughterColor == kBlack) daughterColor = kGray + 1;
                TEveLine* daughterLine = new TEveLine(Form("K0 #%d True Daughter %d", i, d));
                daughterLine->SetPoint(0, sx, sy, sz);
                daughterLine->SetPoint(1, ex, ey, ez);
                daughterLine->SetMainColor(daughterColor);
                daughterLine->SetLineWidth(3);
                daughterLine->SetLineStyle(1);
                addElement(truthDaughterGroup, daughterLine);
                anchorPoint(daughterLine, sx, sy, sz);
                anchorPoint(daughterLine, ex, ey, ez);
            }

            for (Int_t s = 0; s < _k0_trueNSiblings[i] && s < kMaxTrueSiblings; ++s) {
                Float_t sx = _k0_trueSiblingStartPos[i][s*3 + 0];
                Float_t sy = _k0_trueSiblingStartPos[i][s*3 + 1];
                Float_t sz = _k0_trueSiblingStartPos[i][s*3 + 2];
                Float_t ex = _k0_trueSiblingEndPos[i][s*3 + 0];
                Float_t ey = _k0_trueSiblingEndPos[i][s*3 + 1];
                Float_t ez = _k0_trueSiblingEndPos[i][s*3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Int_t siblingColor = GetParticleColor(_k0_trueSiblingPDG[i][s]);
                if (siblingColor == kBlack) siblingColor = kGray + 1;
                TEveLine* siblingLine = new TEveLine(Form("K0 #%d True Sibling %d", i, s));
                siblingLine->SetPoint(0, sx, sy, sz);
                siblingLine->SetPoint(1, ex, ey, ez);
                siblingLine->SetMainColor(siblingColor);
                siblingLine->SetLineWidth(2);
                siblingLine->SetLineStyle(1);
                addElement(truthSiblingGroup, siblingLine);
                anchorPoint(siblingLine, sx, sy, sz);
                anchorPoint(siblingLine, ex, ey, ez);
            }
        }

        if (!_k0_hasTrueObject[i]) {
            const bool hasTrueCreationFromParent = (_k0_trueParentEndPos[i][0] > -900 &&
                                                   _k0_trueParentEndPos[i][1] > -900 &&
                                                   _k0_trueParentEndPos[i][2] > -900);
            const bool hasTrueAnnihFromRecoDau = (_k0_trueDecayVtxFromRecoDaughters[i][0] > -900 &&
                                                  _k0_trueDecayVtxFromRecoDaughters[i][1] > -900 &&
                                                  _k0_trueDecayVtxFromRecoDaughters[i][2] > -900);
            TVector3 fallbackLabelPos(
                hasTrueAnnihFromRecoDau ? _k0_trueDecayVtxFromRecoDaughters[i][0] : _k0_trueParentEndPos[i][0],
                hasTrueAnnihFromRecoDau ? _k0_trueDecayVtxFromRecoDaughters[i][1] : _k0_trueParentEndPos[i][1],
                hasTrueAnnihFromRecoDau ? _k0_trueDecayVtxFromRecoDaughters[i][2] : _k0_trueParentEndPos[i][2]);

            if (hasTrueCreationFromParent) {
                TEvePointSet* trueCreation = new TEvePointSet(Form("K0 #%d True Creation Vertex", i));
                trueCreation->SetNextPoint(_k0_trueParentEndPos[i][0], _k0_trueParentEndPos[i][1], _k0_trueParentEndPos[i][2]);
                trueCreation->SetMarkerStyle(29);
                trueCreation->SetMarkerSize(2.3);
                trueCreation->SetMainColor(kBlue + 1);
                addElement(truthGroup, trueCreation);
                anchorPoint(trueCreation, _k0_trueParentEndPos[i][0], _k0_trueParentEndPos[i][1], _k0_trueParentEndPos[i][2]);
            }

            if (hasTrueAnnihFromRecoDau) {
                TEvePointSet* trueAnnih = new TEvePointSet(Form("K0 #%d True Annihilation Vertex", i));
                trueAnnih->SetNextPoint(_k0_trueDecayVtxFromRecoDaughters[i][0], _k0_trueDecayVtxFromRecoDaughters[i][1], _k0_trueDecayVtxFromRecoDaughters[i][2]);
                trueAnnih->SetMarkerStyle(29);
                trueAnnih->SetMarkerSize(2.3);
                trueAnnih->SetMainColor(kRed + 1);
                addElement(truthGroup, trueAnnih);
                anchorPoint(trueAnnih, _k0_trueDecayVtxFromRecoDaughters[i][0], _k0_trueDecayVtxFromRecoDaughters[i][1], _k0_trueDecayVtxFromRecoDaughters[i][2]);
            }

            std::string procLabel = ProcessEnumToString(_k0_trueProcessEnd[i]);
            Int_t ndRecoTrue = _k0_trueParentNDaughters[i];
            std::string labelText;
            if (!procLabel.empty()) {
                labelText = Form("%s, true daughters=%d", procLabel.c_str(), ndRecoTrue);
            } else {
                labelText = Form("true daughters=%d", ndRecoTrue);
            }
            if (annotationGroup && !labelText.empty()) {
                TEveText* procText = new TEveText(labelText.c_str());
                procText->SetMainColor(kRed + 1);
                procText->SetFontSize(18);
                TEveTrans& procTrans = procText->RefMainTrans();
                procTrans.SetPos(fallbackLabelPos.X() + 5.f,
                                 fallbackLabelPos.Y() + 5.f,
                                 fallbackLabelPos.Z() + 5.f);
                addElement(annotationGroup, procText);
            }

            if (_k0_trueParentPDG[i] != 0 && _k0_trueParentStartPos[i][0] > -900 && _k0_trueParentEndPos[i][0] > -900) {
                Int_t parentColor = GetParticleColor(_k0_trueParentPDG[i]);
                if (parentColor == kBlack) parentColor = kGray + 1;
                TEveLine* parentLine = new TEveLine(Form("K0 #%d True Parent", i));
                parentLine->SetPoint(0, _k0_trueParentStartPos[i][0], _k0_trueParentStartPos[i][1], _k0_trueParentStartPos[i][2]);
                parentLine->SetPoint(1, _k0_trueParentEndPos[i][0], _k0_trueParentEndPos[i][1], _k0_trueParentEndPos[i][2]);
                parentLine->SetMainColor(parentColor);
                parentLine->SetLineWidth(3);
                parentLine->SetLineStyle(1);
                addElement(truthParentGroup, parentLine);
                anchorPoint(parentLine, _k0_trueParentStartPos[i][0], _k0_trueParentStartPos[i][1], _k0_trueParentStartPos[i][2]);
                anchorPoint(parentLine, _k0_trueParentEndPos[i][0], _k0_trueParentEndPos[i][1], _k0_trueParentEndPos[i][2]);
            }

            for (Int_t d = 0; d < _k0_trueParentNDaughters[i] && d < kMaxTrueDaughters; ++d) {
                Float_t sx = _k0_trueParentDaughterStartPos[i][d * 3 + 0];
                Float_t sy = _k0_trueParentDaughterStartPos[i][d * 3 + 1];
                Float_t sz = _k0_trueParentDaughterStartPos[i][d * 3 + 2];
                Float_t ex = _k0_trueParentDaughterEndPos[i][d * 3 + 0];
                Float_t ey = _k0_trueParentDaughterEndPos[i][d * 3 + 1];
                Float_t ez = _k0_trueParentDaughterEndPos[i][d * 3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Int_t daughterColor = GetParticleColor(_k0_trueParentDaughterPDG[i][d]);
                if (daughterColor == kBlack) daughterColor = kGray + 1;
                TEveLine* daughterLine = new TEveLine(Form("K0 #%d True Vertex Daughter %d", i, d));
                daughterLine->SetPoint(0, sx, sy, sz);
                daughterLine->SetPoint(1, ex, ey, ez);
                daughterLine->SetMainColor(daughterColor);
                daughterLine->SetLineWidth(4);
                daughterLine->SetLineStyle(1);
                addElement(truthDaughterGroup, daughterLine);
                anchorPoint(daughterLine, sx, sy, sz);
                anchorPoint(daughterLine, ex, ey, ez);

                TEvePointSet* dauStart = new TEvePointSet(Form("K0 #%d True Vertex Daughter %d Start", i, d));
                dauStart->SetNextPoint(sx, sy, sz);
                dauStart->SetMarkerStyle(29);
                dauStart->SetMarkerSize(2.2);
                dauStart->SetMainColor(daughterColor);
                addElement(truthDaughterGroup, dauStart);
                anchorPoint(dauStart, sx, sy, sz);

                TEvePointSet* dauEnd = new TEvePointSet(Form("K0 #%d True Vertex Daughter %d End", i, d));
                dauEnd->SetNextPoint(ex, ey, ez);
                dauEnd->SetMarkerStyle(20);
                dauEnd->SetMarkerSize(1.8);
                dauEnd->SetMainColor(daughterColor);
                addElement(truthDaughterGroup, dauEnd);
                anchorPoint(dauEnd, ex, ey, ez);
            }
        }

        // Draw creation vertex fitted lines (beam + second particle)
        Float_t fitLength = _k0_fitLineLength[i];
        if (fitLength > 0 && _k0_creationVtx_fitLineBeamStart[i][0] > -900 &&
            _k0_creationVtx_closestPtBeam[i][0] > -900) {
            // Get beam particle ID (parent)
            Int_t beamID = _k0_parentID[i];
            Int_t beamPDG = -1;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == beamID) {
                    beamPDG = _particle_PDG[p];
                    break;
                }
            }
            Int_t beamColor = (beamPDG > 0) ? GetParticleColor(beamPDG) : kBlack;

            // Second particle color (from PDG)
            Int_t secondID = _k0_secondParticleID[i];
            Int_t secondPDG = -1;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == secondID) {
                    secondPDG = _particle_PDG[p];
                    break;
                }
            }
            Int_t secondColor = (secondPDG > 0) ? GetParticleColor(secondPDG) : kGray;

            // Get beam particle actual endPos and verify it has hits (3D)
            Float_t beamPosX = -999, beamPosY = -999, beamPosZ = -999;
            Float_t beamStartX = -999, beamStartY = -999, beamStartZ = -999;
            Int_t beamHits = 0;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == beamID) {
                    beamPosX = _particle_endPos[p][0];  // Beam uses END position
                    beamPosY = _particle_endPos[p][1];
                    beamPosZ = _particle_endPos[p][2];
                    beamStartX = _particle_startPos[p][0];
                    beamStartY = _particle_startPos[p][1];
                    beamStartZ = _particle_startPos[p][2];
                    beamHits = _particle_nHits[p];
                    break;
                }
            }

            // Only draw if beam particle position was found and particle has hits
            if (beamPosX < -900 || beamHits == 0) continue;

            // Beam fitted line: draw bidirectionally with span >= reconstructed particle length.
            TVector3 beamPos(beamPosX, beamPosY, beamPosZ);
            TVector3 closestBeam(_k0_creationVtx_closestPtBeam[i][0],
                                _k0_creationVtx_closestPtBeam[i][1],
                                _k0_creationVtx_closestPtBeam[i][2]);
            TVector3 beamDir = closestBeam - beamPos;
            if (beamDir.Mag() > 0) beamDir = beamDir.Unit();
            Float_t beamRecoLength = 0.f;
            if (beamStartX > -900.f && beamStartY > -900.f && beamStartZ > -900.f) {
                beamRecoLength = TVector3(beamPosX - beamStartX, beamPosY - beamStartY, beamPosZ - beamStartZ).Mag();
            }
            const double beamSpan = std::max<double>(100.0, std::max<double>(beamRecoLength, fitLength));
            TVector3 beamLineStart = beamPos - beamSpan * beamDir;
            TVector3 beamLineEnd = beamPos + beamSpan * beamDir;

            TEveLine* beamLine = new TEveLine(Form("K0 #%d Creation Beam Fit", i));
            beamLine->SetPoint(0, beamLineStart.X(),
                                  beamLineStart.Y(),
                                  beamLineStart.Z());
            beamLine->SetPoint(1, beamLineEnd.X(),
                                  beamLineEnd.Y(),
                                  beamLineEnd.Z());
            beamLine->SetMainColor(beamColor);
            beamLine->SetLineStyle(2); // Dashed (fitted line uses endPos/endDir for beam)
            beamLine->SetLineWidth(2);
            addElement(creationFitGroup, beamLine);
            anchorPoint(beamLine, beamLineStart.X(), beamLineStart.Y(), beamLineStart.Z());
            anchorPoint(beamLine, beamLineEnd.X(), beamLineEnd.Y(), beamLineEnd.Z());

            // Second particle fitted line (uses startPos/startDir)
            if (_k0_creationVtx_fitLineSecondStart[i][0] > -900 &&
                _k0_creationVtx_closestPtSecond[i][0] > -900) {
                // Get second particle actual startPos and verify it has hits (3D)
                Float_t secondPosX = -999, secondPosY = -999, secondPosZ = -999;
                Float_t secondEndX = -999, secondEndY = -999, secondEndZ = -999;
                Int_t secondHits = 0;
                for (Int_t p = 0; p < _nParticles; p++) {
                    if (_particle_uniqueID[p] == secondID) {
                        secondPosX = _particle_startPos[p][0];  // Second uses START position
                        secondPosY = _particle_startPos[p][1];
                        secondPosZ = _particle_startPos[p][2];
                        secondEndX = _particle_endPos[p][0];
                        secondEndY = _particle_endPos[p][1];
                        secondEndZ = _particle_endPos[p][2];
                        secondHits = _particle_nHits[p];
                        break;
                    }
                }

                // Only draw if second particle position was found and particle has hits
                if (secondPosX < -900 || secondHits == 0) continue;

                // Draw bidirectional fitted line with span >= reconstructed particle length.
                TVector3 secondPos(secondPosX, secondPosY, secondPosZ);
                TVector3 closestSecond(_k0_creationVtx_closestPtSecond[i][0],
                                      _k0_creationVtx_closestPtSecond[i][1],
                                      _k0_creationVtx_closestPtSecond[i][2]);
                TVector3 secondDir = closestSecond - secondPos;
                if (secondDir.Mag() > 0) secondDir = secondDir.Unit();
                Float_t secondRecoLength = 0.f;
                if (secondEndX > -900.f && secondEndY > -900.f && secondEndZ > -900.f) {
                    secondRecoLength = TVector3(secondEndX - secondPosX, secondEndY - secondPosY, secondEndZ - secondPosZ).Mag();
                }
                const double secondSpan = std::max<double>(100.0, std::max<double>(secondRecoLength, fitLength));
                TVector3 secondLineStart = secondPos - secondSpan * secondDir;
                TVector3 secondLineEnd = secondPos + secondSpan * secondDir;

                TEveLine* secondLine = new TEveLine(Form("K0 #%d Creation Second Fit", i));
                secondLine->SetPoint(0, secondLineStart.X(),
                                        secondLineStart.Y(),
                                        secondLineStart.Z());
                secondLine->SetPoint(1, secondLineEnd.X(),
                                        secondLineEnd.Y(),
                                        secondLineEnd.Z());
                secondLine->SetMainColor(secondColor);
                secondLine->SetLineStyle(2); // Dashed (fitted line uses startPos/startDir)
                secondLine->SetLineWidth(2);
                addElement(creationFitGroup, secondLine);
                anchorPoint(secondLine, secondLineStart.X(), secondLineStart.Y(), secondLineStart.Z());
                anchorPoint(secondLine, secondLineEnd.X(), secondLineEnd.Y(), secondLineEnd.Z());
            }

            // Draw closest points on creation vertex fitted lines
            if (_k0_creationVtx_closestPtBeam[i][0] > -900) {
                TEvePointSet* cpBeam = new TEvePointSet(Form("K0 #%d Creation Beam Closest", i));
                cpBeam->SetNextPoint(_k0_creationVtx_closestPtBeam[i][0],
                                     _k0_creationVtx_closestPtBeam[i][1],
                                     _k0_creationVtx_closestPtBeam[i][2]);
                cpBeam->SetMarkerStyle(24); // Open circle
                cpBeam->SetMarkerSize(2.0);
                cpBeam->SetMainColor(kOrange);
                addElement(creationFitGroup, cpBeam);
                anchorPoint(cpBeam,
                            _k0_creationVtx_closestPtBeam[i][0],
                            _k0_creationVtx_closestPtBeam[i][1],
                            _k0_creationVtx_closestPtBeam[i][2]);
            }

            if (_k0_creationVtx_closestPtSecond[i][0] > -900) {
                TEvePointSet* cpSecond = new TEvePointSet(Form("K0 #%d Creation Second Closest", i));
                cpSecond->SetNextPoint(_k0_creationVtx_closestPtSecond[i][0],
                                       _k0_creationVtx_closestPtSecond[i][1],
                                       _k0_creationVtx_closestPtSecond[i][2]);
                cpSecond->SetMarkerStyle(24); // Open circle
                cpSecond->SetMarkerSize(2.0);
                cpSecond->SetMainColor(kOrange);
                addElement(creationFitGroup, cpSecond);
                anchorPoint(cpSecond,
                            _k0_creationVtx_closestPtSecond[i][0],
                            _k0_creationVtx_closestPtSecond[i][1],
                            _k0_creationVtx_closestPtSecond[i][2]);
            }

            // Draw white dotted line connecting the two closest points
            if (_k0_creationVtx_closestPtBeam[i][0] > -900 && _k0_creationVtx_closestPtSecond[i][0] > -900) {
                TEveLine* connectLine = new TEveLine(Form("K0 #%d Creation Vtx Connect", i));
                connectLine->SetPoint(0, _k0_creationVtx_closestPtBeam[i][0],
                                         _k0_creationVtx_closestPtBeam[i][1],
                                         _k0_creationVtx_closestPtBeam[i][2]);
                connectLine->SetPoint(1, _k0_creationVtx_closestPtSecond[i][0],
                                         _k0_creationVtx_closestPtSecond[i][1],
                                         _k0_creationVtx_closestPtSecond[i][2]);
                connectLine->SetMainColor(kWhite);
                connectLine->SetLineStyle(3); // Dotted
                connectLine->SetLineWidth(2);
                addElement(creationFitGroup, connectLine);
            }
        }

        // Draw annihilation vertex fitted lines (daughter particles)
        if (fitLength > 0) {
            TEveElementList* annihilationFitCandidateGroup =
                new TEveElementList(Form("K0 #%d", i));
            addElement(annihilationFitGroup, annihilationFitCandidateGroup);

            // Get daughter colors by matching UniqueID
            Int_t d1Color = kCyan;
            Int_t d2Color = kSpring;
            Int_t d1ID = _k0_daughter1ID[i];
            Int_t d2ID = _k0_daughter2ID[i];

            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == d1ID) {
                    d1Color = GetParticleColor(_particle_PDG[p]);
                }
                if (_particle_uniqueID[p] == d2ID) {
                    d2Color = GetParticleColor(_particle_PDG[p]);
                }
            }

            // Daughter 1: Pandora direction (solid)
            if (_k0_annVtx_fitLine1Start[i][0] > -900 && _k0_annVtx_pandoraLine1Dir[i][0] > -900) {
                TVector3 s(_k0_annVtx_fitLine1Start[i][0], _k0_annVtx_fitLine1Start[i][1], _k0_annVtx_fitLine1Start[i][2]);
                TVector3 d(_k0_annVtx_pandoraLine1Dir[i][0], _k0_annVtx_pandoraLine1Dir[i][1], _k0_annVtx_pandoraLine1Dir[i][2]);
                if (d.Mag() > 0) {
                    d = d.Unit();
                    TVector3 lineStart = s;
                    TVector3 lineEnd = s + 100.0 * d;
                    if (_k0_annVtx_pandoraClosestPt1[i][0] > -900) {
                        TVector3 cp(_k0_annVtx_pandoraClosestPt1[i][0], _k0_annVtx_pandoraClosestPt1[i][1], _k0_annVtx_pandoraClosestPt1[i][2]);
                        const double t = (cp - s).Dot(d);
                        if (t >= 0.0) {
                            lineStart = s;
                            lineEnd = cp + 20.0 * d;
                        } else {
                            lineStart = cp;
                            lineEnd = s + 100.0 * d;
                        }
                    }
                    TEveLine* line = new TEveLine(Form("K0 #%d D1 Pandora Dir", i));
                    line->SetPoint(0, lineStart.X(), lineStart.Y(), lineStart.Z());
                    line->SetPoint(1, lineEnd.X(), lineEnd.Y(), lineEnd.Z());
                    line->SetMainColor(d1Color);
                    line->SetLineWidth(1);
                    line->SetLineStyle(1);
                    addElement(momentumArrowGroup, line);
                }
            }

            // Daughter 1: extrapolated-fit direction (dashed)
            if (_k0_annVtx_fitLine1Start[i][0] > -900 && _k0_annVtx_fitLine1Dir[i][0] > -900) {
                TVector3 s(_k0_annVtx_fitLine1Start[i][0], _k0_annVtx_fitLine1Start[i][1], _k0_annVtx_fitLine1Start[i][2]);
                TVector3 d(_k0_annVtx_fitLine1Dir[i][0], _k0_annVtx_fitLine1Dir[i][1], _k0_annVtx_fitLine1Dir[i][2]);
                if (d.Mag() > 0) {
                    d = d.Unit();
                    TVector3 base = s;
                    double d1RecoLength = 0.0;
                    for (Int_t p = 0; p < _nParticles; ++p) {
                        if (_particle_uniqueID[p] != d1ID) continue;
                        const Float_t sx = _particle_startPos[p][0];
                        const Float_t sy = _particle_startPos[p][1];
                        const Float_t sz = _particle_startPos[p][2];
                        const Float_t ex = _particle_endPos[p][0];
                        const Float_t ey = _particle_endPos[p][1];
                        const Float_t ez = _particle_endPos[p][2];
                        if (sx > -900.f && sy > -900.f && sz > -900.f && ex > -900.f && ey > -900.f && ez > -900.f) {
                            d1RecoLength = TVector3(ex - sx, ey - sy, ez - sz).Mag();
                        }
                        break;
                    }
                    const double span = std::max<double>(100.0, std::max<double>(d1RecoLength, fitLength));
                    double tMin = -span;
                    double tMax = span;
                    if (_k0_annVtx_closestPt1[i][0] > -900) {
                        TVector3 cp(_k0_annVtx_closestPt1[i][0], _k0_annVtx_closestPt1[i][1], _k0_annVtx_closestPt1[i][2]);
                        const TVector3 anchorToClosest = cp - s;
                        if (anchorToClosest.Mag2() > 1e-8) {
                            d = anchorToClosest.Unit();
                        }
                        base = cp; // Guarantee rendered fit line passes through closest point.
                        const double tAnchor = (s - base).Dot(d);
                        tMin = std::min(tMin, tAnchor - span);
                        tMax = std::max(tMax, tAnchor + span);
                    }
                    TVector3 lineStart = base + tMin * d;
                    TVector3 lineEnd = base + tMax * d;
                    TEveLine* line = new TEveLine(Form("K0 #%d D1 Fit Dir", i));
                    line->SetPoint(0, lineStart.X(), lineStart.Y(), lineStart.Z());
                    line->SetPoint(1, lineEnd.X(), lineEnd.Y(), lineEnd.Z());
                    line->SetMainColor(d1Color);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    addElement(momentumArrowGroup, line);
                }
            }

            // Daughter 2: Pandora direction (solid)
            if (_k0_annVtx_fitLine2Start[i][0] > -900 && _k0_annVtx_pandoraLine2Dir[i][0] > -900) {
                TVector3 s(_k0_annVtx_fitLine2Start[i][0], _k0_annVtx_fitLine2Start[i][1], _k0_annVtx_fitLine2Start[i][2]);
                TVector3 d(_k0_annVtx_pandoraLine2Dir[i][0], _k0_annVtx_pandoraLine2Dir[i][1], _k0_annVtx_pandoraLine2Dir[i][2]);
                if (d.Mag() > 0) {
                    d = d.Unit();
                    TVector3 lineStart = s;
                    TVector3 lineEnd = s + 100.0 * d;
                    if (_k0_annVtx_pandoraClosestPt2[i][0] > -900) {
                        TVector3 cp(_k0_annVtx_pandoraClosestPt2[i][0], _k0_annVtx_pandoraClosestPt2[i][1], _k0_annVtx_pandoraClosestPt2[i][2]);
                        const double t = (cp - s).Dot(d);
                        if (t >= 0.0) {
                            lineStart = s;
                            lineEnd = cp + 20.0 * d;
                        } else {
                            lineStart = cp;
                            lineEnd = s + 100.0 * d;
                        }
                    }
                    TEveLine* line = new TEveLine(Form("K0 #%d D2 Pandora Dir", i));
                    line->SetPoint(0, lineStart.X(), lineStart.Y(), lineStart.Z());
                    line->SetPoint(1, lineEnd.X(), lineEnd.Y(), lineEnd.Z());
                    line->SetMainColor(d2Color);
                    line->SetLineWidth(1);
                    line->SetLineStyle(1);
                    addElement(momentumArrowGroup, line);
                }
            }

            // Daughter 2: extrapolated-fit direction (dashed)
            if (_k0_annVtx_fitLine2Start[i][0] > -900 && _k0_annVtx_fitLine2Dir[i][0] > -900) {
                TVector3 s(_k0_annVtx_fitLine2Start[i][0], _k0_annVtx_fitLine2Start[i][1], _k0_annVtx_fitLine2Start[i][2]);
                TVector3 d(_k0_annVtx_fitLine2Dir[i][0], _k0_annVtx_fitLine2Dir[i][1], _k0_annVtx_fitLine2Dir[i][2]);
                if (d.Mag() > 0) {
                    d = d.Unit();
                    TVector3 base = s;
                    double d2RecoLength = 0.0;
                    for (Int_t p = 0; p < _nParticles; ++p) {
                        if (_particle_uniqueID[p] != d2ID) continue;
                        const Float_t sx = _particle_startPos[p][0];
                        const Float_t sy = _particle_startPos[p][1];
                        const Float_t sz = _particle_startPos[p][2];
                        const Float_t ex = _particle_endPos[p][0];
                        const Float_t ey = _particle_endPos[p][1];
                        const Float_t ez = _particle_endPos[p][2];
                        if (sx > -900.f && sy > -900.f && sz > -900.f && ex > -900.f && ey > -900.f && ez > -900.f) {
                            d2RecoLength = TVector3(ex - sx, ey - sy, ez - sz).Mag();
                        }
                        break;
                    }
                    const double span = std::max<double>(100.0, std::max<double>(d2RecoLength, fitLength));
                    double tMin = -span;
                    double tMax = span;
                    if (_k0_annVtx_closestPt2[i][0] > -900) {
                        TVector3 cp(_k0_annVtx_closestPt2[i][0], _k0_annVtx_closestPt2[i][1], _k0_annVtx_closestPt2[i][2]);
                        const TVector3 anchorToClosest = cp - s;
                        if (anchorToClosest.Mag2() > 1e-8) {
                            d = anchorToClosest.Unit();
                        }
                        base = cp; // Guarantee rendered fit line passes through closest point.
                        const double tAnchor = (s - base).Dot(d);
                        tMin = std::min(tMin, tAnchor - span);
                        tMax = std::max(tMax, tAnchor + span);
                    }
                    TVector3 lineStart = base + tMin * d;
                    TVector3 lineEnd = base + tMax * d;
                    TEveLine* line = new TEveLine(Form("K0 #%d D2 Fit Dir", i));
                    line->SetPoint(0, lineStart.X(), lineStart.Y(), lineStart.Z());
                    line->SetPoint(1, lineEnd.X(), lineEnd.Y(), lineEnd.Z());
                    line->SetMainColor(d2Color);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    addElement(momentumArrowGroup, line);
                }
            }

            // Closest points (points of minimum distance)
            if (_k0_annVtx_pandoraClosestPt1[i][0] > -900) {
                TEvePointSet* pt1Pand = new TEvePointSet(Form("K0 #%d Pandora Closest Pt D1", i));
                pt1Pand->SetNextPoint(_k0_annVtx_pandoraClosestPt1[i][0],
                                      _k0_annVtx_pandoraClosestPt1[i][1],
                                      _k0_annVtx_pandoraClosestPt1[i][2]);
                pt1Pand->SetMarkerStyle(25); // Open square
                pt1Pand->SetMarkerSize(1.8);
                pt1Pand->SetMainColor(kAzure + 2);
                addElement(annihilationFitCandidateGroup, pt1Pand);
            }

            if (_k0_annVtx_pandoraClosestPt2[i][0] > -900) {
                TEvePointSet* pt2Pand = new TEvePointSet(Form("K0 #%d Pandora Closest Pt D2", i));
                pt2Pand->SetNextPoint(_k0_annVtx_pandoraClosestPt2[i][0],
                                      _k0_annVtx_pandoraClosestPt2[i][1],
                                      _k0_annVtx_pandoraClosestPt2[i][2]);
                pt2Pand->SetMarkerStyle(25); // Open square
                pt2Pand->SetMarkerSize(1.8);
                pt2Pand->SetMainColor(kAzure + 2);
                addElement(annihilationFitCandidateGroup, pt2Pand);
            }

            if (_k0_annVtx_pandoraClosestPt1[i][0] > -900 && _k0_annVtx_pandoraClosestPt2[i][0] > -900) {
                TEveLine* connectPand = new TEveLine(Form("K0 #%d Ann Vtx Connect Pandora", i));
                connectPand->SetPoint(0, _k0_annVtx_pandoraClosestPt1[i][0],
                                      _k0_annVtx_pandoraClosestPt1[i][1],
                                      _k0_annVtx_pandoraClosestPt1[i][2]);
                connectPand->SetPoint(1, _k0_annVtx_pandoraClosestPt2[i][0],
                                      _k0_annVtx_pandoraClosestPt2[i][1],
                                      _k0_annVtx_pandoraClosestPt2[i][2]);
                connectPand->SetMainColor(kAzure + 2);
                connectPand->SetLineStyle(1);
                connectPand->SetLineWidth(2);
                addElement(annihilationFitCandidateGroup, connectPand);
            }

            if (_k0_annVtx_closestPt1[i][0] > -900) {
                TEvePointSet* pt1 = new TEvePointSet(Form("K0 #%d Closest Pt D1", i));
                pt1->SetNextPoint(_k0_annVtx_closestPt1[i][0],
                                 _k0_annVtx_closestPt1[i][1],
                                 _k0_annVtx_closestPt1[i][2]);
                pt1->SetMarkerStyle(24); // Open circle
                pt1->SetMarkerSize(2.0);
                pt1->SetMainColor(kOrange);
                addElement(annihilationFitCandidateGroup, pt1);
            }

            if (_k0_annVtx_closestPt2[i][0] > -900) {
                TEvePointSet* pt2 = new TEvePointSet(Form("K0 #%d Closest Pt D2", i));
                pt2->SetNextPoint(_k0_annVtx_closestPt2[i][0],
                                 _k0_annVtx_closestPt2[i][1],
                                 _k0_annVtx_closestPt2[i][2]);
                pt2->SetMarkerStyle(24); // Open circle
                pt2->SetMarkerSize(2.0);
                pt2->SetMainColor(kOrange);
                addElement(annihilationFitCandidateGroup, pt2);
            }

            // Draw white dotted line connecting the two closest points (annihilation vertex)
            if (_k0_annVtx_closestPt1[i][0] > -900 && _k0_annVtx_closestPt2[i][0] > -900) {
                TEveLine* connectLine = new TEveLine(Form("K0 #%d Ann Vtx Connect", i));
                connectLine->SetPoint(0, _k0_annVtx_closestPt1[i][0],
                                         _k0_annVtx_closestPt1[i][1],
                                         _k0_annVtx_closestPt1[i][2]);
                connectLine->SetPoint(1, _k0_annVtx_closestPt2[i][0],
                                         _k0_annVtx_closestPt2[i][1],
                                         _k0_annVtx_closestPt2[i][2]);
                connectLine->SetMainColor(kWhite);
                connectLine->SetLineStyle(3); // Dotted
                connectLine->SetLineWidth(2);
                addElement(annihilationFitCandidateGroup, connectLine);
            }
        }
    }

    if (_nK0Candidates > 0) {
        std::cout << "Drew " << _nK0Candidates << " K0 candidates in 3D" << std::endl;
    }

    // Draw all true particle trajectories
    for (Int_t i = 0; i < _nAllTrueParticles && i < kMaxAllTrueParticles; ++i) {
        if (_allTrueParticle_startPos[i][0] <= -900 || _allTrueParticle_endPos[i][0] <= -900) continue;

        TEveElementList* trueParticleGroup =
            new TEveElementList(Form("True UID=%d PDG=%d", i, _allTrueParticle_PDG[i]));
        addElement(trueParticlesGroup, trueParticleGroup);

        Int_t trueColor = GetParticleColor(_allTrueParticle_PDG[i]);
        if (trueColor == kBlack) trueColor = kGray + 1;

        TEveLine* trueLine = new TEveLine(Form("True Particle #%d", i));
        trueLine->SetPoint(0, _allTrueParticle_startPos[i][0], _allTrueParticle_startPos[i][1], _allTrueParticle_startPos[i][2]);
        trueLine->SetPoint(1, _allTrueParticle_endPos[i][0], _allTrueParticle_endPos[i][1], _allTrueParticle_endPos[i][2]);
        trueLine->SetMainColor(trueColor);
        trueLine->SetLineWidth(2);
        trueLine->SetLineStyle(1);
        addElement(trueParticleGroup, trueLine);

        TEvePointSet* trueStart = new TEvePointSet(Form("True Particle #%d Start", i));
        trueStart->SetNextPoint(_allTrueParticle_startPos[i][0], _allTrueParticle_startPos[i][1], _allTrueParticle_startPos[i][2]);
        trueStart->SetMarkerStyle(29);
        trueStart->SetMarkerSize(1.8);
        trueStart->SetMainColor(trueColor);
        addElement(trueParticleGroup, trueStart);

        TEvePointSet* trueEnd = new TEvePointSet(Form("True Particle #%d End", i));
        trueEnd->SetNextPoint(_allTrueParticle_endPos[i][0], _allTrueParticle_endPos[i][1], _allTrueParticle_endPos[i][2]);
        trueEnd->SetMarkerStyle(29);
        trueEnd->SetMarkerSize(1.8);
        trueEnd->SetMainColor(trueColor);
        addElement(trueParticleGroup, trueEnd);

        std::string procLabel = ProcessEnumToString(_allTrueParticle_processEnd[i]);
        std::string labelText;
        if (!procLabel.empty()) {
            labelText = Form("%s, PDG=%d", procLabel.c_str(), _allTrueParticle_PDG[i]);
        } else {
            labelText = Form("PDG=%d", _allTrueParticle_PDG[i]);
        }
        if (!labelText.empty()) {
            TEveText* procText = new TEveText(labelText.c_str());
            procText->SetMainColor(trueColor);
            procText->SetFontSize(14);
            TEveTrans& procTrans = procText->RefMainTrans();
            procTrans.SetPos(_allTrueParticle_endPos[i][0], _allTrueParticle_endPos[i][1], _allTrueParticle_endPos[i][2]);
            addElement(trueParticleGroup, procText);
        }
    }
}

//********************************************************************
void neutralKaonEventDisplay::DrawAnalysisContent2D(TEveProjectionManager* manager, const std::string& projection_type) {
//********************************************************************
    // K0 vertices will be automatically projected from 3D via ImportElements
    // No additional 2D-specific drawing needed (TEve projections disabled)
    (void)manager;
    (void)projection_type;
}

//********************************************************************
void neutralKaonEventDisplay::DrawAnalysisContentCanvas2D(TCanvas* canvas, const std::string& projection_type) {
//********************************************************************
    // Draw K0 vertices (creation and annihilation) as circles on TCanvas
    canvas->cd();

    const Bool_t showCreationVertices = kFALSE;
    const Bool_t showAnnihilationVertices = IsGroupVisible("Neutral Particles");
    const Bool_t showStandaloneTruth = IsGroupVisible("True Particles");
    const Bool_t showMomentumArrows = IsGroupVisible("Neutral Particles");

    for (Int_t i = 0; i < _nK0Candidates && i < kMaxK0 && false; i++) {
        // Get vertex positions
        Float_t creationX = _k0_creationVtxPos[i][0];
        Float_t creationY = _k0_creationVtxPos[i][1];
        Float_t creationZ = _k0_creationVtxPos[i][2];

        Float_t annihilationX = _k0_annihilationVtxPos[i][0];
        Float_t annihilationY = _k0_annihilationVtxPos[i][1];
        Float_t annihilationZ = _k0_annihilationVtxPos[i][2];
        Float_t annihilationFitX = _k0_annihilationVtxFitPos[i][0];
        Float_t annihilationFitY = _k0_annihilationVtxFitPos[i][1];
        Float_t annihilationFitZ = _k0_annihilationVtxFitPos[i][2];

        // Get vertex radii and degeneracies
        Float_t creationRadius = _k0_creationVtxRadius[i];
        Float_t annihilationRadius = _k0_annihilationVtxRadius[i];
        Int_t creationDeg = _k0_creationVtxDeg[i];
        Int_t annihilationDeg = _k0_annihilationVtxDeg[i];

        // Get parent color for creation vertex
        Int_t parentColor = kBlue; // Default
        if (_k0_parentID[i] >= 0 && _k0_parentID[i] < kMaxParticles) {
            Int_t parentPDG = _particle_PDG[_k0_parentID[i]];
            parentColor = GetParticleColor(parentPDG);
        }

        // Project to 2D coordinates
        Float_t creationCoord1, creationCoord2;
        Float_t annCoord1, annCoord2;

        if (projection_type == "XY") {
            creationCoord1 = creationX; creationCoord2 = creationY;
            annCoord1 = annihilationX; annCoord2 = annihilationY;
        } else if (projection_type == "XZ") {
            creationCoord1 = creationX; creationCoord2 = creationZ;
            annCoord1 = annihilationX; annCoord2 = annihilationZ;
        } else if (projection_type == "YZ") {
            creationCoord1 = creationY; creationCoord2 = creationZ;
            annCoord1 = annihilationY; annCoord2 = annihilationZ;
        } else {
            continue;
        }

        Float_t annFitCoord1 = -999.f;
        Float_t annFitCoord2 = -999.f;
        if (projection_type == "XY") {
            annFitCoord1 = annihilationFitX; annFitCoord2 = annihilationFitY;
        } else if (projection_type == "XZ") {
            annFitCoord1 = annihilationFitX; annFitCoord2 = annihilationFitZ;
        } else if (projection_type == "YZ") {
            annFitCoord1 = annihilationFitY; annFitCoord2 = annihilationFitZ;
        }

        // Draw creation vertex circle (if valid position)
        if (showCreationVertices && creationX > -900 && creationRadius > 0) {
            TEllipse* creationCircle = new TEllipse(creationCoord1, creationCoord2, creationRadius, creationRadius);
            creationCircle->SetFillStyle(0); // Hollow
            creationCircle->SetLineColor(parentColor);
            creationCircle->SetLineWidth(2);
            creationCircle->Draw("SAME");

            // Creation vertex position marker (star - same as 3D)
            TMarker* creationVtxMarker = new TMarker(creationCoord1, creationCoord2, 29); // Star
            creationVtxMarker->SetMarkerColor(parentColor);
            creationVtxMarker->SetMarkerSize(2.5);
            creationVtxMarker->Draw("SAME");

            // Add degeneracy label
            if (creationDeg > 0) {
                TLatex* degLabel = new TLatex(creationCoord1 + creationRadius, creationCoord2 + creationRadius,
                                             Form("#scale[0.6]{deg=%d}", creationDeg));
                degLabel->SetTextColor(parentColor);
                degLabel->SetTextSize(0.025);
                degLabel->Draw("SAME");
            }
        }

        // Draw annihilation vertex circle (if valid position)
        if (showAnnihilationVertices && annihilationX > -900 && annihilationRadius > 0) {
            TEllipse* annCircle = new TEllipse(annCoord1, annCoord2, annihilationRadius, annihilationRadius);
            annCircle->SetFillStyle(0); // Hollow
            annCircle->SetLineColor(kRed);
            annCircle->SetLineWidth(2);
            annCircle->Draw("SAME");

            // Annihilation vertex position marker (star - same as 3D)
            TMarker* annVtxMarker = new TMarker(annCoord1, annCoord2, 29); // Star
            annVtxMarker->SetMarkerColor(kRed);
            annVtxMarker->SetMarkerSize(2.5);
            annVtxMarker->Draw("SAME");

            // Fit-based annihilation vertex marker (different marker/color from Pandora)
            if (annFitCoord1 > -900 && annFitCoord2 > -900) {
                TMarker* annFitVtxMarker = new TMarker(annFitCoord1, annFitCoord2, 20); // Filled circle
                annFitVtxMarker->SetMarkerColor(kOrange + 1);
                annFitVtxMarker->SetMarkerSize(1.8);
                annFitVtxMarker->Draw("SAME");
            }

            // Add degeneracy label
            if (annihilationDeg > 0) {
                TLatex* degLabel = new TLatex(annCoord1 + annihilationRadius, annCoord2 + annihilationRadius,
                                             Form("#scale[0.6]{deg=%d}", annihilationDeg));
                degLabel->SetTextColor(kRed);
                degLabel->SetTextSize(0.025);
                degLabel->Draw("SAME");
            }
        }

        // Draw reconstructed neutral trajectory from creation to annihilation.
        if (_k0_creationVtxPos[i][0] > -900 && _k0_annihilationVtxPos[i][0] > -900) {
            Float_t startX = _k0_creationVtxPos[i][0];
            Float_t startY = _k0_creationVtxPos[i][1];
            Float_t startZ = _k0_creationVtxPos[i][2];
            Float_t endX = _k0_annihilationVtxPos[i][0];
            Float_t endY = _k0_annihilationVtxPos[i][1];
            Float_t endZ = _k0_annihilationVtxPos[i][2];

            Float_t x1, y1, x2, y2;
            if (projection_type == "XY") {
                x1 = startX; y1 = startY;
                x2 = endX; y2 = endY;
            } else if (projection_type == "XZ") {
                x1 = startX; y1 = startZ;
                x2 = endX; y2 = endZ;
            } else if (projection_type == "YZ") {
                x1 = startY; y1 = startZ;
                x2 = endY; y2 = endZ;
            } else {
                continue;
            }

            Int_t recoColor = kGray + 1;
            if (_k0_hasTrueObject[i] && _k0_truePDG[i] != -999) {
                recoColor = GetParticleColor(_k0_truePDG[i]);
                if (recoColor == kBlack) recoColor = kGray + 1;
            }

            TLine* k0Line = new TLine(x1, y1, x2, y2);
            k0Line->SetLineColor(recoColor);
            k0Line->SetLineStyle(2); // Dashed
            k0Line->SetLineWidth(3);
            k0Line->Draw("SAME");

        }

        // Draw K0 true trajectory (if TrueObject exists - solid line, K0 color)
        if (_k0_hasTrueObject[i] && _k0_trueStartPos[i][0] > -900 && _k0_trueEndPos[i][0] > -900) {
            Int_t trueColor = GetParticleColor(_k0_truePDG[i]);
            if (trueColor == kBlack) trueColor = kGray + 1;

            Float_t trueStartX = _k0_trueStartPos[i][0];
            Float_t trueStartY = _k0_trueStartPos[i][1];
            Float_t trueStartZ = _k0_trueStartPos[i][2];
            Float_t trueEndX = _k0_trueEndPos[i][0];
            Float_t trueEndY = _k0_trueEndPos[i][1];
            Float_t trueEndZ = _k0_trueEndPos[i][2];

            Float_t tx1, ty1, tx2, ty2;
            if (projection_type == "XY") {
                tx1 = trueStartX; ty1 = trueStartY;
                tx2 = trueEndX; ty2 = trueEndY;
            } else if (projection_type == "XZ") {
                tx1 = trueStartX; ty1 = trueStartZ;
                tx2 = trueEndX; ty2 = trueEndZ;
            } else if (projection_type == "YZ") {
                tx1 = trueStartY; ty1 = trueStartZ;
                tx2 = trueEndY; ty2 = trueEndZ;
            } else {
                continue;
            }

            TLine* trueK0Line = new TLine(tx1, ty1, tx2, ty2);
            trueK0Line->SetLineColor(trueColor);
            trueK0Line->SetLineStyle(1); // Solid for truth
            trueK0Line->SetLineWidth(3);
            trueK0Line->Draw("SAME");

            // Add markers at true start and end positions
            TMarker* trueStartMarker = new TMarker(tx1, ty1, 29); // Star at start
            trueStartMarker->SetMarkerColor(trueColor);
            trueStartMarker->SetMarkerSize(2.0);
            trueStartMarker->Draw("SAME");

            TMarker* trueEndMarker = new TMarker(tx2, ty2, 29); // Star at end
            trueEndMarker->SetMarkerColor(trueColor);
            trueEndMarker->SetMarkerSize(2.0);
            trueEndMarker->Draw("SAME");

            std::string procLabel = ProcessEnumToString(_k0_trueProcessEnd[i]);
            Int_t ndRecoTrue = _k0_trueParentNDaughters[i];
            std::string labelText;
            if (!procLabel.empty()) {
                labelText = Form("#scale[0.6]{%s, true daughters=%d}", procLabel.c_str(), ndRecoTrue);
            } else {
                labelText = Form("#scale[0.6]{true daughters=%d}", ndRecoTrue);
            }
            if (!labelText.empty()) {
                TLatex* procLatex = new TLatex(tx2, ty2, labelText.c_str());
                procLatex->SetTextColor(trueColor);
                procLatex->SetTextSize(0.04);
                procLatex->Draw("SAME");
            }

            if (_k0_trueParentPDG[i] != 0 && _k0_trueParentStartPos[i][0] > -900 && _k0_trueParentEndPos[i][0] > -900) {
                Float_t px1, py1, px2, py2;
                if (projection_type == "XY") {
                    px1 = _k0_trueParentStartPos[i][0]; py1 = _k0_trueParentStartPos[i][1];
                    px2 = _k0_trueParentEndPos[i][0]; py2 = _k0_trueParentEndPos[i][1];
                } else if (projection_type == "XZ") {
                    px1 = _k0_trueParentStartPos[i][0]; py1 = _k0_trueParentStartPos[i][2];
                    px2 = _k0_trueParentEndPos[i][0]; py2 = _k0_trueParentEndPos[i][2];
                } else {
                    px1 = _k0_trueParentStartPos[i][1]; py1 = _k0_trueParentStartPos[i][2];
                    px2 = _k0_trueParentEndPos[i][1]; py2 = _k0_trueParentEndPos[i][2];
                }
                Int_t parentColor = GetParticleColor(_k0_trueParentPDG[i]);
                if (parentColor == kBlack) parentColor = kGray + 1;
                TLine* parentLine = new TLine(px1, py1, px2, py2);
                parentLine->SetLineColor(parentColor);
                parentLine->SetLineStyle(1);
                parentLine->SetLineWidth(3);
                parentLine->Draw("SAME");
            }

            for (Int_t d = 0; d < _k0_trueParentNDaughters[i] && d < kMaxTrueDaughters; ++d) {
                Float_t sx = _k0_trueParentDaughterStartPos[i][d*3 + 0];
                Float_t sy = _k0_trueParentDaughterStartPos[i][d*3 + 1];
                Float_t sz = _k0_trueParentDaughterStartPos[i][d*3 + 2];
                Float_t ex = _k0_trueParentDaughterEndPos[i][d*3 + 0];
                Float_t ey = _k0_trueParentDaughterEndPos[i][d*3 + 1];
                Float_t ez = _k0_trueParentDaughterEndPos[i][d*3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Float_t dx1, dy1, dx2, dy2;
                if (projection_type == "XY") {
                    dx1 = sx; dy1 = sy;
                    dx2 = ex; dy2 = ey;
                } else if (projection_type == "XZ") {
                    dx1 = sx; dy1 = sz;
                    dx2 = ex; dy2 = ez;
                } else {
                    dx1 = sy; dy1 = sz;
                    dx2 = ey; dy2 = ez;
                }
                Int_t daughterColor = GetParticleColor(_k0_trueParentDaughterPDG[i][d]);
                if (daughterColor == kBlack) daughterColor = kGray + 1;
                TLine* daughterLine = new TLine(dx1, dy1, dx2, dy2);
                daughterLine->SetLineColor(daughterColor);
                daughterLine->SetLineStyle(1);
                daughterLine->SetLineWidth(3);
                daughterLine->Draw("SAME");
            }

            for (Int_t s = 0; s < _k0_trueNSiblings[i] && s < kMaxTrueSiblings; ++s) {
                Float_t sx = _k0_trueSiblingStartPos[i][s*3 + 0];
                Float_t sy = _k0_trueSiblingStartPos[i][s*3 + 1];
                Float_t sz = _k0_trueSiblingStartPos[i][s*3 + 2];
                Float_t ex = _k0_trueSiblingEndPos[i][s*3 + 0];
                Float_t ey = _k0_trueSiblingEndPos[i][s*3 + 1];
                Float_t ez = _k0_trueSiblingEndPos[i][s*3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Float_t sx2d, sy2d, ex2d, ey2d;
                if (projection_type == "XY") {
                    sx2d = sx; sy2d = sy;
                    ex2d = ex; ey2d = ey;
                } else if (projection_type == "XZ") {
                    sx2d = sx; sy2d = sz;
                    ex2d = ex; ey2d = ez;
                } else {
                    sx2d = sy; sy2d = sz;
                    ex2d = ey; ey2d = ez;
                }
                Int_t siblingColor = GetParticleColor(_k0_trueSiblingPDG[i][s]);
                if (siblingColor == kBlack) siblingColor = kGray + 1;
                TLine* siblingLine = new TLine(sx2d, sy2d, ex2d, ey2d);
                siblingLine->SetLineColor(siblingColor);
                siblingLine->SetLineStyle(1);
                siblingLine->SetLineWidth(2);
                siblingLine->Draw("SAME");
            }
        }

        if (!_k0_hasTrueObject[i]) {
            if (_k0_trueParentEndPos[i][0] > -900) {
                Float_t cx, cy;
                if (projection_type == "XY") {
                    cx = _k0_trueParentEndPos[i][0];
                    cy = _k0_trueParentEndPos[i][1];
                } else if (projection_type == "XZ") {
                    cx = _k0_trueParentEndPos[i][0];
                    cy = _k0_trueParentEndPos[i][2];
                } else {
                    cx = _k0_trueParentEndPos[i][1];
                    cy = _k0_trueParentEndPos[i][2];
                }
                TMarker* trueCreationMarker = new TMarker(cx, cy, 29);
                trueCreationMarker->SetMarkerColor(kBlue + 1);
                trueCreationMarker->SetMarkerSize(2.0);
                trueCreationMarker->Draw("SAME");
            }

            if (_k0_trueDecayVtxFromRecoDaughters[i][0] > -900) {
                Float_t ax, ay;
                if (projection_type == "XY") {
                    ax = _k0_trueDecayVtxFromRecoDaughters[i][0];
                    ay = _k0_trueDecayVtxFromRecoDaughters[i][1];
                } else if (projection_type == "XZ") {
                    ax = _k0_trueDecayVtxFromRecoDaughters[i][0];
                    ay = _k0_trueDecayVtxFromRecoDaughters[i][2];
                } else {
                    ax = _k0_trueDecayVtxFromRecoDaughters[i][1];
                    ay = _k0_trueDecayVtxFromRecoDaughters[i][2];
                }
                TMarker* trueAnnihMarker = new TMarker(ax, ay, 29);
                trueAnnihMarker->SetMarkerColor(kRed + 1);
                trueAnnihMarker->SetMarkerSize(2.0);
                trueAnnihMarker->Draw("SAME");
            }

            if (_k0_trueParentPDG[i] != 0 && _k0_trueParentStartPos[i][0] > -900 && _k0_trueParentEndPos[i][0] > -900) {
                Float_t px1, py1, px2, py2;
                if (projection_type == "XY") {
                    px1 = _k0_trueParentStartPos[i][0]; py1 = _k0_trueParentStartPos[i][1];
                    px2 = _k0_trueParentEndPos[i][0]; py2 = _k0_trueParentEndPos[i][1];
                } else if (projection_type == "XZ") {
                    px1 = _k0_trueParentStartPos[i][0]; py1 = _k0_trueParentStartPos[i][2];
                    px2 = _k0_trueParentEndPos[i][0]; py2 = _k0_trueParentEndPos[i][2];
                } else {
                    px1 = _k0_trueParentStartPos[i][1]; py1 = _k0_trueParentStartPos[i][2];
                    px2 = _k0_trueParentEndPos[i][1]; py2 = _k0_trueParentEndPos[i][2];
                }
                Int_t parentColor = GetParticleColor(_k0_trueParentPDG[i]);
                if (parentColor == kBlack) parentColor = kGray + 1;
                TLine* parentLine = new TLine(px1, py1, px2, py2);
                parentLine->SetLineColor(parentColor);
                parentLine->SetLineStyle(1);
                parentLine->SetLineWidth(3);
                parentLine->Draw("SAME");
            }

            for (Int_t d = 0; d < _k0_trueParentNDaughters[i] && d < kMaxTrueDaughters; ++d) {
                Float_t sx = _k0_trueParentDaughterStartPos[i][d * 3 + 0];
                Float_t sy = _k0_trueParentDaughterStartPos[i][d * 3 + 1];
                Float_t sz = _k0_trueParentDaughterStartPos[i][d * 3 + 2];
                Float_t ex = _k0_trueParentDaughterEndPos[i][d * 3 + 0];
                Float_t ey = _k0_trueParentDaughterEndPos[i][d * 3 + 1];
                Float_t ez = _k0_trueParentDaughterEndPos[i][d * 3 + 2];
                if (sx <= -900 || ex <= -900) continue;
                Float_t dx1, dy1, dx2, dy2;
                if (projection_type == "XY") {
                    dx1 = sx; dy1 = sy;
                    dx2 = ex; dy2 = ey;
                } else if (projection_type == "XZ") {
                    dx1 = sx; dy1 = sz;
                    dx2 = ex; dy2 = ez;
                } else {
                    dx1 = sy; dy1 = sz;
                    dx2 = ey; dy2 = ez;
                }
                Int_t daughterColor = GetParticleColor(_k0_trueParentDaughterPDG[i][d]);
                if (daughterColor == kBlack) daughterColor = kGray + 1;
                TLine* daughterLine = new TLine(dx1, dy1, dx2, dy2);
                daughterLine->SetLineColor(daughterColor);
                daughterLine->SetLineStyle(1);
                daughterLine->SetLineWidth(3);
                daughterLine->Draw("SAME");
            }
        }

        // NOTE: Old fitted line code removed - now using particle positions below

        // Draw closest points from Pandora-line method
        if (_k0_annVtx_pandoraClosestPt1[i][0] > -900) {
            Float_t pt1_x, pt1_y;
            if (projection_type == "XY") {
                pt1_x = _k0_annVtx_pandoraClosestPt1[i][0];
                pt1_y = _k0_annVtx_pandoraClosestPt1[i][1];
            } else if (projection_type == "XZ") {
                pt1_x = _k0_annVtx_pandoraClosestPt1[i][0];
                pt1_y = _k0_annVtx_pandoraClosestPt1[i][2];
            } else if (projection_type == "YZ") {
                pt1_x = _k0_annVtx_pandoraClosestPt1[i][1];
                pt1_y = _k0_annVtx_pandoraClosestPt1[i][2];
            } else {
                continue;
            }
            TMarker* minDistPt1Pand = new TMarker(pt1_x, pt1_y, 25); // Open square
            minDistPt1Pand->SetMarkerSize(1.0);
            minDistPt1Pand->SetMarkerColor(kAzure + 2);
            minDistPt1Pand->Draw("SAME");
        }

        if (_k0_annVtx_pandoraClosestPt2[i][0] > -900) {
            Float_t pt2_x, pt2_y;
            if (projection_type == "XY") {
                pt2_x = _k0_annVtx_pandoraClosestPt2[i][0];
                pt2_y = _k0_annVtx_pandoraClosestPt2[i][1];
            } else if (projection_type == "XZ") {
                pt2_x = _k0_annVtx_pandoraClosestPt2[i][0];
                pt2_y = _k0_annVtx_pandoraClosestPt2[i][2];
            } else if (projection_type == "YZ") {
                pt2_x = _k0_annVtx_pandoraClosestPt2[i][1];
                pt2_y = _k0_annVtx_pandoraClosestPt2[i][2];
            } else {
                continue;
            }
            TMarker* minDistPt2Pand = new TMarker(pt2_x, pt2_y, 25); // Open square
            minDistPt2Pand->SetMarkerSize(1.0);
            minDistPt2Pand->SetMarkerColor(kAzure + 2);
            minDistPt2Pand->Draw("SAME");
        }

        if (_k0_annVtx_pandoraClosestPt1[i][0] > -900 && _k0_annVtx_pandoraClosestPt2[i][0] > -900) {
            Float_t p1x, p1y, p2x, p2y;
            if (projection_type == "XY") {
                p1x = _k0_annVtx_pandoraClosestPt1[i][0];
                p1y = _k0_annVtx_pandoraClosestPt1[i][1];
                p2x = _k0_annVtx_pandoraClosestPt2[i][0];
                p2y = _k0_annVtx_pandoraClosestPt2[i][1];
            } else if (projection_type == "XZ") {
                p1x = _k0_annVtx_pandoraClosestPt1[i][0];
                p1y = _k0_annVtx_pandoraClosestPt1[i][2];
                p2x = _k0_annVtx_pandoraClosestPt2[i][0];
                p2y = _k0_annVtx_pandoraClosestPt2[i][2];
            } else if (projection_type == "YZ") {
                p1x = _k0_annVtx_pandoraClosestPt1[i][1];
                p1y = _k0_annVtx_pandoraClosestPt1[i][2];
                p2x = _k0_annVtx_pandoraClosestPt2[i][1];
                p2y = _k0_annVtx_pandoraClosestPt2[i][2];
            } else {
                continue;
            }
            TLine* pandConnect = new TLine(p1x, p1y, p2x, p2y);
            pandConnect->SetLineColor(kAzure + 2);
            pandConnect->SetLineStyle(1);
            pandConnect->SetLineWidth(2);
            pandConnect->Draw("SAME");
        }

        // Draw closest points on fitted lines (points of minimum distance)
        if (_k0_annVtx_closestPt1[i][0] > -900) {
            Float_t pt1_x, pt1_y;
            if (projection_type == "XY") {
                pt1_x = _k0_annVtx_closestPt1[i][0];
                pt1_y = _k0_annVtx_closestPt1[i][1];
            } else if (projection_type == "XZ") {
                pt1_x = _k0_annVtx_closestPt1[i][0];
                pt1_y = _k0_annVtx_closestPt1[i][2];
            } else if (projection_type == "YZ") {
                pt1_x = _k0_annVtx_closestPt1[i][1];
                pt1_y = _k0_annVtx_closestPt1[i][2];
            } else {
                continue;
            }

            TMarker* minDistPt1 = new TMarker(pt1_x, pt1_y, 24); // Open circle
            minDistPt1->SetMarkerSize(1.0);
            minDistPt1->SetMarkerColor(kGreen+2);
            minDistPt1->Draw("SAME");
        }

        if (_k0_annVtx_closestPt2[i][0] > -900) {
            Float_t pt2_x, pt2_y;
            if (projection_type == "XY") {
                pt2_x = _k0_annVtx_closestPt2[i][0];
                pt2_y = _k0_annVtx_closestPt2[i][1];
            } else if (projection_type == "XZ") {
                pt2_x = _k0_annVtx_closestPt2[i][0];
                pt2_y = _k0_annVtx_closestPt2[i][2];
            } else if (projection_type == "YZ") {
                pt2_x = _k0_annVtx_closestPt2[i][1];
                pt2_y = _k0_annVtx_closestPt2[i][2];
            } else {
                continue;
            }

            TMarker* minDistPt2 = new TMarker(pt2_x, pt2_y, 24); // Open circle
            minDistPt2->SetMarkerSize(1.0);
            minDistPt2->SetMarkerColor(kGreen+2);
            minDistPt2->Draw("SAME");
        }
    }

    // Draw creation vertex fitted lines and closest points
    for (Int_t i = 0; i < _nK0Candidates && i < kMaxK0; i++) {
        Float_t fitLength = _k0_fitLineLength[i];

        if (fitLength > 0 && _k0_creationVtx_fitLineBeamStart[i][0] > -900 &&
            _k0_creationVtx_closestPtBeam[i][0] > -900) {
            // Get beam particle color
            Int_t beamID = _k0_parentID[i];
            Int_t beamPDG = -1;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == beamID) {
                    beamPDG = _particle_PDG[p];
                    break;
                }
            }
            Int_t beamColor = (beamPDG > 0) ? GetParticleColor(beamPDG) : kBlack;

            // Second particle color (from PDG)
            Int_t secondID = _k0_secondParticleID[i];
            Int_t secondPDG = -1;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == secondID) {
                    secondPDG = _particle_PDG[p];
                    break;
                }
            }
            Int_t secondColor = (secondPDG > 0) ? GetParticleColor(secondPDG) : kGray;

            // Get beam particle actual endPos and verify it has hits (2D)
            Float_t beamPosX = -999, beamPosY = -999, beamPosZ = -999;
            Float_t beamStartX = -999, beamStartY = -999, beamStartZ = -999;
            Int_t beamHits = 0;
            for (Int_t p = 0; p < _nParticles; p++) {
                if (_particle_uniqueID[p] == beamID) {
                    beamPosX = _particle_endPos[p][0];  // Beam uses END position
                    beamPosY = _particle_endPos[p][1];
                    beamPosZ = _particle_endPos[p][2];
                    beamStartX = _particle_startPos[p][0];
                    beamStartY = _particle_startPos[p][1];
                    beamStartZ = _particle_startPos[p][2];
                    beamHits = _particle_nHits[p];
                    break;
                }
            }

            // Only draw if beam particle position was found and particle has hits
            if (beamPosX < -900 || beamHits == 0) continue;

            // Beam fitted line: draw bidirectionally with span >= reconstructed particle length.
            Float_t closestBeamX = _k0_creationVtx_closestPtBeam[i][0];
            Float_t closestBeamY = _k0_creationVtx_closestPtBeam[i][1];
            Float_t closestBeamZ = _k0_creationVtx_closestPtBeam[i][2];
            Float_t beamDirExtX = closestBeamX - beamPosX;
            Float_t beamDirExtY = closestBeamY - beamPosY;
            Float_t beamDirExtZ = closestBeamZ - beamPosZ;
            Float_t beamMag = sqrt(beamDirExtX*beamDirExtX + beamDirExtY*beamDirExtY + beamDirExtZ*beamDirExtZ);
            if (beamMag > 0) {
                beamDirExtX /= beamMag;
                beamDirExtY /= beamMag;
                beamDirExtZ /= beamMag;
            }
            Float_t beamRecoLength = 0.f;
            if (beamStartX > -900.f && beamStartY > -900.f && beamStartZ > -900.f) {
                beamRecoLength = TVector3(beamPosX - beamStartX, beamPosY - beamStartY, beamPosZ - beamStartZ).Mag();
            }
            const Float_t beamSpan = static_cast<Float_t>(std::max<double>(100.0, std::max<double>(beamRecoLength, fitLength)));

            Float_t beamX1, beamY1, beamX2, beamY2;
            if (projection_type == "XY") {
                beamX1 = beamPosX - beamSpan * beamDirExtX;
                beamY1 = beamPosY - beamSpan * beamDirExtY;
                beamX2 = beamPosX + beamSpan * beamDirExtX;
                beamY2 = beamPosY + beamSpan * beamDirExtY;
            } else if (projection_type == "XZ") {
                beamX1 = beamPosX - beamSpan * beamDirExtX;
                beamY1 = beamPosZ - beamSpan * beamDirExtZ;
                beamX2 = beamPosX + beamSpan * beamDirExtX;
                beamY2 = beamPosZ + beamSpan * beamDirExtZ;
            } else if (projection_type == "YZ") {
                beamX1 = beamPosY - beamSpan * beamDirExtY;
                beamY1 = beamPosZ - beamSpan * beamDirExtZ;
                beamX2 = beamPosY + beamSpan * beamDirExtY;
                beamY2 = beamPosZ + beamSpan * beamDirExtZ;
            } else {
                continue;
            }

            TLine* beamLine = new TLine(beamX1, beamY1, beamX2, beamY2);
            beamLine->SetLineColor(beamColor);
            beamLine->SetLineStyle(2); // Dashed (beam uses endPos/endDir)
            beamLine->SetLineWidth(2);
            beamLine->Draw("SAME");

            // Second particle fitted line (uses startPos/startDir)
            if (_k0_creationVtx_fitLineSecondStart[i][0] > -900 &&
                _k0_creationVtx_closestPtSecond[i][0] > -900) {
                // Get second particle actual startPos and verify it has hits (2D)
                Float_t secondPosX = -999, secondPosY = -999, secondPosZ = -999;
                Float_t secondEndX = -999, secondEndY = -999, secondEndZ = -999;
                Int_t secondHits = 0;
                for (Int_t p = 0; p < _nParticles; p++) {
                    if (_particle_uniqueID[p] == secondID) {
                        secondPosX = _particle_startPos[p][0];  // Second uses START position
                        secondPosY = _particle_startPos[p][1];
                        secondPosZ = _particle_startPos[p][2];
                        secondEndX = _particle_endPos[p][0];
                        secondEndY = _particle_endPos[p][1];
                        secondEndZ = _particle_endPos[p][2];
                        secondHits = _particle_nHits[p];
                        break;
                    }
                }

                // Only draw if second particle position was found and particle has hits
                if (secondPosX < -900 || secondHits == 0) continue;

                // Draw bidirectional fitted line with span >= reconstructed particle length.
                Float_t closestSecondX = _k0_creationVtx_closestPtSecond[i][0];
                Float_t closestSecondY = _k0_creationVtx_closestPtSecond[i][1];
                Float_t closestSecondZ = _k0_creationVtx_closestPtSecond[i][2];
                Float_t secondDirExtX = closestSecondX - secondPosX;
                Float_t secondDirExtY = closestSecondY - secondPosY;
                Float_t secondDirExtZ = closestSecondZ - secondPosZ;
                Float_t secondMag = sqrt(secondDirExtX*secondDirExtX + secondDirExtY*secondDirExtY + secondDirExtZ*secondDirExtZ);
                if (secondMag > 0) {
                    secondDirExtX /= secondMag;
                    secondDirExtY /= secondMag;
                    secondDirExtZ /= secondMag;
                }
                Float_t secondRecoLength = 0.f;
                if (secondEndX > -900.f && secondEndY > -900.f && secondEndZ > -900.f) {
                    secondRecoLength = TVector3(secondEndX - secondPosX, secondEndY - secondPosY, secondEndZ - secondPosZ).Mag();
                }
                const Float_t secondSpan = static_cast<Float_t>(std::max<double>(100.0, std::max<double>(secondRecoLength, fitLength)));

                Float_t secondX1 = 0.f, secondY1 = 0.f, secondX2 = 0.f, secondY2 = 0.f;
                if (projection_type == "XY") {
                    secondX1 = secondPosX - secondSpan * secondDirExtX;
                    secondY1 = secondPosY - secondSpan * secondDirExtY;
                    secondX2 = secondPosX + secondSpan * secondDirExtX;
                    secondY2 = secondPosY + secondSpan * secondDirExtY;
                } else if (projection_type == "XZ") {
                    secondX1 = secondPosX - secondSpan * secondDirExtX;
                    secondY1 = secondPosZ - secondSpan * secondDirExtZ;
                    secondX2 = secondPosX + secondSpan * secondDirExtX;
                    secondY2 = secondPosZ + secondSpan * secondDirExtZ;
                } else if (projection_type == "YZ") {
                    secondX1 = secondPosY - secondSpan * secondDirExtY;
                    secondY1 = secondPosZ - secondSpan * secondDirExtZ;
                    secondX2 = secondPosY + secondSpan * secondDirExtY;
                    secondY2 = secondPosZ + secondSpan * secondDirExtZ;
                }

                TLine* secondLine = new TLine(secondX1, secondY1, secondX2, secondY2);
                secondLine->SetLineColor(secondColor);
                secondLine->SetLineStyle(2); // Dashed (uses startPos/startDir)
                secondLine->SetLineWidth(2);
                secondLine->Draw("SAME");
            }

            // Draw closest points
            if (_k0_creationVtx_closestPtBeam[i][0] > -900) {
                Float_t cpX, cpY;
                if (projection_type == "XY") {
                    cpX = _k0_creationVtx_closestPtBeam[i][0];
                    cpY = _k0_creationVtx_closestPtBeam[i][1];
                } else if (projection_type == "XZ") {
                    cpX = _k0_creationVtx_closestPtBeam[i][0];
                    cpY = _k0_creationVtx_closestPtBeam[i][2];
                } else if (projection_type == "YZ") {
                    cpX = _k0_creationVtx_closestPtBeam[i][1];
                    cpY = _k0_creationVtx_closestPtBeam[i][2];
                } else {
                    continue;
                }

                TMarker* marker = new TMarker(cpX, cpY, 24); // Open circle
                marker->SetMarkerColor(kOrange);
                marker->SetMarkerSize(1.5);
                marker->Draw("SAME");
            }

            if (_k0_creationVtx_closestPtSecond[i][0] > -900) {
                Float_t cpX, cpY;
                if (projection_type == "XY") {
                    cpX = _k0_creationVtx_closestPtSecond[i][0];
                    cpY = _k0_creationVtx_closestPtSecond[i][1];
                } else if (projection_type == "XZ") {
                    cpX = _k0_creationVtx_closestPtSecond[i][0];
                    cpY = _k0_creationVtx_closestPtSecond[i][2];
                } else if (projection_type == "YZ") {
                    cpX = _k0_creationVtx_closestPtSecond[i][1];
                    cpY = _k0_creationVtx_closestPtSecond[i][2];
                } else {
                    continue;
                }

                TMarker* marker = new TMarker(cpX, cpY, 24); // Open circle
                marker->SetMarkerColor(kOrange);
                marker->SetMarkerSize(1.5);
                marker->Draw("SAME");
            }

            // Draw white dotted line connecting the two closest points (creation vertex)
            if (_k0_creationVtx_closestPtBeam[i][0] > -900 && _k0_creationVtx_closestPtSecond[i][0] > -900) {
                Float_t cp1X, cp1Y, cp2X, cp2Y;
                if (projection_type == "XY") {
                    cp1X = _k0_creationVtx_closestPtBeam[i][0];
                    cp1Y = _k0_creationVtx_closestPtBeam[i][1];
                    cp2X = _k0_creationVtx_closestPtSecond[i][0];
                    cp2Y = _k0_creationVtx_closestPtSecond[i][1];
                } else if (projection_type == "XZ") {
                    cp1X = _k0_creationVtx_closestPtBeam[i][0];
                    cp1Y = _k0_creationVtx_closestPtBeam[i][2];
                    cp2X = _k0_creationVtx_closestPtSecond[i][0];
                    cp2Y = _k0_creationVtx_closestPtSecond[i][2];
                } else if (projection_type == "YZ") {
                    cp1X = _k0_creationVtx_closestPtBeam[i][1];
                    cp1Y = _k0_creationVtx_closestPtBeam[i][2];
                    cp2X = _k0_creationVtx_closestPtSecond[i][1];
                    cp2Y = _k0_creationVtx_closestPtSecond[i][2];
                } else {
                    continue;
                }

                TLine* connectLine = new TLine(cp1X, cp1Y, cp2X, cp2Y);
                connectLine->SetLineColor(kWhite);
                connectLine->SetLineStyle(3); // Dotted
                connectLine->SetLineWidth(2);
                connectLine->Draw("SAME");
            }
        }
    }

    // Draw annihilation vertex fitted lines and closest points for vertex visualization
    for (Int_t i = 0; i < _nK0Candidates && i < kMaxK0; i++) {
        Float_t fitLength = _k0_fitLineLength[i];
        if (fitLength <= 0) continue;

        // Get daughter colors by matching UniqueID
        Int_t d1Color = kCyan;
        Int_t d2Color = kSpring;
        Int_t d1ID = _k0_daughter1ID[i];
        Int_t d2ID = _k0_daughter2ID[i];

        for (Int_t p = 0; p < _nParticles; p++) {
            if (_particle_uniqueID[p] == d1ID) {
                d1Color = GetParticleColor(_particle_PDG[p]);
            }
            if (_particle_uniqueID[p] == d2ID) {
                d2Color = GetParticleColor(_particle_PDG[p]);
            }
        }

        // Draw fitted line 1 using fitted direction from start position
        Float_t d1PosX = _k0_annVtx_fitLine1Start[i][0];
        Float_t d1PosY = _k0_annVtx_fitLine1Start[i][1];
        Float_t d1PosZ = _k0_annVtx_fitLine1Start[i][2];
        Float_t dirX = _k0_annVtx_fitLine1Dir[i][0];
        Float_t dirY = _k0_annVtx_fitLine1Dir[i][1];
        Float_t dirZ = _k0_annVtx_fitLine1Dir[i][2];

        if (d1PosX > -900 && (dirX*dirX + dirY*dirY + dirZ*dirZ) > 1e-10f) {
            Float_t mag = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
            dirX /= mag;
            dirY /= mag;
            dirZ /= mag;

            TVector3 s3(d1PosX, d1PosY, d1PosZ);
            TVector3 d3(dirX, dirY, dirZ);
            TVector3 base3 = s3;
            double d1RecoLength = 0.0;
            for (Int_t p = 0; p < _nParticles; ++p) {
                if (_particle_uniqueID[p] != d1ID) continue;
                const Float_t sx = _particle_startPos[p][0];
                const Float_t sy = _particle_startPos[p][1];
                const Float_t sz = _particle_startPos[p][2];
                const Float_t ex = _particle_endPos[p][0];
                const Float_t ey = _particle_endPos[p][1];
                const Float_t ez = _particle_endPos[p][2];
                if (sx > -900.f && sy > -900.f && sz > -900.f && ex > -900.f && ey > -900.f && ez > -900.f) {
                    d1RecoLength = TVector3(ex - sx, ey - sy, ez - sz).Mag();
                }
                break;
            }
            const double span3 = std::max<double>(100.0, std::max<double>(d1RecoLength, fitLength));
            double tMin3 = -span3;
            double tMax3 = span3;
            if (_k0_annVtx_closestPt1[i][0] > -900) {
                TVector3 cp3(_k0_annVtx_closestPt1[i][0], _k0_annVtx_closestPt1[i][1], _k0_annVtx_closestPt1[i][2]);
                const TVector3 anchorToClosest3 = cp3 - s3;
                if (anchorToClosest3.Mag2() > 1e-8) {
                    d3 = anchorToClosest3.Unit();
                }
                base3 = cp3; // Guarantee rendered fit line passes through closest point.
                const double tAnchor3 = (s3 - base3).Dot(d3);
                tMin3 = std::min(tMin3, tAnchor3 - span3);
                tMax3 = std::max(tMax3, tAnchor3 + span3);
            }
            TVector3 lineStart3 = base3 + tMin3 * d3;
            TVector3 lineEnd3 = base3 + tMax3 * d3;

            Float_t line1_x1, line1_y1, line1_x2, line1_y2;

            if (projection_type == "XY") {
                line1_x1 = lineStart3.X(); line1_y1 = lineStart3.Y();
                line1_x2 = lineEnd3.X(); line1_y2 = lineEnd3.Y();
            } else if (projection_type == "XZ") {
                line1_x1 = lineStart3.X(); line1_y1 = lineStart3.Z();
                line1_x2 = lineEnd3.X(); line1_y2 = lineEnd3.Z();
            } else if (projection_type == "YZ") {
                line1_x1 = lineStart3.Y(); line1_y1 = lineStart3.Z();
                line1_x2 = lineEnd3.Y(); line1_y2 = lineEnd3.Z();
            } else {
                continue;
            }

            TLine* fitLine1 = new TLine(line1_x1, line1_y1, line1_x2, line1_y2);
            fitLine1->SetLineColor(d1Color);
            fitLine1->SetLineStyle(2); // Dashed (daughter uses startPos/startDir)
            fitLine1->SetLineWidth(1);
            fitLine1->Draw("SAME");

            Float_t anchor1_x, anchor1_y;
            if (projection_type == "XY") {
                anchor1_x = d1PosX; anchor1_y = d1PosY;
            } else if (projection_type == "XZ") {
                anchor1_x = d1PosX; anchor1_y = d1PosZ;
            } else {
                anchor1_x = d1PosY; anchor1_y = d1PosZ;
            }
            TMarker* anchor1 = new TMarker(anchor1_x, anchor1_y, 29);
            anchor1->SetMarkerColor(kYellow + 1);
            anchor1->SetMarkerSize(1.6);
            anchor1->Draw("SAME");
        }

        // Draw fitted line 2 using fitted direction from start position
        Float_t d2PosX = _k0_annVtx_fitLine2Start[i][0];
        Float_t d2PosY = _k0_annVtx_fitLine2Start[i][1];
        Float_t d2PosZ = _k0_annVtx_fitLine2Start[i][2];
        Float_t dir2X = _k0_annVtx_fitLine2Dir[i][0];
        Float_t dir2Y = _k0_annVtx_fitLine2Dir[i][1];
        Float_t dir2Z = _k0_annVtx_fitLine2Dir[i][2];

        if (d2PosX > -900 && (dir2X*dir2X + dir2Y*dir2Y + dir2Z*dir2Z) > 1e-10f) {
            Float_t mag2 = sqrt(dir2X*dir2X + dir2Y*dir2Y + dir2Z*dir2Z);
            dir2X /= mag2;
            dir2Y /= mag2;
            dir2Z /= mag2;

            TVector3 s3(d2PosX, d2PosY, d2PosZ);
            TVector3 d3(dir2X, dir2Y, dir2Z);
            TVector3 base3 = s3;
            double d2RecoLength = 0.0;
            for (Int_t p = 0; p < _nParticles; ++p) {
                if (_particle_uniqueID[p] != d2ID) continue;
                const Float_t sx = _particle_startPos[p][0];
                const Float_t sy = _particle_startPos[p][1];
                const Float_t sz = _particle_startPos[p][2];
                const Float_t ex = _particle_endPos[p][0];
                const Float_t ey = _particle_endPos[p][1];
                const Float_t ez = _particle_endPos[p][2];
                if (sx > -900.f && sy > -900.f && sz > -900.f && ex > -900.f && ey > -900.f && ez > -900.f) {
                    d2RecoLength = TVector3(ex - sx, ey - sy, ez - sz).Mag();
                }
                break;
            }
            const double span3 = std::max<double>(100.0, std::max<double>(d2RecoLength, fitLength));
            double tMin3 = -span3;
            double tMax3 = span3;
            if (_k0_annVtx_closestPt2[i][0] > -900) {
                TVector3 cp3(_k0_annVtx_closestPt2[i][0], _k0_annVtx_closestPt2[i][1], _k0_annVtx_closestPt2[i][2]);
                const TVector3 anchorToClosest3 = cp3 - s3;
                if (anchorToClosest3.Mag2() > 1e-8) {
                    d3 = anchorToClosest3.Unit();
                }
                base3 = cp3; // Guarantee rendered fit line passes through closest point.
                const double tAnchor3 = (s3 - base3).Dot(d3);
                tMin3 = std::min(tMin3, tAnchor3 - span3);
                tMax3 = std::max(tMax3, tAnchor3 + span3);
            }
            TVector3 lineStart3 = base3 + tMin3 * d3;
            TVector3 lineEnd3 = base3 + tMax3 * d3;

            Float_t line2_x1, line2_y1, line2_x2, line2_y2;

            if (projection_type == "XY") {
                line2_x1 = lineStart3.X(); line2_y1 = lineStart3.Y();
                line2_x2 = lineEnd3.X(); line2_y2 = lineEnd3.Y();
            } else if (projection_type == "XZ") {
                line2_x1 = lineStart3.X(); line2_y1 = lineStart3.Z();
                line2_x2 = lineEnd3.X(); line2_y2 = lineEnd3.Z();
            } else if (projection_type == "YZ") {
                line2_x1 = lineStart3.Y(); line2_y1 = lineStart3.Z();
                line2_x2 = lineEnd3.Y(); line2_y2 = lineEnd3.Z();
            } else {
                continue;
            }

            TLine* fitLine2 = new TLine(line2_x1, line2_y1, line2_x2, line2_y2);
            fitLine2->SetLineColor(d2Color);
            fitLine2->SetLineStyle(2); // Dashed (daughter uses startPos/startDir)
            fitLine2->SetLineWidth(1);
            fitLine2->Draw("SAME");

            Float_t anchor2_x, anchor2_y;
            if (projection_type == "XY") {
                anchor2_x = d2PosX; anchor2_y = d2PosY;
            } else if (projection_type == "XZ") {
                anchor2_x = d2PosX; anchor2_y = d2PosZ;
            } else {
                anchor2_x = d2PosY; anchor2_y = d2PosZ;
            }
            TMarker* anchor2 = new TMarker(anchor2_x, anchor2_y, 29);
            anchor2->SetMarkerColor(kYellow + 1);
            anchor2->SetMarkerSize(1.6);
            anchor2->Draw("SAME");
        }

        // Draw reconstructed vertex-momentum and true-K0 direction (respects "K0 Momentum Arrows" visibility).
        const Float_t annX = _k0_annihilationVtxPos[i][0];
        const Float_t annY = _k0_annihilationVtxPos[i][1];
        const Float_t annZ = _k0_annihilationVtxPos[i][2];
        if (showMomentumArrows && annX > -900 && annY > -900 && annZ > -900) {
            const Float_t vecLen = std::max(12.f, std::min(38.f, _k0_annihilationVtxRadius[i] * 0.5f));
            auto drawProjectedDir = [&](const Float_t dirArr[3], Float_t ax, Float_t ay, Float_t az, Int_t color, Int_t style, const char* name) {
                if (dirArr[0] <= -900) return;
                if (ax <= -900 || ay <= -900 || az <= -900) return;
                TVector3 d(dirArr[0], dirArr[1], dirArr[2]);
                if (d.Mag2() <= 1e-10) return;
                d = d.Unit();
                const TVector3 p0(ax, ay, az);
                const TVector3 p1 = p0 + vecLen * d;
                Float_t x0, y0, x1, y1;
                if (projection_type == "XY") {
                    x0 = p0.X(); y0 = p0.Y(); x1 = p1.X(); y1 = p1.Y();
                } else if (projection_type == "XZ") {
                    x0 = p0.X(); y0 = p0.Z(); x1 = p1.X(); y1 = p1.Z();
                } else if (projection_type == "YZ") {
                    x0 = p0.Y(); y0 = p0.Z(); x1 = p1.Y(); y1 = p1.Z();
                } else {
                    return;
                }
                TArrow* arrow = new TArrow(x0, y0, x1, y1, 0.01, "|>");
                arrow->SetLineColor(color);
                arrow->SetFillColor(color);
                arrow->SetLineStyle(style);
                arrow->SetLineWidth(2);
                arrow->Draw("SAME");
                (void)name;
            };
            drawProjectedDir(_k0_annVtx_momentumDirFit[i], annX, annY, annZ, kMagenta + 1, 1, "RecoVtxMomDir");
            if (_k0_hasTrueObject[i]) {
                Float_t trueMcVtx[3];
                bool haveTrueMcVtx = false;
                if (_k0_trueDecayVtxFromRecoDaughters[i][0] > -900.f &&
                    _k0_trueDecayVtxFromRecoDaughters[i][1] > -900.f &&
                    _k0_trueDecayVtxFromRecoDaughters[i][2] > -900.f) {
                    trueMcVtx[0] = _k0_trueDecayVtxFromRecoDaughters[i][0];
                    trueMcVtx[1] = _k0_trueDecayVtxFromRecoDaughters[i][1];
                    trueMcVtx[2] = _k0_trueDecayVtxFromRecoDaughters[i][2];
                    haveTrueMcVtx = true;
                }
                if (!haveTrueMcVtx) {
                    haveTrueMcVtx = TrueK0McDecayVertex(trueMcVtx, _k0_trueNDaughters[i], _k0_trueDaughterPDG[i],
                                                        _k0_trueDaughterStartPos[i]);
                }
                if (!haveTrueMcVtx && _k0_trueEndPos[i][0] > -900 && _k0_trueEndPos[i][1] > -900 &&
                    _k0_trueEndPos[i][2] > -900) {
                    trueMcVtx[0] = _k0_trueEndPos[i][0];
                    trueMcVtx[1] = _k0_trueEndPos[i][1];
                    trueMcVtx[2] = _k0_trueEndPos[i][2];
                    haveTrueMcVtx = true;
                }
                if (haveTrueMcVtx) {
                    drawProjectedDir(_k0_trueK0Dir[i], trueMcVtx[0], trueMcVtx[1], trueMcVtx[2], kBlue + 1, 2,
                                     "TrueK0Dir");
                }
            }
        }

        // Draw closest points on annihilation fitted lines (points of minimum distance)
        Float_t closestPt1X = _k0_annVtx_closestPt1[i][0];
        Float_t closestPt1Y = _k0_annVtx_closestPt1[i][1];
        Float_t closestPt1Z = _k0_annVtx_closestPt1[i][2];
        Float_t closestPt2X = _k0_annVtx_closestPt2[i][0];
        Float_t closestPt2Y = _k0_annVtx_closestPt2[i][1];
        Float_t closestPt2Z = _k0_annVtx_closestPt2[i][2];

        if (closestPt1X > -900) {
            Float_t pt1_x, pt1_y;
            if (projection_type == "XY") {
                pt1_x = closestPt1X; pt1_y = closestPt1Y;
            } else if (projection_type == "XZ") {
                pt1_x = closestPt1X; pt1_y = closestPt1Z;
            } else if (projection_type == "YZ") {
                pt1_x = closestPt1Y; pt1_y = closestPt1Z;
            } else {
                continue;
            }

            TMarker* marker1 = new TMarker(pt1_x, pt1_y, 24); // Open circle
            marker1->SetMarkerColor(kOrange);
            marker1->SetMarkerSize(1.5);
            marker1->Draw("SAME");
        }

        if (closestPt2X > -900) {
            Float_t pt2_x, pt2_y;
            if (projection_type == "XY") {
                pt2_x = closestPt2X; pt2_y = closestPt2Y;
            } else if (projection_type == "XZ") {
                pt2_x = closestPt2X; pt2_y = closestPt2Z;
            } else if (projection_type == "YZ") {
                pt2_x = closestPt2Y; pt2_y = closestPt2Z;
            } else {
                continue;
            }

            TMarker* marker2 = new TMarker(pt2_x, pt2_y, 24); // Open circle
            marker2->SetMarkerColor(kOrange);
            marker2->SetMarkerSize(1.5);
            marker2->Draw("SAME");
        }

        // Draw white dotted line connecting the two closest points (annihilation vertex)
        if (closestPt1X > -900 && closestPt2X > -900) {
            Float_t cp1X, cp1Y, cp2X, cp2Y;
            if (projection_type == "XY") {
                cp1X = closestPt1X; cp1Y = closestPt1Y;
                cp2X = closestPt2X; cp2Y = closestPt2Y;
            } else if (projection_type == "XZ") {
                cp1X = closestPt1X; cp1Y = closestPt1Z;
                cp2X = closestPt2X; cp2Y = closestPt2Z;
            } else if (projection_type == "YZ") {
                cp1X = closestPt1Y; cp1Y = closestPt1Z;
                cp2X = closestPt2Y; cp2Y = closestPt2Z;
            } else {
                continue;
            }

            TLine* connectLine = new TLine(cp1X, cp1Y, cp2X, cp2Y);
            connectLine->SetLineColor(kWhite);
            connectLine->SetLineStyle(3); // Dotted
            connectLine->SetLineWidth(2);
            connectLine->Draw("SAME");
        }
    }

    if (_nK0Candidates > 0) {
        std::cout << "Drew " << _nK0Candidates << " K0 candidates on " << projection_type << " canvas" << std::endl;
    }

    // Draw all true particle trajectories on 2D canvas
    for (Int_t i = 0; i < _nAllTrueParticles && i < kMaxAllTrueParticles; ++i) {
        if (_allTrueParticle_startPos[i][0] <= -900 || _allTrueParticle_endPos[i][0] <= -900) continue;

        Float_t x1 = 0.f, y1 = 0.f, x2 = 0.f, y2 = 0.f;

        if (projection_type == "XY") {
            x1 = _allTrueParticle_startPos[i][0]; y1 = _allTrueParticle_startPos[i][1];
            x2 = _allTrueParticle_endPos[i][0];   y2 = _allTrueParticle_endPos[i][1];
        } else if (projection_type == "XZ") {
            x1 = _allTrueParticle_startPos[i][0]; y1 = _allTrueParticle_startPos[i][2];
            x2 = _allTrueParticle_endPos[i][0];   y2 = _allTrueParticle_endPos[i][2];
        } else if (projection_type == "YZ") {
            x1 = _allTrueParticle_startPos[i][1]; y1 = _allTrueParticle_startPos[i][2];
            x2 = _allTrueParticle_endPos[i][1];   y2 = _allTrueParticle_endPos[i][2];
        } else {
            continue;
        }

        Int_t trueColor = GetParticleColor(_allTrueParticle_PDG[i]);
        if (trueColor == kBlack) trueColor = kGray + 1;

        if (showStandaloneTruth) {
            TLine* trueLine = new TLine(x1, y1, x2, y2);
            trueLine->SetLineColor(trueColor);
            trueLine->SetLineStyle(1);
            trueLine->SetLineWidth(3);
            trueLine->Draw("SAME");

            TMarker* trueStartMarker = new TMarker(x1, y1, 29);
            trueStartMarker->SetMarkerColor(trueColor);
            trueStartMarker->SetMarkerSize(2.0);
            trueStartMarker->Draw("SAME");

            TMarker* trueEndMarker = new TMarker(x2, y2, 29);
            trueEndMarker->SetMarkerColor(trueColor);
            trueEndMarker->SetMarkerSize(2.0);
            trueEndMarker->Draw("SAME");

            std::string procLabel = ProcessEnumToString(_allTrueParticle_processEnd[i]);
            std::string standaloneLabel;
            if (!procLabel.empty()) {
                standaloneLabel = Form("%s (PDG=%d)", procLabel.c_str(), _allTrueParticle_PDG[i]);
            } else {
                standaloneLabel = Form("PDG=%d", _allTrueParticle_PDG[i]);
            }
            if (!standaloneLabel.empty()) {
                TLatex* procLatex = new TLatex(x2, y2, Form("#scale[0.6]{%s}", standaloneLabel.c_str()));
                procLatex->SetTextColor(trueColor);
                procLatex->SetTextSize(0.03);
                procLatex->Draw("SAME");
            }
        }
    }
}

