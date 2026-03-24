#ifndef neutralKaonEventDisplay_h
#define neutralKaonEventDisplay_h

#include "pdEventDisplay.hxx"
#include <unordered_map>
#include <vector>

// Forward declarations
class ToyBoxNeutralKaon;
class TH1F;
class TLine;
class TLegend;
class TCanvas;

/// Neutral kaon analysis-specific event display
/// Extends pdEventDisplay with K0 candidate visualization (vertices, radii)
class neutralKaonEventDisplay : public pdEventDisplay {
public:
    neutralKaonEventDisplay();
    virtual ~neutralKaonEventDisplay();

    static constexpr Int_t kParentTrajHistBins = 60;
    static constexpr Float_t kParentTrajHistMin = -1.f;
    static constexpr Float_t kParentTrajHistMax = 1.f;

protected:
    // Override analysis-specific virtual methods from pdEventDisplay

    /// Add K0 candidate variables to EventDisplayData tree
    virtual void AddAnalysisVariables(OutputManager& output, Int_t tree_index) override;

    /// Fill K0 candidate data from ToyBoxNeutralKaon
    virtual void FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box) override;

    /// Read K0 candidate data from tree
    virtual bool ReadAnalysisData(TTree* tree) override;

    /// Draw K0 vertices and trajectories in 3D
    virtual void DrawAnalysisContent3D(TEveScene* scene) override;

    /// Draw K0 vertices in 2D projections
    virtual void DrawAnalysisContent2D(TEveProjectionManager* manager, const std::string& projection_type) override;

    /// Draw K0 vertices on TCanvas 2D view (with circles and degeneracy labels)
    virtual void DrawAnalysisContentCanvas2D(TCanvas* canvas, const std::string& projection_type) override;

private:
    // K0 candidate data (read from tree during visualization)
    static constexpr Int_t kMaxK0 = 10;
    static constexpr Int_t kMaxTrueK0 = 20;
    static constexpr Int_t kMaxTrueDaughters = 8;
    static constexpr Int_t kMaxTrueSiblings = 8;
    Int_t _nK0Candidates;
    Int_t _k0_daughter1ID[kMaxK0];
    Int_t _k0_daughter2ID[kMaxK0];
    Int_t _k0_parentID[kMaxK0];
    Int_t _k0_secondParticleID[kMaxK0];  // Second particle in creation vertex
    Float_t _k0_creationVtxRadius[kMaxK0];
    Float_t _k0_annihilationVtxRadius[kMaxK0];
    Int_t _k0_creationVtxDeg[kMaxK0];
    Int_t _k0_annihilationVtxDeg[kMaxK0];
    Float_t _k0_creationVtxPos[kMaxK0][3];
    Float_t _k0_annihilationVtxPos[kMaxK0][3];
    Float_t _k0_startPos[kMaxK0][3];
    Float_t _k0_endPos[kMaxK0][3];

    // Fitted lines for annihilation vertex
    Float_t _k0_annVtx_fitLine1Start[kMaxK0][3];  // Daughter 1 fit line start
    Float_t _k0_annVtx_fitLine1Dir[kMaxK0][3];    // Daughter 1 fit line direction
    Float_t _k0_annVtx_fitLine2Start[kMaxK0][3];  // Daughter 2 fit line start
    Float_t _k0_annVtx_fitLine2Dir[kMaxK0][3];    // Daughter 2 fit line direction
    Float_t _k0_annVtx_closestPt1[kMaxK0][3];     // Closest point on line 1
    Float_t _k0_annVtx_closestPt2[kMaxK0][3];     // Closest point on line 2

    // Fitted lines for creation vertex
    Float_t _k0_creationVtx_fitLineBeamStart[kMaxK0][3];   // Beam particle fit line start
    Float_t _k0_creationVtx_fitLineBeamDir[kMaxK0][3];     // Beam particle fit line direction
    Float_t _k0_creationVtx_fitLineSecondStart[kMaxK0][3]; // Second particle fit line start
    Float_t _k0_creationVtx_fitLineSecondDir[kMaxK0][3];   // Second particle fit line direction
    Float_t _k0_creationVtx_closestPtBeam[kMaxK0][3];      // Closest point on beam line
    Float_t _k0_creationVtx_closestPtSecond[kMaxK0][3];    // Closest point on second line

    Float_t _k0_fitLineLength[kMaxK0];            // Fit line length for drawing

    // Degeneracy distances (distances from vertex to particles within sphere)
    Float_t _k0_creationVtxDegDist[kMaxK0][5];   // Creation vertex degeneracy distances
    Float_t _k0_annihilationVtxDegDist[kMaxK0][5]; // Annihilation vertex degeneracy distances

    // True K0 trajectory (if TrueObject exists)
    Int_t _k0_hasTrueObject[kMaxK0];              // Has associated true neutral particle
    Float_t _k0_trueStartPos[kMaxK0][3];          // True K0 start position
    Float_t _k0_trueEndPos[kMaxK0][3];            // True K0 end position
    Int_t _k0_truePDG[kMaxK0];                    // True K0 PDG code
    Int_t _k0_trueProcessEnd[kMaxK0];
    Int_t _k0_trueParentPDG[kMaxK0];
    Float_t _k0_trueParentStartPos[kMaxK0][3];
    Float_t _k0_trueParentEndPos[kMaxK0][3];
    Int_t _k0_trueNDaughters[kMaxK0];
    Float_t _k0_trueDaughterStartPos[kMaxK0][kMaxTrueDaughters*3];
    Float_t _k0_trueDaughterEndPos[kMaxK0][kMaxTrueDaughters*3];
    Int_t _k0_trueDaughterPDG[kMaxK0][kMaxTrueDaughters];
    Int_t _k0_trueNSiblings[kMaxK0];
    Float_t _k0_trueSiblingStartPos[kMaxK0][kMaxTrueSiblings*3];
    Float_t _k0_trueSiblingEndPos[kMaxK0][kMaxTrueSiblings*3];
    Int_t _k0_trueSiblingPDG[kMaxK0][kMaxTrueSiblings];
    Float_t _k0_parentStartPos[kMaxK0][3];
    Float_t _k0_parentEndPos[kMaxK0][3];
    Float_t _k0_parentLength[kMaxK0];
    Float_t _k0_parentTrajDir[kMaxK0][3];
    Float_t _k0_parentTrajDirHist[kMaxK0][kParentTrajHistBins*3];
    Int_t _k0_parentTrajDirNPts[kMaxK0];
    Float_t _k0_secondTrajDir[kMaxK0][3];
    Int_t _k0_secondTrajDirNPts[kMaxK0];
    Float_t _k0_dau1TrajDir[kMaxK0][3];
    Int_t _k0_dau1TrajDirNPts[kMaxK0];
    Float_t _k0_dau2TrajDir[kMaxK0][3];
    Int_t _k0_dau2TrajDirNPts[kMaxK0];

    // Standalone true K0 data (without reconstructed neutral particle)
    Int_t _nTrueK0;
    Float_t _trueK0_startPos[kMaxTrueK0][3];
    Float_t _trueK0_endPos[kMaxTrueK0][3];
    Int_t _trueK0_PDG[kMaxTrueK0];
    Int_t _trueK0_processEnd[kMaxTrueK0];
    Int_t _trueK0_parentPDG[kMaxTrueK0];
    Float_t _trueK0_parentStartPos[kMaxTrueK0][3];
    Float_t _trueK0_parentEndPos[kMaxTrueK0][3];
    Int_t _trueK0_nDaughters[kMaxTrueK0];
    Float_t _trueK0_daughterStartPos[kMaxTrueK0][kMaxTrueDaughters*3];
    Float_t _trueK0_daughterEndPos[kMaxTrueK0][kMaxTrueDaughters*3];
    Int_t _trueK0_daughterPDG[kMaxTrueK0][kMaxTrueDaughters];
    Int_t _trueK0_nSiblings[kMaxTrueK0];
    Float_t _trueK0_siblingStartPos[kMaxTrueK0][kMaxTrueSiblings*3];
    Float_t _trueK0_siblingEndPos[kMaxTrueK0][kMaxTrueSiblings*3];
    Int_t _trueK0_siblingPDG[kMaxTrueK0][kMaxTrueSiblings];

    // Variable indices for OutputManager (extend pdEventDisplay enum)
    enum enumNK0EventDisplayVars {
        ednK0Candidates = pdEventDisplay::edmaxvars,
        edk0_daughter1ID,
        edk0_daughter2ID,
        edk0_parentID,
        edk0_secondParticleID,
        edk0_creationVtxRadius,
        edk0_annihilationVtxRadius,
        edk0_creationVtxDeg,
        edk0_annihilationVtxDeg,
        edk0_creationVtxPos,
        edk0_annihilationVtxPos,
        edk0_startPos,
        edk0_endPos,
        edk0_annVtx_fitLine1Start,
        edk0_annVtx_fitLine1Dir,
        edk0_annVtx_fitLine2Start,
        edk0_annVtx_fitLine2Dir,
        edk0_annVtx_closestPt1,
        edk0_annVtx_closestPt2,
        edk0_creationVtx_fitLineBeamStart,
        edk0_creationVtx_fitLineBeamDir,
        edk0_creationVtx_fitLineSecondStart,
        edk0_creationVtx_fitLineSecondDir,
        edk0_creationVtx_closestPtBeam,
        edk0_creationVtx_closestPtSecond,
        edk0_fitLineLength,
        edk0_creationVtxDegDist,
        edk0_annihilationVtxDegDist,
        edk0_hasTrueObject,
        edk0_trueStartPos,
        edk0_trueEndPos,
        edk0_truePDG,
        edk0_trueProcessEnd,
        edk0_trueParentPDG,
        edk0_trueParentStartPos,
        edk0_trueParentEndPos,
        edk0_trueNDaughters,
        edk0_trueDaughterStartPos,
        edk0_trueDaughterEndPos,
        edk0_trueDaughterPDG,
        edk0_trueNSiblings,
        edk0_trueSiblingStartPos,
        edk0_trueSiblingEndPos,
        edk0_trueSiblingPDG,
        edk0_parentStartPos,
        edk0_parentEndPos,
        edk0_parentLength,
        edk0_parentTrajDir,
        edk0_parentTrajDirHist,
        edk0_parentTrajDirNPts,
        edk0_secondTrajDir,
        edk0_secondTrajDirNPts,
        edk0_dau1TrajDir,
        edk0_dau1TrajDirNPts,
        edk0_dau2TrajDir,
        edk0_dau2TrajDirNPts,
        ednTrueK0,
        edtrueK0_startPos,
        edtrueK0_endPos,
        edtrueK0_PDG,
        edtrueK0_processEnd,
        edtrueK0_parentPDG,
        edtrueK0_parentStartPos,
        edtrueK0_parentEndPos,
        edtrueK0_nDaughters,
        edtrueK0_daughterStartPos,
        edtrueK0_daughterEndPos,
        edtrueK0_daughterPDG,
        edtrueK0_nSiblings,
        edtrueK0_siblingStartPos,
        edtrueK0_siblingEndPos,
        edtrueK0_siblingPDG,
        nk0maxvars
    };

    ClassDef(neutralKaonEventDisplay, 0);

    void EnsureSelectionHooks();
    void OnSelectionAdded(TEveElement* element);
    void OnSelectionRepeated(TEveElement* element);
    void ShowParentTrajectoryHistograms(Int_t index);
    void ClearParentDirectionVisuals();

    std::unordered_map<TEveElement*, Int_t> _parentDirElementToIndex;
    std::vector<TEveLine*> _parentDirectionLines;
    Bool_t _parentDirSelectionHooked = kFALSE;
    TCanvas* _parentDirCanvas = nullptr;
    TH1F* _parentDirHists[3] = {nullptr, nullptr, nullptr};
    TLine* _parentDirTrueLines[3] = {nullptr, nullptr, nullptr};
    TLegend* _parentDirLegends[3] = {nullptr, nullptr, nullptr};
};

#endif

