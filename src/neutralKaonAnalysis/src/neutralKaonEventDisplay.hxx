#ifndef neutralKaonEventDisplay_h
#define neutralKaonEventDisplay_h

#include "pdEventDisplay.hxx"

// Forward declarations
class ToyBoxNeutralKaon;

/// Neutral kaon analysis-specific event display
/// Extends pdEventDisplay with K0 candidate visualization (vertices, radii)
class neutralKaonEventDisplay : public pdEventDisplay {
public:
    neutralKaonEventDisplay();
    virtual ~neutralKaonEventDisplay();

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

private:
    // K0 candidate data (read from tree during visualization)
    static const Int_t kMaxK0 = 10;
    Int_t _nK0Candidates;
    Int_t _k0_daughter1ID[kMaxK0];
    Int_t _k0_daughter2ID[kMaxK0];
    Int_t _k0_parentID[kMaxK0];
    Float_t _k0_creationVtxRadius[kMaxK0];
    Float_t _k0_annihilationVtxRadius[kMaxK0];
    Float_t _k0_creationVtxPos[kMaxK0][3];
    Float_t _k0_annihilationVtxPos[kMaxK0][3];
    Float_t _k0_startPos[kMaxK0][3];
    Float_t _k0_endPos[kMaxK0][3];

    // Variable indices for OutputManager (extend pdEventDisplay enum)
    enum enumNK0EventDisplayVars {
        ednK0Candidates = pdEventDisplay::edmaxvars,
        edk0_daughter1ID,
        edk0_daughter2ID,
        edk0_parentID,
        edk0_creationVtxRadius,
        edk0_annihilationVtxRadius,
        edk0_creationVtxPos,
        edk0_annihilationVtxPos,
        edk0_startPos,
        edk0_endPos,
        nk0maxvars
    };
};

#endif

