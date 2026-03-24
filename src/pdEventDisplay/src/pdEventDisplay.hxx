#ifndef pdEventDisplay_h
#define pdEventDisplay_h

#include "EventDisplayBase.hxx"
#include "pdDataClasses.hxx"
#include <string>
#include <vector>
#include <map>

// Forward declarations
class ToyBoxNeutralKaon;
class TEveScene;
class TEveProjectionManager;
class OutputManager;
class AnaEventB;

/// ProtoDUNE-specific implementation of EventDisplayBase
/// Handles ProtoDUNE detector geometry and K0 candidate visualization
class pdEventDisplay : public EventDisplayBase {
public:
    pdEventDisplay();
    virtual ~pdEventDisplay();

protected:
    // ========== Pure virtual methods from EventDisplayBase ==========

    /// Add ProtoDUNE-specific variables to EventDisplayData tree
    /// Includes: particles, hits, K0 candidates
    virtual void AddExperimentVariables(OutputManager& output) override;

    /// Fill ProtoDUNE-specific event data
    /// Fills particle positions, PDG codes, hits, K0 candidates
    virtual void FillExperimentData(OutputManager& output, const AnaEventB& event, void* box) override;

    /// Draw ProtoDUNE-SP detector geometry in 3D
    /// Draws APAs, CPAs, and field cage
    virtual void DrawDetectorGeometry(TEveScene* scene) override;

    /// Draw ProtoDUNE-SP detector geometry for 2D projections
    virtual void DrawDetectorGeometry2D(TEveProjectionManager* manager, const std::string& projection_type) override;

    /// Get color for particle based on PDG code (ProtoDUNE color scheme)
    virtual int GetParticleColor(int pdg) override;

    // ========== Override virtual methods for custom behavior ==========

    /// Read ProtoDUNE-specific event data from tree
    virtual bool ReadEventData(TTree* tree, Int_t run, Int_t subrun, Int_t event) override;

    /// Draw particles and hits in 3D (ProtoDUNE-specific)
    virtual void DrawParticles3D(TEveScene* scene) override;

    /// Draw particles and hits in 2D projections (ProtoDUNE-specific)
    virtual void DrawParticles2D(TEveProjectionManager* manager, const std::string& projection_type) override;

    // ========== Analysis-Specific Virtual Methods ==========

    /// Add analysis-specific variables to EventDisplayData tree
    /// Called by AddExperimentVariables() after adding standard ProtoDUNE variables
    /// Override in analysis-specific derived classes (e.g., neutralKaonEventDisplay)
    virtual void AddAnalysisVariables(OutputManager& output, Int_t tree_index);

    /// Fill analysis-specific data
    /// Called by FillExperimentData() after filling standard ProtoDUNE data
    /// Override in analysis-specific derived classes
    virtual void FillAnalysisData(OutputManager& output, const AnaEventB& event, void* box);

    /// Read analysis-specific data from tree
    /// Called by ReadEventData() after reading standard ProtoDUNE data
    /// Override in analysis-specific derived classes
    virtual bool ReadAnalysisData(TTree* tree);

    /// Draw analysis-specific reconstruction products in 3D
    /// Called by DrawParticles3D() after drawing particles and hits
    /// Override to draw vertices, secondary particles, etc.
    virtual void DrawAnalysisContent3D(TEveScene* scene);

    /// Draw analysis-specific reconstruction products in 2D
    /// Called by DrawParticles2D() after drawing particles
    /// Override to draw vertices, secondary particles, etc.
    virtual void DrawAnalysisContent2D(TEveProjectionManager* manager, const std::string& projection_type);

protected:
    // ========== Protected Data Members for Analysis-Specific Variables ==========

    /// Enum counter for derived classes to extend variable indices
    enum enumPDEventDisplayVars {
        edmaxvars = 200  // Reserve space for ProtoDUNE variables
    };

protected:
    // ========== Helper Methods ==========

    /// Initialize particle color map
    void InitializeParticleColors();

    /// Get particle name from PDG code
    std::string GetParticleName(int pdg);

    // ========== Data Members ==========

    /// Particle color map (PDG code -> ROOT color)
    std::map<int, int> _particleColors;

    // Data read from tree (stored for visualization)
    Int_t _nParticles;
    Int_t _totalHits;
    Int_t _nK0Candidates;
    Int_t _hasTrueK0;

    // Particle data arrays
    static const Int_t kMaxParticles = 200;
    Int_t _particle_uniqueID[kMaxParticles];
    Int_t _particle_PDG[kMaxParticles];
    Float_t _particle_startPos[kMaxParticles][3];
    Float_t _particle_endPos[kMaxParticles][3];
    Float_t _particle_charge[kMaxParticles];
    Float_t _particle_momentum[kMaxParticles];
    Int_t _particle_nHits[kMaxParticles];
    Int_t _particle_hitIndex[kMaxParticles];

    // Hit data arrays
    static const Int_t kMaxHits = 50000;
    Float_t* _hit_X;
    Float_t* _hit_Y;
    Float_t* _hit_Z;

    // K0 candidate data arrays
    static const Int_t kMaxK0 = 10;
    Float_t _k0_creationVtxPos[kMaxK0][3];
    Float_t _k0_annihilationVtxPos[kMaxK0][3];
    Float_t _k0_startPos[kMaxK0][3];
    Float_t _k0_endPos[kMaxK0][3];
    Int_t _k0_daughter1ID[kMaxK0];
    Int_t _k0_daughter2ID[kMaxK0];
    Int_t _k0_parentID[kMaxK0];
    Float_t _k0_creationVtxRadius[kMaxK0];
    Float_t _k0_annihilationVtxRadius[kMaxK0];

    // Variable indices for EventDisplayData tree
    // Start after EventDisplayBase variables
    enum enumEventDisplayVar {
        edrun = EventDisplayBase::edbasemaxvars,
        edsubrun,
        edevent,
        edhasTrueK0,

        // Particle data
        ednParticles,
        edparticle_startPos,
        edparticle_endPos,
        edparticle_PDG,
        edparticle_momentum,
        edparticle_uniqueID,
        edparticle_charge,

        // Particle hit data
        edparticle_nHits,
        edparticle_hitIndex,
        edtotalHits,
        edhit_X,
        edhit_Y,
        edhit_Z,

        // K0 candidate data
        ednK0Candidates,
        edk0_creationVtxPos,
        edk0_annihilationVtxPos,
        edk0_startPos,
        edk0_endPos,
        edk0_creationVtxRadius,
        edk0_annihilationVtxRadius,
        edk0_daughter1ID,
        edk0_daughter2ID,
        edk0_parentID,

        NEventDisplayVars
    };
};

#endif
