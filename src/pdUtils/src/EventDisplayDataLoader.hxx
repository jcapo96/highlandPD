#ifndef EventDisplayDataLoader_h
#define EventDisplayDataLoader_h

#include <TTree.h>
#include <vector>

// Forward declarations
class AnaEventB;
class ToyBoxNeutralKaon;
class AnaParticlePD;
class AnaNeutralParticlePD;

/// Helper class to load event display data from the EventDisplayData tree
/// and reconstruct objects for visualization
class EventDisplayDataLoader {
public:
    EventDisplayDataLoader(TTree* tree);
    virtual ~EventDisplayDataLoader();

    /// Find and load data for a specific event
    bool LoadEvent(Int_t run, Int_t subrun, Int_t event);

    /// Get reconstructed event object (valid after LoadEvent)
    AnaEventB* GetEvent() { return _event; }

    /// Get reconstructed ToyBox object (valid after LoadEvent)
    ToyBoxNeutralKaon* GetToyBox() { return _toyBox; }

private:
    /// Build index on the tree for fast lookup
    void BuildIndex();

    /// Reconstruct event and toybox from loaded data
    void ReconstructObjects();

    /// Clean up allocated memory
    void ClearObjects();

private:
    TTree* _tree;
    bool _indexBuilt;

    // Reconstructed objects
    AnaEventB* _event;
    ToyBoxNeutralKaon* _toyBox;
    std::vector<AnaParticlePD*> _particles;
    std::vector<AnaNeutralParticlePD*> _neutralParticles;

    // Data buffers for reading from tree
    Int_t _run, _subrun, _evt, _hasTrueK0;
    Int_t _nParticles, _totalHits, _nK0Candidates;

    // Arrays for particle data (max 200 particles)
    Int_t _particle_uniqueID[200];
    Int_t _particle_PDG[200];
    Int_t _particle_charge[200];
    Float_t _particle_momentum[200];
    Int_t _particle_nHits[200];
    Int_t _particle_hitIndex[200];
    Float_t _particle_startPos[200][3];
    Float_t _particle_endPos[200][3];
    Float_t _particle_startDir[200][3];
    Float_t _particle_endDir[200][3];

    // Arrays for hit data (max 50000 hits)
    Float_t* _hit_X;
    Float_t* _hit_Y;
    Float_t* _hit_Z;

    // Arrays for K0 candidate data (max 10 K0s)
    Int_t _k0_daughter1ID[10];
    Int_t _k0_daughter2ID[10];
    Int_t _k0_parentID[10];
    Float_t _k0_creationVtxRadius[10];
    Float_t _k0_creationVtxPos[10][3];
    Float_t _k0_annihilationVtxPos[10][3];
    Float_t _k0_startPos[10][3];
    Float_t _k0_endPos[10][3];
    Float_t _k0_startDir[10][3];
    Float_t _k0_endDir[10][3];
};

#endif

