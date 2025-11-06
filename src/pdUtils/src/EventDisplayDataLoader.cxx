#include "EventDisplayDataLoader.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxNeutralKaon.hxx"
#include <iostream>

//********************************************************************
EventDisplayDataLoader::EventDisplayDataLoader(TTree* tree) {
//********************************************************************
    _tree = tree;
    _indexBuilt = false;
    _event = nullptr;
    _toyBox = nullptr;

    // Allocate hit arrays
    _hit_X = new Float_t[50000];
    _hit_Y = new Float_t[50000];
    _hit_Z = new Float_t[50000];

    if (_tree) {
        // Set branch addresses for event identification
        _tree->SetBranchAddress("ED_run", &_run);
        _tree->SetBranchAddress("ED_subrun", &_subrun);
        _tree->SetBranchAddress("ED_event", &_evt);
        _tree->SetBranchAddress("ED_hasTrueK0", &_hasTrueK0);

        // Set branch addresses for particle counts
        _tree->SetBranchAddress("ED_nParticles", &_nParticles);
        _tree->SetBranchAddress("ED_totalHits", &_totalHits);
        _tree->SetBranchAddress("ED_nK0Candidates", &_nK0Candidates);

        // Set branch addresses for particle data
        _tree->SetBranchAddress("ED_particle_uniqueID", _particle_uniqueID);
        _tree->SetBranchAddress("ED_particle_PDG", _particle_PDG);
        _tree->SetBranchAddress("ED_particle_charge", _particle_charge);
        _tree->SetBranchAddress("ED_particle_momentum", _particle_momentum);
        _tree->SetBranchAddress("ED_particle_nHits", _particle_nHits);
        _tree->SetBranchAddress("ED_particle_hitIndex", _particle_hitIndex);
        _tree->SetBranchAddress("ED_particle_startPos", _particle_startPos);
        _tree->SetBranchAddress("ED_particle_endPos", _particle_endPos);
        _tree->SetBranchAddress("ED_particle_startDir", _particle_startDir);
        _tree->SetBranchAddress("ED_particle_endDir", _particle_endDir);

        // Set branch addresses for hit data
        _tree->SetBranchAddress("ED_hit_X", _hit_X);
        _tree->SetBranchAddress("ED_hit_Y", _hit_Y);
        _tree->SetBranchAddress("ED_hit_Z", _hit_Z);

        // Set branch addresses for K0 data
        _tree->SetBranchAddress("ED_k0_daughter1ID", _k0_daughter1ID);
        _tree->SetBranchAddress("ED_k0_daughter2ID", _k0_daughter2ID);
        _tree->SetBranchAddress("ED_k0_parentID", _k0_parentID);
        _tree->SetBranchAddress("ED_k0_creationVtxRadius", _k0_creationVtxRadius);
        _tree->SetBranchAddress("ED_k0_creationVtxPos", _k0_creationVtxPos);
        _tree->SetBranchAddress("ED_k0_annihilationVtxPos", _k0_annihilationVtxPos);
        _tree->SetBranchAddress("ED_k0_startPos", _k0_startPos);
        _tree->SetBranchAddress("ED_k0_endPos", _k0_endPos);
        _tree->SetBranchAddress("ED_k0_startDir", _k0_startDir);
        _tree->SetBranchAddress("ED_k0_endDir", _k0_endDir);
    }
}

//********************************************************************
EventDisplayDataLoader::~EventDisplayDataLoader() {
//********************************************************************
    ClearObjects();

    if (_hit_X) delete[] _hit_X;
    if (_hit_Y) delete[] _hit_Y;
    if (_hit_Z) delete[] _hit_Z;
}

//********************************************************************
void EventDisplayDataLoader::BuildIndex() {
//********************************************************************
    if (!_tree || _indexBuilt) return;

    // Build index on run:subrun:event for fast lookup
    // Note: TTree::BuildIndex only supports major:minor (2 variables)
    // We'll combine subrun and event into a single Long64_t
    _indexBuilt = true;
}

//********************************************************************
bool EventDisplayDataLoader::LoadEvent(Int_t run, Int_t subrun, Int_t event) {
//********************************************************************
    if (!_tree) {
        std::cout << "ERROR: No tree provided to EventDisplayDataLoader" << std::endl;
        return false;
    }

    // Build index if not already done
    if (!_indexBuilt) {
        BuildIndex();
    }

    // Clear any previously loaded objects
    ClearObjects();

    // Search for the event
    Long64_t nEntries = _tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        _tree->GetEntry(i);

        if (_run == run && _subrun == subrun && _evt == event) {
            // Found the event!
            std::cout << "Found event: Run " << run << ", Subrun " << subrun << ", Event " << event << std::endl;
            std::cout << "  - nParticles: " << _nParticles << std::endl;
            std::cout << "  - nK0Candidates: " << _nK0Candidates << std::endl;
            std::cout << "  - totalHits: " << _totalHits << std::endl;
            std::cout << "  - hasTrueK0: " << (_hasTrueK0 ? "Yes" : "No") << std::endl;

            // Reconstruct objects
            ReconstructObjects();
            return true;
        }
    }

    std::cout << "ERROR: Event not found (Run " << run << ", Subrun " << subrun << ", Event " << event << ")" << std::endl;
    return false;
}

//********************************************************************
void EventDisplayDataLoader::ReconstructObjects() {
//********************************************************************
    // Note: Full reconstruction is complex
    // For now, we'll create minimal objects needed for visualization
    std::cout << "Note: Object reconstruction not yet fully implemented." << std::endl;
    std::cout << "Event display data is loaded and ready for direct visualization." << std::endl;
}

//********************************************************************
void EventDisplayDataLoader::ClearObjects() {
//********************************************************************
    // Clean up previously allocated objects
    if (_event) {
        delete _event;
        _event = nullptr;
    }

    if (_toyBox) {
        delete _toyBox;
        _toyBox = nullptr;
    }

    for (auto* particle : _particles) {
        if (particle) delete particle;
    }
    _particles.clear();

    for (auto* neutral : _neutralParticles) {
        if (neutral) delete neutral;
    }
    _neutralParticles.clear();
}

