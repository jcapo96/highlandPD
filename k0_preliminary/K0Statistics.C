// File: K0Statistics.C
// Analyzes all events to find K0 particles that are daughters of K+ beam particles
// Based on truth information
void K0Statistics(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
    // Open the file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // Get the tree
    TTree *MiniTree = (TTree*)f->Get("MiniTree");
    if (!MiniTree) {
        std::cerr << "Could not find MiniTree in file" << std::endl;
        return;
    }

    // Set up the object for the branch
    AnaSpill* spill = nullptr;
    MiniTree->SetBranchAddress("Spill", &spill);

    Long64_t nEntries = MiniTree->GetEntries();
    std::cout << "Total events to process: " << nEntries << std::endl;

    // Statistics counters
    int totalEvents = 0;
    int eventsWithKPlusBeam = 0;
    int eventsWithK0FromKPlusBeam = 0;
    int totalK0FromKPlusBeam = 0;
    int eventsWithK0FromKPlusNonBeam = 0;
    int totalK0FromKPlusNonBeam = 0;
    int totalK0WithPiDaughters = 0;
    int eventsWithK0PiDaughters = 0;
    int totalK0FromAnyKPlus = 0;
    int eventsWithK0FromAnyKPlus = 0;

    // Detailed statistics
    std::map<int, int> k0DecayModes; // PDG -> count
    std::vector<int> k0DecayChains; // Store decay chain information

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);
        if (!spill) continue;

        totalEvents++;

        // Check if there's a beam particle
        AnaParticleB* beamParticle = nullptr;
        if (spill->Beam) {
            AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
            if (beam && beam->BeamParticle) {
                beamParticle = beam->BeamParticle;
            }
        }

        // Check if beam particle is K+ (PDG = 321)
        bool hasKPlusBeam = false;
        AnaTrueParticleB* beamTrue = nullptr;
        if (beamParticle) {
            beamTrue = beamParticle->GetTrueParticle();
            if (beamTrue && beamTrue->PDG == 321) {
                hasKPlusBeam = true;
                eventsWithKPlusBeam++;
            }
        }

        UInt_t nTrueParticles = spill->TrueParticles.size();

        // Find all K+ particles in the event (both beam and non-beam)
        std::vector<AnaTrueParticleB*> kPlusParticles;
        for (UInt_t j = 0; j < nTrueParticles; j++) {
            AnaTrueParticleB* truePart = spill->TrueParticles[j];
            if (!truePart) continue;
            if (truePart->PDG == 321) { // K+
                kPlusParticles.push_back(truePart);
            }
        }

        // Find all K0 particles that are daughters of any K+ particle
        std::vector<AnaTrueParticleB*> k0FromKPlusBeam;
        std::vector<AnaTrueParticleB*> k0FromKPlusNonBeam;
        std::vector<AnaTrueParticleB*> k0FromAnyKPlus;

        for (UInt_t j = 0; j < nTrueParticles; j++) {
            AnaTrueParticleB* truePart = spill->TrueParticles[j];
            if (!truePart) continue;

            // Check if this is a K0 (PDG = 310)
            if (truePart->PDG == 310) {
                // Check if this K0 is a daughter of any K+ particle
                for (UInt_t k = 0; k < kPlusParticles.size(); k++) {
                    AnaTrueParticleB* kPlus = kPlusParticles[k];
                    if (truePart->ParentID == kPlus->ID) {
                        k0FromAnyKPlus.push_back(truePart);
                        totalK0FromAnyKPlus++;

                        // Check if this K+ is the beam particle
                        if (beamTrue && kPlus->ID == beamTrue->ID) {
                            k0FromKPlusBeam.push_back(truePart);
                            totalK0FromKPlusBeam++;
                        } else {
                            k0FromKPlusNonBeam.push_back(truePart);
                            totalK0FromKPlusNonBeam++;
                        }
                        break; // Found the parent, no need to check other K+ particles
                    }
                }
            }
        }

        // Update event counters
        if (k0FromKPlusBeam.size() > 0) {
            eventsWithK0FromKPlusBeam++;
        }
        if (k0FromKPlusNonBeam.size() > 0) {
            eventsWithK0FromKPlusNonBeam++;
        }
        if (k0FromAnyKPlus.size() > 0) {
            eventsWithK0FromAnyKPlus++;
        }

        // Analyze decays of all K0 particles from any K+ parent
        for (UInt_t k0Idx = 0; k0Idx < k0FromAnyKPlus.size(); k0Idx++) {
            AnaTrueParticleB* k0Particle = k0FromAnyKPlus[k0Idx];

            // Find daughters of this K0
            std::vector<AnaTrueParticleB*> k0Daughters;
            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Check if this particle is a daughter of the K0
                if (truePart->ParentID == k0Particle->ID) {
                    k0Daughters.push_back(truePart);
                }
            }

            // Analyze decay mode
            if (k0Daughters.size() > 0) {
                // Count different daughter types
                int piPlusCount = 0;
                int piMinusCount = 0;
                int piZeroCount = 0;
                int electronCount = 0;
                int muonCount = 0;
                int otherCount = 0;

                for (UInt_t dIdx = 0; dIdx < k0Daughters.size(); dIdx++) {
                    AnaTrueParticleB* daughter = k0Daughters[dIdx];
                    int pdg = daughter->PDG;

                    if (pdg == 211) piPlusCount++;
                    else if (pdg == -211) piMinusCount++;
                    else if (pdg == 111) piZeroCount++;
                    else if (pdg == 11 || pdg == -11) electronCount++;
                    else if (pdg == 13 || pdg == -13) muonCount++;
                    else otherCount++;
                }

                // Determine decay mode
                int decayMode = 0;
                if (piPlusCount == 1 && piMinusCount == 1) {
                    decayMode = 1; // π+π-
                } else if (piZeroCount == 2) {
                    decayMode = 2; // π0π0
                } else if (piZeroCount == 1 && electronCount == 1) {
                    decayMode = 3; // π0e+νe
                } else if (piZeroCount == 1 && muonCount == 1) {
                    decayMode = 4; // π0μ+νμ
                } else if (electronCount == 2) {
                    decayMode = 5; // e+e-
                } else if (muonCount == 2) {
                    decayMode = 6; // μ+μ-
                } else {
                    decayMode = 99; // Other
                }

                k0DecayModes[decayMode]++;

                // Store decay chain information
                if (decayMode == 1) { // π+π- decay
                    totalK0WithPiDaughters++;
                    if (eventsWithK0PiDaughters == 0 ||
                        std::find(k0DecayChains.begin(), k0DecayChains.end(), i) == k0DecayChains.end()) {
                        eventsWithK0PiDaughters++;
                        k0DecayChains.push_back(i);
                    }
                }

                // Print detailed decay information for first few events
                if (i < 10) {
                    std::cout << "Event " << i << ": K+ -> K0(ID:" << k0Particle->ID << ") -> ";
                    for (UInt_t dIdx = 0; dIdx < k0Daughters.size(); dIdx++) {
                        AnaTrueParticleB* daughter = k0Daughters[dIdx];
                        if (dIdx > 0) std::cout << " + ";
                        std::cout << "PDG:" << daughter->PDG;
                    }
                    std::cout << " (Mode: " << decayMode << ")" << std::endl;
                }
            }
        }

        // Progress indicator
        if (i % 5000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries
                      << " (K+ beam events: " << eventsWithKPlusBeam
                      << ", K0 from K+ beam: " << eventsWithK0FromKPlusBeam
                      << ", K0 from K+ non-beam: " << eventsWithK0FromKPlusNonBeam << ")" << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== K0 Statistics Analysis ===" << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << "Events with K+ beam: " << eventsWithKPlusBeam << std::endl;
    std::cout << "Events with K0 daughters of K+ beam: " << eventsWithK0FromKPlusBeam << std::endl;
    std::cout << "Events with K0 daughters of K+ non-beam: " << eventsWithK0FromKPlusNonBeam << std::endl;
    std::cout << "Events with K0 daughters of any K+: " << eventsWithK0FromAnyKPlus << std::endl;
    std::cout << "Total K0 particles from K+ beam: " << totalK0FromKPlusBeam << std::endl;
    std::cout << "Total K0 particles from K+ non-beam: " << totalK0FromKPlusNonBeam << std::endl;
    std::cout << "Total K0 particles from any K+: " << totalK0FromAnyKPlus << std::endl;
    std::cout << "K0 particles with π+π- daughters: " << totalK0WithPiDaughters << std::endl;
    std::cout << "Events with K0→π+π- decays: " << eventsWithK0PiDaughters << std::endl;

    std::cout << "\n=== K0 Decay Modes ===" << std::endl;
    std::cout << "Mode 1 (π+π-): " << k0DecayModes[1] << std::endl;
    std::cout << "Mode 2 (π0π0): " << k0DecayModes[2] << std::endl;
    std::cout << "Mode 3 (π0e+νe): " << k0DecayModes[3] << std::endl;
    std::cout << "Mode 4 (π0μ+νμ): " << k0DecayModes[4] << std::endl;
    std::cout << "Mode 5 (e+e-): " << k0DecayModes[5] << std::endl;
    std::cout << "Mode 6 (μ+μ-): " << k0DecayModes[6] << std::endl;
    std::cout << "Mode 99 (Other): " << k0DecayModes[99] << std::endl;

    // Calculate fractions
    if (totalK0FromAnyKPlus > 0) {
        std::cout << "\n=== Decay Fractions (All K0 from K+ parents) ===" << std::endl;
        std::cout << "π+π- fraction: " << (double)k0DecayModes[1] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "π0π0 fraction: " << (double)k0DecayModes[2] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "π0e+νe fraction: " << (double)k0DecayModes[3] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "π0μ+νμ fraction: " << (double)k0DecayModes[4] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "e+e- fraction: " << (double)k0DecayModes[5] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "μ+μ- fraction: " << (double)k0DecayModes[6] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
        std::cout << "Other fraction: " << (double)k0DecayModes[99] / totalK0FromAnyKPlus * 100.0 << "%" << std::endl;
    }

    if (eventsWithKPlusBeam > 0) {
        std::cout << "\n=== Beam Interaction Fractions ===" << std::endl;
        std::cout << "K+ beam events with K0 production: " << (double)eventsWithK0FromKPlusBeam / eventsWithKPlusBeam * 100.0 << "%" << std::endl;
        std::cout << "K+ beam events with K0→π+π-: " << (double)eventsWithK0PiDaughters / eventsWithKPlusBeam * 100.0 << "%" << std::endl;
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Found " << totalK0FromKPlusBeam << " K0 particles from K+ beam interactions" << std::endl;
    std::cout << "Found " << totalK0FromKPlusNonBeam << " K0 particles from K+ non-beam interactions" << std::endl;
    std::cout << "Total K0 particles from any K+ parent: " << totalK0FromAnyKPlus << std::endl;
    std::cout << "Of all K0 from K+ parents, " << totalK0WithPiDaughters << " decayed to π+π- ("
              << (totalK0FromAnyKPlus > 0 ? (double)totalK0WithPiDaughters / totalK0FromAnyKPlus * 100.0 : 0.0) << "%)" << std::endl;

    // Close the file
    f->Close();
    delete f;
}
