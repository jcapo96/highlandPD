// File: PreliminaryK0SelectionStatistic.C
// Statistical analysis of preliminary K0 selection criteria
// Applies the same selection criteria and checks for K0 in truth information
// Selection criteria:
// - Two daughters of the beam particle whose origin is 10-50cm away from the end of the beam particle
// - These two daughters start positions are closer than 1 cm
// - Then checks if there is a K0 in truth that is a daughter of the beam particle

void PreliminaryK0SelectionStatistic(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
    // Open the file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // Get the MiniTree
    TTree *MiniTree = (TTree*)f->Get("MiniTree");
    if (!MiniTree) {
        std::cerr << "Error getting MiniTree from file" << std::endl;
        f->Close();
        return;
    }

    // Set up branch addresses
    AnaSpillB *spill = nullptr;
    MiniTree->SetBranchAddress("Spill", &spill);

    std::cout << "=== Preliminary K0 Selection Statistics ===" << std::endl;
    std::cout << "File: " << filename << std::endl;
    std::cout << "Selection criteria:" << std::endl;
    std::cout << "  - Two daughters of beam particle whose origin is 10-50cm away from beam particle end" << std::endl;
    std::cout << "  - These two daughters start positions are closer than 1cm" << std::endl;
    std::cout << "  - Check for K0 in truth that is a daughter of beam particle" << std::endl;
    std::cout << std::endl;

    Long64_t nEntries = MiniTree->GetEntries();
    std::cout << "Total events in file: " << nEntries << std::endl;

    // Statistics counters
    int totalEvents = 0;
    int eventsWithBeamParticle = 0;
    int eventsWithEnoughDaughters = 0;
    int eventsPassingDistanceCriteria = 0;
    int eventsPassingAllCriteria = 0;
    int eventsWithK0InTruth = 0;
    int eventsPassingCriteriaWithK0 = 0;

    // Detailed statistics
    std::vector<int> eventsWithK0Truth;
    std::vector<int> eventsPassingWithK0Truth;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        MiniTree->GetEntry(entry);
        totalEvents++;

        if (entry % 10000 == 0) {
            std::cout << "Processing event " << entry << " / " << nEntries << std::endl;
        }

        // Check if we have a spill and beam particle
        if (!spill || !spill->Beam) {
            continue;
        }

        AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
        if (!beam || !beam->BeamParticle) {
            continue;
        }

        AnaParticleB* beamParticle = beam->BeamParticle;
        AnaTrueParticleB* beamTrue = beamParticle->GetTrueParticle();
        if (!beamTrue) {
            continue;
        }

        eventsWithBeamParticle++;

        // Check for K0 in truth information that is a daughter of beam particle
        bool hasK0Daughter = false;
        for (UInt_t i = 0; i < spill->TrueParticles.size(); i++) {
            AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(spill->TrueParticles[i]);
            if (truePart && truePart->ParentID == beamTrue->ID && truePart->PDG == 310) {
                hasK0Daughter = true;
                break;
            }
        }

        if (hasK0Daughter) {
            eventsWithK0InTruth++;
            eventsWithK0Truth.push_back(entry);
        }

        // Find daughters of the beam particle in RECONSTRUCTED particles
        std::vector<AnaParticlePD*> beamDaughters;

        for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
            AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
            if (!bunch) continue;

            for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                AnaParticleB* particle = bunch->Particles[j];
                if (!particle) continue;

                AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                if (!particlePD) continue;

                AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                if (!truePart) continue;

                // Check if this reconstructed particle is a daughter of the beam particle
                if (truePart->ParentID == beamTrue->ID) {
                    beamDaughters.push_back(particlePD);
                }
            }
        }

        if (beamDaughters.size() < 2) {
            continue;
        }

        eventsWithEnoughDaughters++;

        // Get beam particle end position
        float beamEndX = beamTrue->PositionEnd[0];
        float beamEndY = beamTrue->PositionEnd[1];
        float beamEndZ = beamTrue->PositionEnd[2];

        // Check all pairs of daughters
        bool foundValidPair = false;
        for (size_t i = 0; i < beamDaughters.size(); i++) {
            if (foundValidPair) break;

            for (size_t j = i + 1; j < beamDaughters.size(); j++) {
                AnaParticlePD* dau1 = beamDaughters[i];
                AnaParticlePD* dau2 = beamDaughters[j];

                // Check for valid start positions
                if (dau1->PositionStart[0] == -999 || dau1->PositionStart[1] == -999 || dau1->PositionStart[2] == -999 ||
                    dau2->PositionStart[0] == -999 || dau2->PositionStart[1] == -999 || dau2->PositionStart[2] == -999) {
                    continue;
                }

                // Calculate distance from beam particle end to each daughter start
                float dist1 = sqrt((dau1->PositionStart[0] - beamEndX) * (dau1->PositionStart[0] - beamEndX) +
                                 (dau1->PositionStart[1] - beamEndY) * (dau1->PositionStart[1] - beamEndY) +
                                 (dau1->PositionStart[2] - beamEndZ) * (dau1->PositionStart[2] - beamEndZ));

                float dist2 = sqrt((dau2->PositionStart[0] - beamEndX) * (dau2->PositionStart[0] - beamEndX) +
                                 (dau2->PositionStart[1] - beamEndY) * (dau2->PositionStart[1] - beamEndY) +
                                 (dau2->PositionStart[2] - beamEndZ) * (dau2->PositionStart[2] - beamEndZ));

                // Check if both daughters are reasonably away from beam particle end
                // Tightened criteria: at least 10cm away and at most 50cm away
                if (dist1 >= 3 && dist2 >= 3) {
                    eventsPassingDistanceCriteria++;

                    // Calculate distance between the two daughters' start positions
                    float distBetweenDaughters = sqrt((dau1->PositionStart[0] - dau2->PositionStart[0]) * (dau1->PositionStart[0] - dau2->PositionStart[0]) +
                                                     (dau1->PositionStart[1] - dau2->PositionStart[1]) * (dau1->PositionStart[1] - dau2->PositionStart[1]) +
                                                     (dau1->PositionStart[2] - dau2->PositionStart[2]) * (dau1->PositionStart[2] - dau2->PositionStart[2]));

                    // Check if daughters start positions are close (tightened to 1cm)
                    if (distBetweenDaughters < 3.0) {
                        eventsPassingAllCriteria++;
                        foundValidPair = true;

                        // Check if this event also has K0 in truth
                        if (hasK0Daughter) {
                            eventsPassingCriteriaWithK0++;
                            eventsPassingWithK0Truth.push_back(entry);
                        }
                        break;
                    }
                }
            }
        }
    }

    // Print results
    std::cout << std::endl;
    std::cout << "=== RESULTS ===" << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << "Events with beam particle: " << eventsWithBeamParticle << " ("
              << (100.0 * eventsWithBeamParticle / totalEvents) << "%)" << std::endl;
    std::cout << "Events with K0 in truth (daughter of beam): " << eventsWithK0InTruth << " ("
              << (100.0 * eventsWithK0InTruth / totalEvents) << "%)" << std::endl;
    std::cout << "Events with enough daughters (>=2): " << eventsWithEnoughDaughters << " ("
              << (100.0 * eventsWithEnoughDaughters / totalEvents) << "%)" << std::endl;
    std::cout << "Events passing distance criteria (3cm from beam end): " << eventsPassingDistanceCriteria << " ("
              << (100.0 * eventsPassingDistanceCriteria / totalEvents) << "%)" << std::endl;
    std::cout << "Events passing ALL selection criteria: " << eventsPassingAllCriteria << " ("
              << (100.0 * eventsPassingAllCriteria / totalEvents) << "%)" << std::endl;
    std::cout << "Events passing criteria AND have K0 in truth: " << eventsPassingCriteriaWithK0 << " ("
              << (100.0 * eventsPassingCriteriaWithK0 / totalEvents) << "%)" << std::endl;

    std::cout << std::endl;
    std::cout << "=== EFFICIENCY ANALYSIS ===" << std::endl;
    if (eventsWithK0InTruth > 0) {
        std::cout << "K0 efficiency (events with K0 that pass criteria): "
                  << (100.0 * eventsPassingCriteriaWithK0 / eventsWithK0InTruth) << "%" << std::endl;
    }
    if (eventsPassingAllCriteria > 0) {
        std::cout << "Purity (events passing criteria that have K0): "
                  << (100.0 * eventsPassingCriteriaWithK0 / eventsPassingAllCriteria) << "%" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "=== EVENT LISTS ===" << std::endl;
    std::cout << "Events with K0 in truth (total: " << eventsWithK0Truth.size() << "): ";
    for (size_t i = 0; i < std::min(eventsWithK0Truth.size(), (size_t)10); i++) {
        std::cout << eventsWithK0Truth[i];
        if (i < std::min(eventsWithK0Truth.size(), (size_t)10) - 1) std::cout << ", ";
    }
    if (eventsWithK0Truth.size() > 10) std::cout << "...";
    std::cout << std::endl;

    std::cout << "Events passing criteria with K0 (total: " << eventsPassingWithK0Truth.size() << "): ";
    for (size_t i = 0; i < std::min(eventsPassingWithK0Truth.size(), (size_t)10); i++) {
        std::cout << eventsPassingWithK0Truth[i];
        if (i < std::min(eventsPassingWithK0Truth.size(), (size_t)10) - 1) std::cout << ", ";
    }
    if (eventsPassingWithK0Truth.size() > 10) std::cout << "...";
    std::cout << std::endl;

    f->Close();
    std::cout << std::endl;
    std::cout << "Analysis complete!" << std::endl;
}
