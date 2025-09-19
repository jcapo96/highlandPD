// File: PreliminaryK0Selection.C
// Modified selection criteria: loop over beam particle daughters
// Keep events where there are two daughters of the beam particle whose origin is 10-50cm away
// from the end of the beam particle and that these two daughters start positions are closer than 5 cm
void PreliminaryK0Selection(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
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

    // Define particle types and their colors
    std::map<int, std::string> particleNames;
    std::map<int, int> particleColors;

    // Pions
    particleNames[211] = "#pi^{+}";
    particleColors[211] = kBlue;
    particleNames[-211] = "#pi^{-}";
    particleColors[-211] = kOrange;
    particleNames[111] = "#pi^{0}";
    particleColors[111] = kCyan;

    // Kaons
    particleNames[321] = "K^{+}";
    particleColors[321] = kRed;
    particleNames[-321] = "K^{-}";
    particleColors[-321] = kTeal;
    particleNames[311] = "K^{0}";
    particleColors[311] = kYellow+2;
    particleNames[-311] = "#bar{K}^{0}";
    particleColors[-311] = kTeal+2;
    particleNames[310] = "K^{0}_{S}";
    particleColors[310] = kMagenta;

    // Leptons
    particleNames[13] = "#mu^{-}";
    particleColors[13] = kBlack;
    particleNames[-13] = "#mu^{+}";
    particleColors[-13] = kViolet+2;
    particleNames[11] = "e^{-}";
    particleColors[11] = kPink;
    particleNames[-11] = "e^{+}";
    particleColors[-11] = kPink+2;
    particleNames[15] = "#tau^{-}";
    particleColors[15] = kMagenta+1;
    particleNames[-15] = "#tau^{+}";
    particleColors[-15] = kMagenta+2;

    // Neutrinos
    particleNames[12] = "#nu_{e}";
    particleColors[12] = kGreen+3;
    particleNames[-12] = "#bar{#nu}_{e}";
    particleColors[-12] = kGreen+3;
    particleNames[14] = "#nu_{#mu}";
    particleColors[14] = kGreen+3;
    particleNames[-14] = "#bar{#nu}_{#mu}";
    particleColors[-14] = kGreen+3;
    particleNames[16] = "#nu_{#tau}";
    particleColors[16] = kGreen+3;
    particleNames[-16] = "#bar{#nu}_{#tau}";
    particleColors[-16] = kGreen+3;

    // Bosons
    particleNames[22] = "#gamma";
    particleColors[22] = kGreen+2;
    particleNames[23] = "Z^{0}";
    particleColors[23] = kRed+2;
    particleNames[24] = "W^{+}";
    particleColors[24] = kRed+3;
    particleNames[-24] = "W^{-}";
    particleColors[-24] = kRed+3;
    particleNames[25] = "H^{0}";
    particleColors[25] = kOrange+3;

    // Baryons
    particleNames[2212] = "p";
    particleColors[2212] = kGray;
    particleNames[-2212] = "#bar{p}";
    particleColors[-2212] = kGray+1;
    particleNames[2112] = "n";
    particleColors[2112] = kGray+2;
    particleNames[-2112] = "#bar{n}";
    particleColors[-2112] = kGray+3;

    // Lambda particles
    particleNames[3122] = "#Lambda";
    particleColors[3122] = kCyan+2;
    particleNames[-3122] = "#bar{#Lambda}";
    particleColors[-3122] = kCyan+3;

    // Sigma particles
    particleNames[3222] = "#Sigma^{+}";
    particleColors[3222] = kBlue+2;
    particleNames[-3222] = "#bar{#Sigma}^{-}";
    particleColors[-3222] = kBlue+3;
    particleNames[3212] = "#Sigma^{0}";
    particleColors[3212] = kBlue+1;
    particleNames[-3212] = "#bar{#Sigma}^{0}";
    particleColors[-3212] = kBlue+1;
    particleNames[3112] = "#Sigma^{-}";
    particleColors[3112] = kBlue+3;
    particleNames[-3112] = "#bar{#Sigma}^{+}";
    particleColors[-3112] = kBlue+2;

    // Xi particles
    particleNames[3322] = "#Xi^{0}";
    particleColors[3322] = kYellow+3;
    particleNames[-3322] = "#bar{#Xi}^{0}";
    particleColors[-3322] = kYellow+3;
    particleNames[3312] = "#Xi^{-}";
    particleColors[3312] = kYellow+4;
    particleNames[-3312] = "#bar{#Xi}^{+}";
    particleColors[-3312] = kYellow+4;

    // Omega particles
    particleNames[3334] = "#Omega^{-}";
    particleColors[3334] = kMagenta+3;
    particleNames[-3334] = "#bar{#Omega}^{+}";
    particleColors[-3334] = kMagenta+3;

    // Quarks (for completeness, though rarely seen as free particles)
    particleNames[1] = "d";
    particleColors[1] = kRed+1;
    particleNames[-1] = "#bar{d}";
    particleColors[-1] = kRed+1;
    particleNames[2] = "u";
    particleColors[2] = kBlue+1;
    particleNames[-2] = "#bar{u}";
    particleColors[-2] = kBlue+1;
    particleNames[3] = "s";
    particleColors[3] = kGreen+1;
    particleNames[-3] = "#bar{s}";
    particleColors[-3] = kGreen+1;
    particleNames[4] = "c";
    particleColors[4] = kMagenta+1;
    particleNames[-4] = "#bar{c}";
    particleColors[-4] = kMagenta+1;
    particleNames[5] = "b";
    particleColors[5] = kOrange+1;
    particleNames[-5] = "#bar{b}";
    particleColors[-5] = kOrange+1;
    particleNames[6] = "t";
    particleColors[6] = kCyan+1;
    particleNames[-6] = "#bar{t}";
    particleColors[-6] = kCyan+1;

    // Nuclear particles (nuclei)
    particleNames[1000010020] = "^{2}H";  // Deuterium
    particleColors[1000010020] = kGray+1;
    particleNames[1000010030] = "^{3}H";  // Tritium
    particleColors[1000010030] = kGray+2;
    particleNames[1000020030] = "^{3}He"; // Helium-3
    particleColors[1000020030] = kGray+3;
    particleNames[1000020040] = "^{4}He"; // Helium-4
    particleColors[1000020040] = kGray+4;
    particleNames[1000120240] = "^{24}Mg"; // Magnesium-24
    particleColors[1000120240] = kGray+5;
    particleNames[1000120270] = "^{27}Al"; // Aluminum-27
    particleColors[1000120270] = kGray+6;
    particleNames[1000140310] = "^{31}P"; // Phosphorus-31
    particleColors[1000140310] = kGray+8;
    particleNames[1000180360] = "^{36}Ar"; // Argon-36
    particleColors[1000180360] = kGray+7;
    particleNames[1000180400] = "^{40}Ar"; // Argon-40
    particleColors[1000180400] = kGray+9;
    particleNames[1000190390] = "^{39}K"; // Potassium-39
    particleColors[1000190390] = kGray+10;
    particleNames[1000180380] = "^{38}Ar"; // Argon-38
    particleColors[1000180380] = kGray+11;
    particleNames[1000140280] = "^{28}Si"; // Silicon-28
    particleColors[1000140280] = kGray+12;

    // Loop over entries and plot one event at a time
    Long64_t nEntries = MiniTree->GetEntries();
    // DEBUG: Limit to first 3 events for debugging
    // nEntries = std::min(nEntries, (Long64_t)3);
    int selectedEventCount = 0;

    // Histogram for reconstructed invariant masses
    TH1F* hRecoInvariantMass = new TH1F("hRecoInvariantMass", "Selected Events Invariant Mass;Mass [MeV/c^{2}];Events", 50, 400, 600);
    hRecoInvariantMass->SetLineColor(kBlue);
    hRecoInvariantMass->SetLineWidth(2);
    hRecoInvariantMass->SetFillColor(kBlue);
    hRecoInvariantMass->SetFillStyle(3004);

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);

        std::cout << "\n=== DEBUG: Processing Event " << i << " ===" << std::endl;

        // Get beam particle
        AnaParticleB* beamParticle = nullptr;
        AnaTrueParticleB* beamTrue = nullptr;
        if (spill->Beam) {
            AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
            if (beam && beam->BeamParticle) {
                beamParticle = beam->BeamParticle;
                beamTrue = beamParticle->GetTrueParticle();
                std::cout << "DEBUG: Found beam particle, True particle: " << (beamTrue ? "YES" : "NO") << std::endl;
                if (beamTrue) {
                    std::cout << "DEBUG: Beam True Particle - ID: " << beamTrue->ID << ", PDG: " << beamTrue->PDG << std::endl;
                    std::cout << "DEBUG: Beam True Position: (" << beamTrue->Position[0] << ", " << beamTrue->Position[1] << ", " << beamTrue->Position[2] << ")" << std::endl;
                    std::cout << "DEBUG: Beam True PositionEnd: (" << beamTrue->PositionEnd[0] << ", " << beamTrue->PositionEnd[1] << ", " << beamTrue->PositionEnd[2] << ")" << std::endl;
                }
            } else {
                std::cout << "DEBUG: No beam particle found" << std::endl;
            }
        } else {
            std::cout << "DEBUG: No beam found in spill" << std::endl;
        }

        if (!beamParticle || !beamTrue) {
            std::cout << "DEBUG: Skipping event - no beam particle or true particle" << std::endl;
            continue; // Skip events without beam particle or its true particle
        }

        // Get beam particle end position
        float beamEndX = beamTrue->PositionEnd[0];
        float beamEndY = beamTrue->PositionEnd[1];
        float beamEndZ = beamTrue->PositionEnd[2];

        // Check if beam particle end position is valid
        if (beamEndX < -900 || beamEndY < -900 || beamEndZ < -900) {
            std::cout << "DEBUG: Skipping event - invalid beam particle end position: (" << beamEndX << ", " << beamEndY << ", " << beamEndZ << ")" << std::endl;
            continue; // Skip events with invalid beam particle end position
        }

        // Find daughters of the beam particle in RECONSTRUCTED particles
        std::vector<AnaParticlePD*> beamDaughters;
        UInt_t totalRecoParticles = 0;

        std::cout << "DEBUG: Looking for reconstructed daughters of beam particle ID: " << beamTrue->ID << std::endl;

        // Loop through all bunches to find reconstructed particles
        for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
            AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
            if (!bunch) continue;

            for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                AnaParticleB* particle = bunch->Particles[j];
                if (!particle) continue;

                AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                if (!particlePD) continue;

                totalRecoParticles++;

                // Get associated true particle
                AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                if (!truePart) continue;

                // Check if this reconstructed particle is a daughter of the beam particle
                if (truePart->ParentID == beamTrue->ID) {
                    beamDaughters.push_back(particlePD);
                    std::cout << "DEBUG: Found beam daughter (reco) - TrueID: " << truePart->ID << ", PDG: " << truePart->PDG << std::endl;
                    std::cout << "DEBUG: Daughter start position: (" << particlePD->PositionStart[0] << ", " << particlePD->PositionStart[1] << ", " << particlePD->PositionStart[2] << ")" << std::endl;
                    std::cout << "DEBUG: Daughter end position: (" << particlePD->PositionEnd[0] << ", " << particlePD->PositionEnd[1] << ", " << particlePD->PositionEnd[2] << ")" << std::endl;
                }
            }
        }

        std::cout << "DEBUG: Total reconstructed particles in event: " << totalRecoParticles << std::endl;
        std::cout << "DEBUG: Total beam daughters found (reco): " << beamDaughters.size() << std::endl;

        // Need at least two daughters
        if (beamDaughters.size() < 2) {
            std::cout << "DEBUG: Skipping event - not enough daughters (" << beamDaughters.size() << ")" << std::endl;
            continue;
        }

        // Check if there are exactly two daughters whose origin is ~10cm away from beam particle end
        // and whose start positions are closer than 10cm to each other
        bool foundValidPair = false;
        AnaParticlePD* daughter1 = nullptr;
        AnaParticlePD* daughter2 = nullptr;

        std::cout << "DEBUG: Checking all daughter pairs for distance criteria..." << std::endl;
        for (UInt_t d1 = 0; d1 < beamDaughters.size(); d1++) {
            for (UInt_t d2 = d1 + 1; d2 < beamDaughters.size(); d2++) {
                AnaParticlePD* dau1 = beamDaughters[d1];
                AnaParticlePD* dau2 = beamDaughters[d2];

                if (!dau1 || !dau2) continue;

                // Get true particles for PDG info
                AnaTrueParticleB* true1 = dau1->GetTrueParticle();
                AnaTrueParticleB* true2 = dau2->GetTrueParticle();

                std::cout << "DEBUG: Checking pair - Daughter1 (TrueID:" << (true1 ? true1->ID : -1) << ", PDG:" << (true1 ? true1->PDG : -1) << "), Daughter2 (TrueID:" << (true2 ? true2->ID : -1) << ", PDG:" << (true2 ? true2->PDG : -1) << ")" << std::endl;

                // Check if daughters have valid start positions
                if (dau1->PositionStart[0] < -900 || dau1->PositionStart[1] < -900 || dau1->PositionStart[2] < -900 ||
                    dau2->PositionStart[0] < -900 || dau2->PositionStart[1] < -900 || dau2->PositionStart[2] < -900) {
                    std::cout << "DEBUG: Skipping pair - invalid start positions" << std::endl;
                    continue;
                }

                // Calculate distance from beam particle end to daughter1 start position
                float dist1 = sqrt((dau1->PositionStart[0] - beamEndX) * (dau1->PositionStart[0] - beamEndX) +
                                 (dau1->PositionStart[1] - beamEndY) * (dau1->PositionStart[1] - beamEndY) +
                                 (dau1->PositionStart[2] - beamEndZ) * (dau1->PositionStart[2] - beamEndZ));

                // Calculate distance from beam particle end to daughter2 start position
                float dist2 = sqrt((dau2->PositionStart[0] - beamEndX) * (dau2->PositionStart[0] - beamEndX) +
                                 (dau2->PositionStart[1] - beamEndY) * (dau2->PositionStart[1] - beamEndY) +
                                 (dau2->PositionStart[2] - beamEndZ) * (dau2->PositionStart[2] - beamEndZ));

                std::cout << "DEBUG: Distance from beam end - Daughter1: " << dist1 << "cm, Daughter2: " << dist2 << "cm" << std::endl;

                // Check if both daughters are reasonably away from beam particle end
                // Tightened criteria: at least 10cm away and at most 50cm away
                if (dist1 >= 3.0 && dist2 >= 3.0) {
                    std::cout << "DEBUG: Both daughters within distance range from beam end" << std::endl;

                    // Calculate distance between the two daughters' start positions
                    float distBetweenDaughters = sqrt((dau1->PositionStart[0] - dau2->PositionStart[0]) * (dau1->PositionStart[0] - dau2->PositionStart[0]) +
                                                     (dau1->PositionStart[1] - dau2->PositionStart[1]) * (dau1->PositionStart[1] - dau2->PositionStart[1]) +
                                                     (dau1->PositionStart[2] - dau2->PositionStart[2]) * (dau1->PositionStart[2] - dau2->PositionStart[2]));

                    std::cout << "DEBUG: Distance between daughters: " << distBetweenDaughters << "cm" << std::endl;

                    // Check if daughters start positions are close (tightened to 1cm)
                    if (distBetweenDaughters < 3.0) {
                        std::cout << "DEBUG: Found valid pair meeting all criteria!" << std::endl;
                        foundValidPair = true;
                        daughter1 = dau1;
                        daughter2 = dau2;

                        std::cout << "Event " << i << ": Found valid daughter pair - ";
                        std::cout << "Daughter1 (PDG:" << (true1 ? true1->PDG : -1) << ", TrueID:" << (true1 ? true1->ID : -1) << ") ";
                        std::cout << "Daughter2 (PDG:" << (true2 ? true2->PDG : -1) << ", TrueID:" << (true2 ? true2->ID : -1) << ") ";
                        std::cout << "Dist from beam end: " << dist1 << "cm, " << dist2 << "cm ";
                        std::cout << "Dist between daughters: " << distBetweenDaughters << "cm" << std::endl;
                        break;
                    }
                }
            }
            if (foundValidPair) break;
        }

        if (!foundValidPair) {
            std::cout << "DEBUG: No valid daughter pair found for this event" << std::endl;
        }

        // Only display events with valid daughter pair
        if (foundValidPair) {
            selectedEventCount++;
            std::cout << "=== Entry " << i << " (Selected Event #" << selectedEventCount << ") ===" << std::endl;
            AnaEventInfoB* eventInfo = spill->EventInfo;

            // We already have the reconstructed particles (daughter1 and daughter2)
            std::vector<AnaParticleB*> recoDaughter1;
            std::vector<AnaParticleB*> recoDaughter2;

            // Add the selected daughters to the vectors
            recoDaughter1.push_back(daughter1);
            recoDaughter2.push_back(daughter2);

            // Create new graphs for this event only
            std::map<int, TGraph*> eventParticleGraphs;
            std::set<AnaParticleB*> selectedParticles; // Track which particles are selected for highlighting

            // Add the selected daughters to selected particles for highlighting
            for (UInt_t idx = 0; idx < recoDaughter1.size(); idx++) {
                selectedParticles.insert(recoDaughter1[idx]);
            }
            for (UInt_t idx = 0; idx < recoDaughter2.size(); idx++) {
                selectedParticles.insert(recoDaughter2[idx]);
            }

            // Process reconstructed particles from all bunches
            for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
                AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
                if (!bunch) continue;

                for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                    AnaParticleB* particle = bunch->Particles[j];
                    if (!particle) continue;

                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                    if (!particlePD) continue;

                    // Get associated true particle
                    AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                    if (!truePart) continue;

                    int particleType = truePart->PDG;

                    // Create graph for this particle type if it doesn't exist
                    if (eventParticleGraphs.find(particleType) == eventParticleGraphs.end()) {
                        eventParticleGraphs[particleType] = new TGraph();
                        eventParticleGraphs[particleType]->SetMarkerStyle(20);
                        eventParticleGraphs[particleType]->SetMarkerSize(0.5);
                        if (particleColors.find(particleType) != particleColors.end()) {
                            eventParticleGraphs[particleType]->SetMarkerColor(particleColors[particleType]);
                        } else {
                            eventParticleGraphs[particleType]->SetMarkerColor(kGray);
                        }
                    }

                    // Add hits from all planes
                    for (int plane = 0; plane < 3; plane++) {
                        for (UInt_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];
                            if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                                int nPoints = eventParticleGraphs[particleType]->GetN();
                                eventParticleGraphs[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Y());
                            }
                        }
                    }
                }
            }

            // Create canvas title with daughter information
            AnaTrueParticleB* trueDaughter1 = daughter1->GetTrueParticle();
            AnaTrueParticleB* trueDaughter2 = daughter2->GetTrueParticle();
            std::string canvasTitle = "Event " + std::to_string(eventInfo->Event) + " (Selected Event #" + std::to_string(selectedEventCount) + ")";
            canvasTitle += " - Daughters: " + std::to_string(recoDaughter1.size()) + " (PDG:" + std::to_string(trueDaughter1 ? trueDaughter1->PDG : -1) + "), " + std::to_string(recoDaughter2.size()) + " (PDG:" + std::to_string(trueDaughter2 ? trueDaughter2->PDG : -1) + ")";

            TCanvas *c1 = new TCanvas("c1", canvasTitle.c_str(), 1400, 700);
            c1->Divide(2, 1); // Two subplots side by side

            // Create graphs for both projections
            std::map<int, TGraph*> eventParticleGraphsXY;
            std::map<int, TGraph*> eventParticleGraphsXZ;

            // Special graphs for highlighted daughters
            TGraph* highlightedDaughter1XY = nullptr;
            TGraph* highlightedDaughter2XY = nullptr;
            TGraph* highlightedDaughter1XZ = nullptr;
            TGraph* highlightedDaughter2XZ = nullptr;

            // Copy graphs for XY projection
            for (auto& pair : eventParticleGraphs) {
                eventParticleGraphsXY[pair.first] = new TGraph();
                eventParticleGraphsXY[pair.first]->SetMarkerStyle(20);
                eventParticleGraphsXY[pair.first]->SetMarkerSize(0.5);
                eventParticleGraphsXY[pair.first]->SetMarkerColor(pair.second->GetMarkerColor());
            }

            // Copy graphs for XZ projection
            for (auto& pair : eventParticleGraphs) {
                eventParticleGraphsXZ[pair.first] = new TGraph();
                eventParticleGraphsXZ[pair.first]->SetMarkerStyle(20);
                eventParticleGraphsXZ[pair.first]->SetMarkerSize(0.5);
                eventParticleGraphsXZ[pair.first]->SetMarkerColor(pair.second->GetMarkerColor());
            }

            // Create highlighted graphs for selected daughters
            if (recoDaughter1.size() > 0) {
                highlightedDaughter1XY = new TGraph();
                highlightedDaughter1XY->SetMarkerStyle(24); // Open circle
                highlightedDaughter1XY->SetMarkerSize(1.2);
                highlightedDaughter1XY->SetMarkerColor(kRed);
                highlightedDaughter1XY->SetLineColor(kRed);
                highlightedDaughter1XY->SetLineWidth(2);

                highlightedDaughter1XZ = new TGraph();
                highlightedDaughter1XZ->SetMarkerStyle(24);
                highlightedDaughter1XZ->SetMarkerSize(1.2);
                highlightedDaughter1XZ->SetMarkerColor(kRed);
                highlightedDaughter1XZ->SetLineColor(kRed);
                highlightedDaughter1XZ->SetLineWidth(2);
            }

            if (recoDaughter2.size() > 0) {
                highlightedDaughter2XY = new TGraph();
                highlightedDaughter2XY->SetMarkerStyle(25); // Open square
                highlightedDaughter2XY->SetMarkerSize(1.2);
                highlightedDaughter2XY->SetMarkerColor(kGreen);
                highlightedDaughter2XY->SetLineColor(kGreen);
                highlightedDaughter2XY->SetLineWidth(2);

                highlightedDaughter2XZ = new TGraph();
                highlightedDaughter2XZ->SetMarkerStyle(25);
                highlightedDaughter2XZ->SetMarkerSize(1.2);
                highlightedDaughter2XZ->SetMarkerColor(kGreen);
                highlightedDaughter2XZ->SetLineColor(kGreen);
                highlightedDaughter2XZ->SetLineWidth(2);
            }

            // Fill XY and XZ graphs with hit data
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

                    int particleType = truePart->PDG;

                    // Add hits to XY projection
                    for (int plane = 0; plane < 3; plane++) {
                        for (UInt_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];
                            if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                                int nPoints = eventParticleGraphsXY[particleType]->GetN();
                                eventParticleGraphsXY[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Y());

                                // Add to highlighted graphs if this is a selected daughter
                                if (selectedParticles.find(particle) != selectedParticles.end() && particle == daughter1 && highlightedDaughter1XY) {
                                    int nHighlighted = highlightedDaughter1XY->GetN();
                                    highlightedDaughter1XY->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Y());
                                } else if (selectedParticles.find(particle) != selectedParticles.end() && particle == daughter2 && highlightedDaughter2XY) {
                                    int nHighlighted = highlightedDaughter2XY->GetN();
                                    highlightedDaughter2XY->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Y());
                                }
                            }
                        }
                    }

                    // Add hits to XZ projection
                    for (int plane = 0; plane < 3; plane++) {
                        for (UInt_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];
                            if (hit.Position.X() > -900 && hit.Position.Z() > -900) {
                                int nPoints = eventParticleGraphsXZ[particleType]->GetN();
                                eventParticleGraphsXZ[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Z());

                                // Add to highlighted graphs if this is a selected daughter
                                if (selectedParticles.find(particle) != selectedParticles.end() && particle == daughter1 && highlightedDaughter1XZ) {
                                    int nHighlighted = highlightedDaughter1XZ->GetN();
                                    highlightedDaughter1XZ->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Z());
                                } else if (selectedParticles.find(particle) != selectedParticles.end() && particle == daughter2 && highlightedDaughter2XZ) {
                                    int nHighlighted = highlightedDaughter2XZ->GetN();
                                    highlightedDaughter2XZ->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Z());
                                }
                            }
                        }
                    }
                }
            }

            // Draw XY projection (left subplot)
            c1->cd(1);
            std::string histTitleXY = "Event " + std::to_string(eventInfo->Event) + " - XY Projection; X [cm]; Y [cm]";
            TH2F *dummyXY = new TH2F("dummyXY", histTitleXY.c_str(), 1000, -360, 360, 1000, 0, 700);
            dummyXY->SetStats(0); // Disable statistics box
            dummyXY->Draw();

            // Draw all particle graphs in XY projection
            for (auto& pair : eventParticleGraphsXY) {
                if (pair.second->GetN() > 0) {
                    pair.second->Draw("P SAME");
                }
            }

            // Draw highlighted daughters in XY projection
            if (highlightedDaughter1XY && highlightedDaughter1XY->GetN() > 0) {
                highlightedDaughter1XY->Draw("P SAME");
            }
            if (highlightedDaughter2XY && highlightedDaughter2XY->GetN() > 0) {
                highlightedDaughter2XY->Draw("P SAME");
            }

            // Draw beam particle trajectory in XY projection
            TLine* beamTrajectoryXY = new TLine(beamTrue->Position[0], beamTrue->Position[1],
                                               beamTrue->PositionEnd[0], beamTrue->PositionEnd[1]);
            beamTrajectoryXY->SetLineColor(kBlack);
            beamTrajectoryXY->SetLineWidth(3);
            beamTrajectoryXY->SetLineStyle(1);
            beamTrajectoryXY->Draw();

            // Draw XZ projection (right subplot)
            c1->cd(2);
            std::string histTitleXZ = "Event " + std::to_string(eventInfo->Event) + " - XZ Projection; X [cm]; Z [cm]";
            TH2F *dummyXZ = new TH2F("dummyXZ", histTitleXZ.c_str(), 1000, -360, 360, 1000, 0, 700);
            dummyXZ->SetStats(0); // Disable statistics box
            dummyXZ->Draw();

            // Draw all particle graphs in XZ projection
            for (auto& pair : eventParticleGraphsXZ) {
                if (pair.second->GetN() > 0) {
                    pair.second->Draw("P SAME");
                }
            }

            // Draw highlighted daughters in XZ projection
            if (highlightedDaughter1XZ && highlightedDaughter1XZ->GetN() > 0) {
                highlightedDaughter1XZ->Draw("P SAME");
            }
            if (highlightedDaughter2XZ && highlightedDaughter2XZ->GetN() > 0) {
                highlightedDaughter2XZ->Draw("P SAME");
            }

            // Draw beam particle trajectory in XZ projection
            TLine* beamTrajectoryXZ = new TLine(beamTrue->Position[0], beamTrue->Position[2],
                                               beamTrue->PositionEnd[0], beamTrue->PositionEnd[2]);
            beamTrajectoryXZ->SetLineColor(kBlack);
            beamTrajectoryXZ->SetLineWidth(3);
            beamTrajectoryXZ->SetLineStyle(1);
            beamTrajectoryXZ->Draw();

            // Add legend to the right subplot
            c1->cd(2);
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetTextSize(0.02);
            legend->SetBorderSize(1);
            legend->SetFillColor(kWhite);

            // Add entries for each particle type found
            for (auto& pair : eventParticleGraphsXY) {
                if (pair.second->GetN() > 0) {
                    std::string particleName = "PDG " + std::to_string(pair.first);
                    if (particleNames.find(pair.first) != particleNames.end()) {
                        particleName = particleNames[pair.first];
                    }
                    legend->AddEntry(pair.second, particleName.c_str(), "P");
                }
            }

            // Add highlighted daughters to legend
            if (highlightedDaughter1XY && highlightedDaughter1XY->GetN() > 0) {
                std::string daughter1Name = "Daughter1 (PDG:" + std::to_string(trueDaughter1 ? trueDaughter1->PDG : -1) + ")";
                legend->AddEntry(highlightedDaughter1XY, daughter1Name.c_str(), "P");
            }
            if (highlightedDaughter2XY && highlightedDaughter2XY->GetN() > 0) {
                std::string daughter2Name = "Daughter2 (PDG:" + std::to_string(trueDaughter2 ? trueDaughter2->PDG : -1) + ")";
                legend->AddEntry(highlightedDaughter2XY, daughter2Name.c_str(), "P");
            }

            // Add beam particle trajectory to legend
            TLine* beamLegendLine = new TLine(0, 0, 1, 1);
            beamLegendLine->SetLineColor(kBlack);
            beamLegendLine->SetLineWidth(3);
            beamLegendLine->SetLineStyle(1);
            legend->AddEntry(beamLegendLine, "Beam Particle", "l");

            legend->Draw();

            c1->Update();
            std::cout << "Canvas is ready for interaction. Press any key in the canvas to continue..." << std::endl;
            c1->WaitPrimitive("k");

            // Clean up graphs
            for (auto& pair : eventParticleGraphs) {
                delete pair.second;
            }
            for (auto& pair : eventParticleGraphsXY) {
                delete pair.second;
            }
            for (auto& pair : eventParticleGraphsXZ) {
                delete pair.second;
            }

            // Clean up highlighted graphs
            if (highlightedDaughter1XY) delete highlightedDaughter1XY;
            if (highlightedDaughter2XY) delete highlightedDaughter2XY;
            if (highlightedDaughter1XZ) delete highlightedDaughter1XZ;
            if (highlightedDaughter2XZ) delete highlightedDaughter2XZ;

            // Clean up beam trajectory lines
            // Note: TLine objects drawn on canvas are automatically cleaned up by ROOT
            // when the canvas is destroyed, so we don't need to manually delete them

        }

        std::cout << "DEBUG: Finished processing event " << i << std::endl;

        // Progress output every 1000 events (but we're only doing 3 events for debugging)
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << " (Selected events found: " << selectedEventCount << ")" << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total selected events found and displayed: " << selectedEventCount << std::endl;

    f->Close();
}
