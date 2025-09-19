// File: ReadSpill.C
// Modified to plot reconstructed particle hits (AnaParticlePD) instead of true particle trajectories
void ReadSpill(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
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
    int k0EventCount = 0;

    // Histogram for reconstructed invariant masses
    TH1F* hRecoInvariantMass = new TH1F("hRecoInvariantMass", "K^{0} Reconstructed Invariant Mass;Mass [MeV/c^{2}];Events", 50, 400, 600);
    hRecoInvariantMass->SetLineColor(kBlue);
    hRecoInvariantMass->SetLineWidth(2);
    hRecoInvariantMass->SetFillColor(kBlue);
    hRecoInvariantMass->SetFillStyle(3004);

    for (Long64_t i = 0; i < std::min<Long64_t>(1000, nEntries); i++) {
        MiniTree->GetEntry(i);

        // Get beam particle for highlighting
        AnaParticleB* beamParticle = nullptr;
        if (spill->Beam) {
            AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
            std::cout << "Beam object found, BeamParticle pointer: " << (beam ? beam->BeamParticle : nullptr) << std::endl;
            if (beam && beam->BeamParticle) {
                beamParticle = beam->BeamParticle;
                std::cout << "Found beam particle with UniqueID: " << beamParticle->UniqueID << std::endl;
            } else {
                std::cout << "No beam particle found in this event" << std::endl;
            }
        } else {
            std::cout << "No beam object found in this event" << std::endl;
        }

        UInt_t nTrueParticles = spill->TrueParticles.size();

        // Check if there's a K0 in this event
        bool hasK0 = false;
        for (UInt_t j = 0; j < nTrueParticles; j++) {
            AnaTrueParticleB* part = spill->TrueParticles[j];
            if (!part) continue;
            if (part->PDG == 310) {
                hasK0 = true;
                break;
            }
        }

        // If event has K0, plot this event individually
        if (hasK0) {
            k0EventCount++;
            std::cout << "=== Entry " << i << " (K0 Event #" << k0EventCount << ") ===" << std::endl;
            AnaEventInfoB* eventInfo = spill->EventInfo;
            std::cout << "Run: " << eventInfo->Run
                    << "  SubRun: " << eventInfo->SubRun
                    << "  Event: " << eventInfo->Event
                    << "  IsMC: " << eventInfo->IsMC
                    << "  IsSand: " << eventInfo->IsSand
                    << std::endl;
            std::cout << "Number of true particles: " << nTrueParticles << std::endl;

            // Create new graphs for this event only
            std::map<int, TGraph*> eventParticleGraphs;
            int particleCount = 0;
            int piPlusCount = 0;
            int piMinusCount = 0;

            // Process reconstructed particles from all bunches
            std::cout << "Processing reconstructed particles..." << std::endl;

            for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
                AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
                if (!bunch) continue;

                std::cout << "Bunch " << bunchIdx << " has " << bunch->Particles.size() << " reconstructed particles" << std::endl;

                for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                    AnaParticleB* particle = bunch->Particles[j];
                    if (!particle) continue;

                    // Try to cast to AnaParticlePD to access ProtoDUNE-specific data
                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                    if (!particlePD) {
                        std::cout << "   Particle " << j << " is not AnaParticlePD, skipping..." << std::endl;
                        continue;
                    }

                    // Print particle information
                    std::cout << "   Reconstructed Particle " << j
                                << " UniqueID=" << particlePD->UniqueID
                                << " NHits=" << particlePD->NHits
                                << " Length=" << particlePD->Length
                                << " Type=" << particlePD->Type
                                << " TrackID=" << particlePD->TrackID
                                << " ShowerID=" << particlePD->ShowerID;

                    // Get associated true particle if available
                    AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                    if (truePart) {
                        std::cout << " TruePDG=" << truePart->PDG;
                }
                std::cout << std::endl;

                    // Process hits for this particle
                    int totalHits = 0;
                    for (int plane = 0; plane < 3; plane++) {
                        totalHits += particlePD->Hits[plane].size();
                    }

                    if (totalHits == 0) {
                        std::cout << "     No hits found for this particle, skipping..." << std::endl;
                        continue;
                    }

                    // Create a unique identifier for this particle type
                    int particleType = 0; // Default for unknown
                    if (truePart) {
                        particleType = truePart->PDG;
                    } else {
                        // Use a combination of type and IDs for unidentified particles
                        particleType = 1000 + particlePD->Type + particlePD->TrackID;
                    }

                    // Include neutrons and protons now
                    // if (truePart && (truePart->PDG == 2112 || truePart->PDG == 2212)) {
                    //     std::cout << "     [SKIPPED - neutron/proton]" << std::endl;
                    //     continue;
                    // }

                    // Create or get TGraph for this particle type in this event
                    if (eventParticleGraphs.find(particleType) == eventParticleGraphs.end()) {
                        eventParticleGraphs[particleType] = new TGraph();
                        eventParticleGraphs[particleType]->SetMarkerStyle(20);
                        eventParticleGraphs[particleType]->SetMarkerSize(0.5);
                        eventParticleGraphs[particleType]->SetLineWidth(1);

                        // Set color based on true particle type or default
                        int particleColor = kGray; // Default color
                        if (truePart && particleColors.find(truePart->PDG) != particleColors.end()) {
                            particleColor = particleColors[truePart->PDG];
                        } else {
                            // Default colors for unidentified particles
                            int colors[] = {kGray, kGray+1, kGray+2, kGray+3};
                            particleColor = colors[particleType % 4];
                        }

                        // Check if this is the beam particle and highlight it
                        if (beamParticle && particlePD->UniqueID == beamParticle->UniqueID) {
                            particleColor = kOrange; // Orange color for beam particle
                            std::cout << "Highlighting beam particle (UniqueID=" << particlePD->UniqueID << ", PDG=" << (truePart ? truePart->PDG : -1) << ")" << std::endl;
                        }

                        eventParticleGraphs[particleType]->SetLineColor(particleColor);
                        eventParticleGraphs[particleType]->SetMarkerColor(particleColor);
                    }

                    // Add hits from all planes to the graph
                    int pointsAdded = 0;
                    for (int plane = 0; plane < 3; plane++) {
                        for (size_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];

                            // Use the 3D position coordinates from the hit
                            float x = hit.Position.X();    // X coordinate in cm
                            float y = hit.Position.Y();    // Y coordinate in cm

                            // Add point if we have valid coordinates
                            if (x != -999 && y != -999) {
                                int nPoints = eventParticleGraphs[particleType]->GetN();
                                eventParticleGraphs[particleType]->SetPoint(nPoints, x, y);
                                pointsAdded++;
                            }
                        }
                    }

                    // If no hits were added, try using the particle's position information as fallback
                    if (pointsAdded == 0) {
                        std::cout << "     No valid hit data found, using particle position as fallback" << std::endl;

                        // Use the particle's start and end positions
                        float startX = particlePD->PositionStart[0];
                        float startY = particlePD->PositionStart[1];
                        float endX = particlePD->PositionEnd[0];
                        float endY = particlePD->PositionEnd[1];

                        if (startX > -900 && startY > -900 && endX > -900 && endY > -900) {
                            int nPoints = eventParticleGraphs[particleType]->GetN();
                            eventParticleGraphs[particleType]->SetPoint(nPoints, startX, startY);
                            eventParticleGraphs[particleType]->SetPoint(nPoints + 1, endX, endY);
                            pointsAdded = 2;
                            std::cout << "     Added particle trajectory: (" << startX << "," << startY << ") to (" << endX << "," << endY << ")" << std::endl;
                        }
                    }

                    std::cout << "     Actually added " << pointsAdded << " points to graph" << std::endl;

                    // Count particles for verification
                    if (truePart) {
                        if (truePart->PDG == 211) piPlusCount++;
                        if (truePart->PDG == -211) piMinusCount++;
                    }
                    particleCount++;

                    std::cout << "     Added " << totalHits << " hits from " << 3 << " planes" << std::endl;

                    // Print some hit details for debugging
                    if (totalHits > 0) {
                        std::cout << "     Hit details (first few hits):" << std::endl;
                        int hitCount = 0;
                        for (int plane = 0; plane < 3 && hitCount < 5; plane++) {
                            for (size_t h = 0; h < particlePD->Hits[plane].size() && hitCount < 5; h++) {
                                AnaHitPD hit = particlePD->Hits[plane][h];
                                std::cout << "       Plane " << plane
                                          << " Position (" << hit.Position.X() << ", " << hit.Position.Y() << ", " << hit.Position.Z() << ")"
                                          << " Wire " << hit.WireID.Wire
                                          << " PeakTime " << hit.PeakTime
                                          << " PeakAmp " << hit.PeakAmplitude
                                          << " Integral " << hit.Integral << std::endl;
                                hitCount++;
                            }
                        }
                    }
                }
            }

        // Count K0 daughters and their PDG codes
        std::map<int, int> k0DaughterCounts;
        int totalK0Daughters = 0;

        // Count K0 parent daughters and their PDG codes
        std::map<int, int> k0ParentDaughterCounts;
        int totalK0ParentDaughters = 0;
        int k0ParentPDG = 0;

        // Count number of K0s in this event
        int k0Count = 0;

        // Calculate K0 invariant mass
        double k0InvariantMass = -1.0; // Default value if no K0 found

        // Calculate K0 invariant mass from reconstructed daughters
        double k0RecoInvariantMass = -1.0; // Default value if no K0 daughters found

            for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

            // Count K0s in this event and calculate invariant mass
            if (truePart->PDG == 310) {
                k0Count++;

                // Calculate invariant mass from momentum and direction
                if (truePart->Momentum > 0 &&
                    truePart->Direction[0] != 0 && truePart->Direction[1] != 0 && truePart->Direction[2] != 0) {

                    // Calculate momentum components
                    double px = truePart->Momentum * truePart->Direction[0];
                    double py = truePart->Momentum * truePart->Direction[1];
                    double pz = truePart->Momentum * truePart->Direction[2];

                    // Calculate momentum magnitude
                    double p = sqrt(px*px + py*py + pz*pz);

                    // For true particles, we can calculate the invariant mass directly
                    // using the 4-momentum: m² = E² - p²
                    // Since this is a true particle, we know its mass should be the K0 mass
                    double k0Mass = 497.611; // MeV/c² (expected K0 mass)
                    double energy = sqrt(p*p + k0Mass*k0Mass);

                    // Calculate invariant mass: m² = E² - p²
                    k0InvariantMass = sqrt(energy*energy - p*p);

                    std::cout << "K0 momentum: " << truePart->Momentum
                              << " MeV/c, direction: (" << truePart->Direction[0]
                              << ", " << truePart->Direction[1]
                              << ", " << truePart->Direction[2] << ")" << std::endl;
                    std::cout << "K0 calculated momentum magnitude: " << p << " MeV/c" << std::endl;
                    std::cout << "K0 energy: " << energy << " MeV" << std::endl;
                    std::cout << "K0 invariant mass: " << k0InvariantMass << " MeV/c² (expected: " << k0Mass << " MeV/c²)" << std::endl;
                }
            }

                // Check if this particle is a daughter of a K0
                if (truePart->ParentPDG == 310) {
                    k0DaughterCounts[truePart->PDG]++;
                    totalK0Daughters++;
                }

                // Check if this particle is a daughter of a K0 parent
                if (truePart->PDG == 310 && truePart->ParentPDG != 0) {
                    k0ParentPDG = truePart->ParentPDG;
                }
            }

            // Count daughters of the K0 parent
            if (k0ParentPDG != 0) {
                for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
                    AnaTrueParticleB* truePart = spill->TrueParticles[j];
                    if (!truePart) continue;

                    if (truePart->ParentPDG == k0ParentPDG) {
                        k0ParentDaughterCounts[truePart->PDG]++;
                        totalK0ParentDaughters++;
                    }
                }
            }

        // Create canvas title with K0 count, invariant mass, and daughter information
        std::string canvasTitle = "Event " + std::to_string(eventInfo->Event) + " (K0 Event #" + std::to_string(k0EventCount) + ", " + std::to_string(k0Count) + " K^{0}s";

        // Add invariant mass if calculated
        if (k0InvariantMass > 0) {
            // Format the invariant mass to 2 decimal places
            char massStr[50];
            snprintf(massStr, sizeof(massStr), "%.2f", k0InvariantMass);
            canvasTitle += ", m_{K^{0}}^{true} = " + std::string(massStr) + " MeV/c^{2}";
        }

        // Add reconstructed invariant mass if calculated
        if (k0RecoInvariantMass > 0) {
            // Format the reconstructed invariant mass to 2 decimal places
            char recoMassStr[50];
            snprintf(recoMassStr, sizeof(recoMassStr), "%.2f", k0RecoInvariantMass);
            canvasTitle += ", m_{K^{0}}^{reco} = " + std::string(recoMassStr) + " MeV/c^{2}";
        }
        canvasTitle += ")";

        // Add K0 daughters information
        if (totalK0Daughters > 0) {
            canvasTitle += " - K0 Daughters: " + std::to_string(totalK0Daughters) + " (";
            bool first = true;
            for (auto& pair : k0DaughterCounts) {
                if (!first) canvasTitle += ", ";
                if (particleNames.find(pair.first) != particleNames.end()) {
                    canvasTitle += particleNames[pair.first] + ":" + std::to_string(pair.second);
                } else {
                    canvasTitle += "PDG" + std::to_string(pair.first) + ":" + std::to_string(pair.second);
                }
                first = false;
            }
            canvasTitle += ")";
        }

            TCanvas *c1 = new TCanvas("c1", canvasTitle.c_str(), 1400, 700);
            c1->Divide(2, 1); // Two subplots side by side

            // Create graphs for both projections
            std::map<int, TGraph*> eventParticleGraphsXY;
            std::map<int, TGraph*> eventParticleGraphsXZ;

            // Copy the existing graphs for both projections
            for (auto& pair : eventParticleGraphs) {
                int pdg = pair.first;
                TGraph* originalGraph = pair.second;

                // Create XY projection graph
                eventParticleGraphsXY[pdg] = new TGraph();
                eventParticleGraphsXY[pdg]->SetMarkerStyle(20);
                eventParticleGraphsXY[pdg]->SetMarkerSize(0.5);
                eventParticleGraphsXY[pdg]->SetLineWidth(1);
                eventParticleGraphsXY[pdg]->SetLineColor(originalGraph->GetLineColor());
                eventParticleGraphsXY[pdg]->SetMarkerColor(originalGraph->GetMarkerColor());

                // Create XZ projection graph
                eventParticleGraphsXZ[pdg] = new TGraph();
                eventParticleGraphsXZ[pdg]->SetMarkerStyle(20);
                eventParticleGraphsXZ[pdg]->SetMarkerSize(0.5);
                eventParticleGraphsXZ[pdg]->SetLineWidth(1);
                eventParticleGraphsXZ[pdg]->SetLineColor(originalGraph->GetLineColor());
                eventParticleGraphsXZ[pdg]->SetMarkerColor(originalGraph->GetMarkerColor());

                // Copy points for both projections
                for (int i = 0; i < originalGraph->GetN(); i++) {
                    double x, y;
                    originalGraph->GetPoint(i, x, y);

                    // For XY projection, we need to get the Z coordinate from the hits
                    // We'll need to find the corresponding hit to get Z
                    // For now, let's use a simplified approach
                    eventParticleGraphsXY[pdg]->SetPoint(i, x, y);
                    eventParticleGraphsXZ[pdg]->SetPoint(i, x, y); // This will be updated below
                }
            }

            // Now we need to recreate the points with proper coordinates for both projections
            // Clear the XZ graphs and rebuild them
            for (auto& pair : eventParticleGraphsXZ) {
                pair.second->Set(0); // Clear all points
            }

            // Rebuild both projections with correct coordinates
            for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
                AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
                if (!bunch) continue;

                for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
                    AnaParticleB* particle = bunch->Particles[j];
                    if (!particle) continue;

                    AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
                    if (!particlePD) continue;

                    // Get associated true particle if available
                    AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
                    int particleType = 0;
                    if (truePart) {
                        particleType = truePart->PDG;
                    } else {
                        particleType = 1000 + particlePD->Type + particlePD->TrackID;
                    }

                    // Skip if this particle type should be skipped
                    if (truePart && (truePart->PDG == 2112 || truePart->PDG == 2212)) {
                        continue;
                    }

                    // Add hits to both projections
                    for (int plane = 0; plane < 3; plane++) {
                        for (size_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];

                            float x = hit.Position.X();
                            float y = hit.Position.Y();
                            float z = hit.Position.Z();

                            // Add point if we have valid coordinates
                            if (x != -999 && y != -999 && z != -999) {
                                // XY projection (X vs Y)
                                int nPointsXY = eventParticleGraphsXY[particleType]->GetN();
                                eventParticleGraphsXY[particleType]->SetPoint(nPointsXY, x, y);

                                // XZ projection (X vs Z)
                                int nPointsXZ = eventParticleGraphsXZ[particleType]->GetN();
                                eventParticleGraphsXZ[particleType]->SetPoint(nPointsXZ, x, z);
                            }
                        }
                    }
                }
            }

        // Add K0 trajectories and parent-K0 connection lines from true particles
        std::vector<TLine*> k0LinesXY;
        std::vector<TLine*> k0LinesXZ;
        std::vector<TLine*> parentK0LinesXY;
        std::vector<TLine*> parentK0LinesXZ;
        std::vector<TMarker*> k0StartMarkersXY;
        std::vector<TMarker*> k0StartMarkersXZ;
        std::vector<TMarker*> k0EndMarkersXY;
        std::vector<TMarker*> k0EndMarkersXZ;
        std::set<Int_t> k0ParentIDs; // Track unique parent IDs
        std::map<Int_t, std::vector<float>> k0EndPositions; // Map parent ID to K0 end positions

    // Add true particles as markers
    std::vector<TMarker*> trueParticleMarkersXY;
    std::vector<TMarker*> trueParticleMarkersXZ;

    for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
        AnaTrueParticleB* truePart = spill->TrueParticles[j];
        if (!truePart) continue;

        if (truePart->PDG == 310) { // K0 particle
            float startX = truePart->Position[0];
            float startY = truePart->Position[1];
            float startZ = truePart->Position[2];
            float endX = truePart->PositionEnd[0];
            float endY = truePart->PositionEnd[1];
            float endZ = truePart->PositionEnd[2];

                // Add K0 trajectory line if both start and end coordinates are valid
                if (startX > -900 && startY > -900 && startZ > -900 &&
                    endX > -900 && endY > -900 && endZ > -900) {

                    // Create magenta lines for K0 trajectories
                    TLine* lineXY = new TLine(startX, startY, endX, endY);
                    lineXY->SetLineColor(kMagenta);
                    lineXY->SetLineWidth(4);
                    k0LinesXY.push_back(lineXY);

                    TLine* lineXZ = new TLine(startX, startZ, endX, endZ);
                    lineXZ->SetLineColor(kMagenta);
                    lineXZ->SetLineWidth(4);
                    k0LinesXZ.push_back(lineXZ);

                    // Add large, distinctive markers for K0 start (production) and end (decay) points
                    // Start point (production) - Green square
                    TMarker* startMarkerXY = new TMarker(startX, startY, 22); // Square marker
                    startMarkerXY->SetMarkerColor(kGreen);
                    startMarkerXY->SetMarkerSize(2.0);
                    k0StartMarkersXY.push_back(startMarkerXY);

                    TMarker* startMarkerXZ = new TMarker(startX, startZ, 22); // Square marker
                    startMarkerXZ->SetMarkerColor(kGreen);
                    startMarkerXZ->SetMarkerSize(2.0);
                    k0StartMarkersXZ.push_back(startMarkerXZ);

                    // End point (decay) - Red circle
                    TMarker* endMarkerXY = new TMarker(endX, endY, 20); // Circle marker
                    endMarkerXY->SetMarkerColor(kRed);
                    endMarkerXY->SetMarkerSize(2.0);
                    k0EndMarkersXY.push_back(endMarkerXY);

                    TMarker* endMarkerXZ = new TMarker(endX, endZ, 20); // Circle marker
                    endMarkerXZ->SetMarkerColor(kRed);
                    endMarkerXZ->SetMarkerSize(2.0);
                    k0EndMarkersXZ.push_back(endMarkerXZ);

                    std::cout << "Added K0 trajectory: (" << startX << ", " << startY << ", " << startZ << ") to ("
                              << endX << ", " << endY << ", " << endZ << ")" << std::endl;
                    std::cout << "K0 Parent ID: " << truePart->ParentID << ", Parent PDG: " << truePart->ParentPDG << std::endl;

                    // Store parent ID and K0 end position for connection lines
                    if (truePart->ParentID > 0) {
                        k0ParentIDs.insert(truePart->ParentID);
                        // Store K0 end position for parent connection
                        k0EndPositions[truePart->ParentID] = {endX, endY, endZ};
                    }
                }
        }
    }

    // Find and create parent-K0 connection lines
    for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
        AnaTrueParticleB* truePart = spill->TrueParticles[j];
        if (!truePart) continue;

        // Check if this particle is a parent of a K0
        if (k0ParentIDs.find(truePart->ID) != k0ParentIDs.end()) {
            float parentX = truePart->PositionEnd[0];
            float parentY = truePart->PositionEnd[1];
            float parentZ = truePart->PositionEnd[2];

                // Get the K0 end position for this parent
                if (k0EndPositions.find(truePart->ID) != k0EndPositions.end()) {
                    std::vector<float> k0End = k0EndPositions[truePart->ID];
                    float k0EndX = k0End[0];
                    float k0EndY = k0End[1];
                    float k0EndZ = k0End[2];

                    // Create connection lines if coordinates are valid
                    if (parentX > -900 && parentY > -900 && parentZ > -900 &&
                        k0EndX > -900 && k0EndY > -900 && k0EndZ > -900) {

                        // Create blue lines connecting parent end to K0 end
                        TLine* lineXY = new TLine(parentX, parentY, k0EndX, k0EndY);
                        lineXY->SetLineColor(kBlue);
                        lineXY->SetLineWidth(2);
                        lineXY->SetLineStyle(2); // Dashed line
                        parentK0LinesXY.push_back(lineXY);

                        TLine* lineXZ = new TLine(parentX, parentZ, k0EndX, k0EndZ);
                        lineXZ->SetLineColor(kBlue);
                        lineXZ->SetLineWidth(2);
                        lineXZ->SetLineStyle(2); // Dashed line
                        parentK0LinesXZ.push_back(lineXZ);

                        std::cout << "Added parent-K0 connection: parent end (" << parentX << ", " << parentY << ", " << parentZ
                                 << ") to K0 end (" << k0EndX << ", " << k0EndY << ", " << k0EndZ << ")" << std::endl;
                    }
                }
        }
        }

        // Calculate K0 invariant mass from reconstructed daughters
        std::vector<AnaParticlePD*> k0Daughters;
        AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[0]);
        if (bunch) {
            for (UInt_t i = 0; i < bunch->Particles.size(); i++) {
                AnaParticleB* particle = bunch->Particles[i];
                if (!particle) continue;

            AnaParticlePD* particlePD = dynamic_cast<AnaParticlePD*>(particle);
            if (!particlePD) continue;

            // Get associated true particle
            AnaTrueParticleB* truePart = particlePD->GetTrueParticle();
            if (!truePart) continue;

            // Check if this is a π+ or π- that is a daughter of a K0
            if ((truePart->PDG == 211 || truePart->PDG == -211) && truePart->ParentPDG == 310) {
                // Add all K0 daughters, we'll handle momentum validation later
                k0Daughters.push_back(particlePD);
                std::cout << "Found K0 daughter: PDG=" << truePart->PDG
                          << ", UniqueID=" << particlePD->UniqueID
                          << ", Momentum=" << particlePD->Momentum << " MeV/c" << std::endl;
            }
        }
        }

        // Calculate K0 invariant mass from daughters
        if (k0Daughters.size() >= 2) {
            // If we have more than 2 daughters, select the first 2 with valid momentum
            if (k0Daughters.size() > 2) {
                std::cout << "Found " << k0Daughters.size() << " K0 daughters, selecting first 2 with valid momentum" << std::endl;
                k0Daughters.resize(2); // Keep only the first 2
            }
            AnaParticlePD* daughter1 = k0Daughters[0];
            AnaParticlePD* daughter2 = k0Daughters[1];

            // Get momentum and direction for both daughters
            double p1 = daughter1->Momentum;
            double p2 = daughter2->Momentum;

            // If momentum is invalid, use default values
            if (p1 <= 0 || p1 >= 10000) {
                std::cout << "Invalid momentum for daughter 1 (" << p1 << "), using default 200 MeV/c" << std::endl;
                p1 = 200.0; // Default value for pion momentum
            }

            if (p2 <= 0 || p2 >= 10000) {
                std::cout << "Invalid momentum for daughter 2 (" << p2 << "), using default 200 MeV/c" << std::endl;
                p2 = 200.0; // Default value for pion momentum
            }

            double dir1[3] = {daughter1->DirectionStart[0], daughter1->DirectionStart[1], daughter1->DirectionStart[2]};
            double dir2[3] = {daughter2->DirectionStart[0], daughter2->DirectionStart[1], daughter2->DirectionStart[2]};

            // Calculate 4-momentum components for each daughter
            double px1 = p1 * dir1[0];
            double py1 = p1 * dir1[1];
            double pz1 = p1 * dir1[2];
            double px2 = p2 * dir2[0];
            double py2 = p2 * dir2[1];
            double pz2 = p2 * dir2[2];

            // Calculate energy for each daughter (assuming pion mass = 139.57 MeV/c²)
            double pionMass = 139.57; // MeV/c²
            double E1 = sqrt(p1*p1 + pionMass*pionMass);
            double E2 = sqrt(p2*p2 + pionMass*pionMass);

            // Calculate total 4-momentum of K0
            double totalPx = px1 + px2;
            double totalPy = py1 + py2;
            double totalPz = pz1 + pz2;
            double totalE = E1 + E2;

            // Calculate invariant mass: m² = E² - p²
            double totalP2 = totalPx*totalPx + totalPy*totalPy + totalPz*totalPz;
            k0RecoInvariantMass = sqrt(totalE*totalE - totalP2);

            std::cout << "K0 reconstructed invariant mass calculation:" << std::endl;
            std::cout << "  Daughter 1: p=" << p1 << " MeV/c, E=" << E1 << " MeV" << std::endl;
            std::cout << "  Daughter 2: p=" << p2 << " MeV/c, E=" << E2 << " MeV" << std::endl;
            std::cout << "  Total 4-momentum: (" << totalPx << ", " << totalPy << ", " << totalPz << ", " << totalE << ")" << std::endl;
            std::cout << "  K0 reconstructed invariant mass: " << k0RecoInvariantMass << " MeV/c²" << std::endl;

            // Fill histogram with reconstructed invariant mass
            hRecoInvariantMass->Fill(k0RecoInvariantMass);
        } else if (k0Daughters.size() > 0) {
            std::cout << "Found " << k0Daughters.size() << " K0 daughters, expected exactly 2 for invariant mass calculation" << std::endl;
        }

        // Add true particles as markers (excluding K0 which is already shown as trajectory)
        for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
            AnaTrueParticleB* truePart = spill->TrueParticles[j];
            if (!truePart) continue;

            // Skip K0 particles as they are shown as trajectories
            if (truePart->PDG == 310) continue;

        // Get particle position (use start position)
        float x = truePart->Position[0];
        float y = truePart->Position[1];
        float z = truePart->Position[2];

        // Only add if coordinates are valid
        if (x > -900 && y > -900 && z > -900) {
            // Get particle color and name
            int color = kGray; // Default color
            if (particleColors.find(truePart->PDG) != particleColors.end()) {
                color = particleColors[truePart->PDG];
            }

            // Create markers for both projections
            TMarker* markerXY = new TMarker(x, y, 20); // Circle marker
            markerXY->SetMarkerColor(color);
            markerXY->SetMarkerSize(0.8);
            trueParticleMarkersXY.push_back(markerXY);

            TMarker* markerXZ = new TMarker(x, z, 20); // Circle marker
            markerXZ->SetMarkerColor(color);
            markerXZ->SetMarkerSize(0.8);
            trueParticleMarkersXZ.push_back(markerXZ);

            std::cout << "Added true particle (PDG=" << truePart->PDG << ") at: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }
    }

            // Draw XY projection (left subplot)
            c1->cd(1);
            std::string histTitleXY = "Event " + std::to_string(eventInfo->Event) + " - XY Projection; X [cm]; Y [cm]";
            TH2F *dummyXY = new TH2F("dummyXY", histTitleXY.c_str(), 1000, -360, 360, 1000, 0, 700);
            dummyXY->SetStats(0);
            dummyXY->Draw("AXIS");

            for (auto& pair : eventParticleGraphsXY) {
                int pdg = pair.first;
                TGraph* graph = pair.second;

                if (graph->GetN() > 0) {
                    graph->Draw("P SAME");
                    std::cout << "Drew XY graph for PDG " << pdg << " with " << graph->GetN() << " points" << std::endl;
                }
            }

        // Draw K0 trajectories as red lines in XY projection
        for (auto* line : k0LinesXY) {
            line->Draw("SAME");
        }

        // Draw K0 production and decay point markers in XY projection
        for (auto* marker : k0StartMarkersXY) {
            marker->Draw("SAME");
        }
        for (auto* marker : k0EndMarkersXY) {
            marker->Draw("SAME");
        }

        // Draw parent-K0 connection lines in XY projection
        for (auto* line : parentK0LinesXY) {
            line->Draw("SAME");
        }

        // Draw true particle markers in XY projection
        for (auto* marker : trueParticleMarkersXY) {
            marker->Draw("SAME");
        }

            // Draw XZ projection (right subplot)
            c1->cd(2);
            std::string histTitleXZ = "Event " + std::to_string(eventInfo->Event) + " - XZ Projection; X [cm]; Z [cm]";
            TH2F *dummyXZ = new TH2F("dummyXZ", histTitleXZ.c_str(), 1000, -360, 360, 1000, 0, 700);
            dummyXZ->SetStats(0);
            dummyXZ->Draw("AXIS");

            for (auto& pair : eventParticleGraphsXZ) {
                int pdg = pair.first;
                TGraph* graph = pair.second;

                if (graph->GetN() > 0) {
                    graph->Draw("P SAME");
                    std::cout << "Drew XZ graph for PDG " << pdg << " with " << graph->GetN() << " points" << std::endl;
                }
            }

        // Draw K0 trajectories as red lines in XZ projection
        for (auto* line : k0LinesXZ) {
            line->Draw("SAME");
        }

        // Draw K0 production and decay point markers in XZ projection
        for (auto* marker : k0StartMarkersXZ) {
            marker->Draw("SAME");
        }
        for (auto* marker : k0EndMarkersXZ) {
            marker->Draw("SAME");
        }

        // Draw parent-K0 connection lines in XZ projection
        for (auto* line : parentK0LinesXZ) {
            line->Draw("SAME");
        }

        // Draw true particle markers in XZ projection
        for (auto* marker : trueParticleMarkersXZ) {
            marker->Draw("SAME");
        }

            // Add legend to the right subplot
            c1->cd(2);
            TLegend *legend = new TLegend(0.70, 0.65, 0.95, 0.95);
            legend->SetBorderSize(1);
            legend->SetTextSize(0.030);
            legend->SetFillColor(kWhite);
            legend->SetFillStyle(1001);
            legend->SetTextFont(42);
            legend->SetNColumns(2);

            for (auto& pair : eventParticleGraphs) {
                int pdg = pair.first;
                TGraph* graph = pair.second;

                if (graph->GetN() > 0) {
                    std::string name;
                    if (particleNames.find(pdg) != particleNames.end()) {
                        name = particleNames[pdg];
                    } else if (pdg >= 1000) {
                        name = "Unidentified " + std::to_string(pdg - 1000);
                    } else {
                        name = "PDG " + std::to_string(pdg);
                    }
                legend->AddEntry(graph, name.c_str(), "LP");
                }
            }

        // Add K0 trajectory to legend if any K0 particles exist
        if (!k0LinesXY.empty()) {
            TLine* k0Line = new TLine(0, 0, 1, 1);
            k0Line->SetLineColor(kMagenta);
            k0Line->SetLineWidth(4);
            legend->AddEntry(k0Line, "K^{0}", "L");
        }

        // Add K0 production and decay point markers to legend
        if (!k0StartMarkersXY.empty()) {
            TMarker* startMarker = new TMarker(0, 0, 22);
            startMarker->SetMarkerColor(kGreen);
            startMarker->SetMarkerSize(2.0);
            legend->AddEntry(startMarker, "K^{0} Start", "P");
        }
        if (!k0EndMarkersXY.empty()) {
            TMarker* endMarker = new TMarker(0, 0, 20);
            endMarker->SetMarkerColor(kRed);
            endMarker->SetMarkerSize(2.0);
            legend->AddEntry(endMarker, "K^{0} End", "P");
        }

        // Add beam particle to legend if found
        if (beamParticle) {
            TMarker* beamMarker = new TMarker(0, 0, 20);
            beamMarker->SetMarkerColor(kOrange);
            beamMarker->SetMarkerSize(1.0);
            legend->AddEntry(beamMarker, "Beam", "P");
        }

        // Add parent-K0 connection to legend if any connection lines exist
        if (!parentK0LinesXY.empty()) {
            TLine* parentLine = new TLine(0, 0, 1, 1);
            parentLine->SetLineColor(kBlue);
            parentLine->SetLineWidth(2);
            parentLine->SetLineStyle(2);

            // Get the parent PDG from the first K0 found and use the particle name
            std::string parentLabel = "K^{0} parent connection";
            for (UInt_t j = 0; j < spill->TrueParticles.size(); j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;
                if (truePart->PDG == 310 && truePart->ParentPDG != 0) {
                    // Look up the particle name for the parent PDG
                    if (particleNames.find(truePart->ParentPDG) != particleNames.end()) {
                        parentLabel = particleNames[truePart->ParentPDG] + " #rightarrow K^{0}";
                    } else {
                        parentLabel = "Parent #rightarrow K^{0}";
                    }
                    break;
                }
            }
            legend->AddEntry(parentLine, parentLabel.c_str(), "L");
        }

            legend->Draw();

            // Add grid for better visualization
            c1->SetGrid();

            // Update canvas to show the plot
            c1->Update();

            // Force canvas to be displayed
            gSystem->ProcessEvents();

            std::cout << "Reconstructed particles plotted in this event: " << particleCount << std::endl;
            std::cout << "π⁺ reconstructed particles plotted: " << piPlusCount << std::endl;
            std::cout << "π⁻ reconstructed particles plotted: " << piMinusCount << std::endl;
            std::cout << "Canvas is ready for interaction. Press any key in the canvas to continue..." << std::endl;

            // Wait for any key press in the canvas (this is the proper ROOT way)
            c1->WaitPrimitive("k");

            // Clean up this event's graphs
            for (auto& pair : eventParticleGraphs) {
                delete pair.second;
            }
            eventParticleGraphs.clear();

            for (auto& pair : eventParticleGraphsXY) {
                delete pair.second;
            }
            eventParticleGraphsXY.clear();

            for (auto& pair : eventParticleGraphsXZ) {
                delete pair.second;
            }
            eventParticleGraphsXZ.clear();

        // Clean up K0 lines
        for (auto* line : k0LinesXY) {
            delete line;
        }
        k0LinesXY.clear();

        for (auto* line : k0LinesXZ) {
            delete line;
        }
        k0LinesXZ.clear();

        // Clean up parent-K0 connection lines
        for (auto* line : parentK0LinesXY) {
            delete line;
        }
        parentK0LinesXY.clear();

        for (auto* line : parentK0LinesXZ) {
            delete line;
        }
        parentK0LinesXZ.clear();

        // Clean up K0 production and decay markers
        for (auto* marker : k0StartMarkersXY) {
            delete marker;
        }
        k0StartMarkersXY.clear();

        for (auto* marker : k0StartMarkersXZ) {
            delete marker;
        }
        k0StartMarkersXZ.clear();

        for (auto* marker : k0EndMarkersXY) {
            delete marker;
        }
        k0EndMarkersXY.clear();

        for (auto* marker : k0EndMarkersXZ) {
            delete marker;
        }
        k0EndMarkersXZ.clear();

        // Clean up true particle markers
        for (auto* marker : trueParticleMarkersXY) {
            delete marker;
        }
        trueParticleMarkersXY.clear();

        for (auto* marker : trueParticleMarkersXZ) {
            delete marker;
        }
        trueParticleMarkersXZ.clear();

            std::cout << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total K0 events found and displayed: " << k0EventCount << std::endl;
    std::cout << "Total reconstructed invariant masses filled: " << hRecoInvariantMass->GetEntries() << std::endl;

    // Display the reconstructed invariant mass histogram
    TCanvas *c2 = new TCanvas("c2", "K^{0} Reconstructed Invariant Mass Distribution", 800, 600);
    c2->cd();
    hRecoInvariantMass->Draw("HIST");

    // Add statistics box
    gStyle->SetOptStat(1111);
    gPad->Update();

    // Add a line for the expected K0 mass
    TLine *expectedMass = new TLine(497.611, 0, 497.611, hRecoInvariantMass->GetMaximum() * 1.1);
    expectedMass->SetLineColor(kRed);
    expectedMass->SetLineWidth(2);
    expectedMass->SetLineStyle(2);
    expectedMass->Draw();

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hRecoInvariantMass, "Reconstructed Mass", "f");
    leg->AddEntry(expectedMass, "Expected K^{0} Mass", "l");
    leg->Draw();

    c2->Update();
    std::cout << "Reconstructed invariant mass histogram displayed. Press any key to continue..." << std::endl;
    c2->WaitPrimitive("k");

    f->Close();
}
