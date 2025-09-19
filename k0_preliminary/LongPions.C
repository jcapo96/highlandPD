// File: LongPions.C
// Modified to plot reconstructed particle hits (AnaParticlePD) instead of true particle trajectories
// Only displays events with K0→π+π- decay based on TRUTH information, but shows RECO information
void LongPions(const char* filename = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/6GeV_prod4a_01_minitree_2023-01-27.root") {
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
    int longPionEventCount = 0;

    // Histogram for reconstructed invariant masses
    TH1F* hRecoInvariantMass = new TH1F("hRecoInvariantMass", "K^{0} Reconstructed Invariant Mass;Mass [MeV/c^{2}];Events", 50, 400, 600);
    hRecoInvariantMass->SetLineColor(kBlue);
    hRecoInvariantMass->SetLineWidth(2);
    hRecoInvariantMass->SetFillColor(kBlue);
    hRecoInvariantMass->SetFillStyle(3004);

    for (Long64_t i = 0; i < nEntries; i++) {
        MiniTree->GetEntry(i);

        // Get beam particle for highlighting
        AnaParticleB* beamParticle = nullptr;
        if (spill->Beam) {
            AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
            if (beam && beam->BeamParticle) {
                beamParticle = beam->BeamParticle;
            }
        }

        UInt_t nTrueParticles = spill->TrueParticles.size();

        // Check if there's a K0 with exactly two daughters (one π+ and one π-) based on TRUTH information
        bool hasK0WithExactPiDaughters = false;
        std::vector<AnaTrueParticleB*> k0PiPlusDaughters;
        std::vector<AnaTrueParticleB*> k0PiMinusDaughters;

        // Find all K0 particles
        std::vector<AnaTrueParticleB*> k0Particles;
        for (UInt_t j = 0; j < nTrueParticles; j++) {
            AnaTrueParticleB* truePart = spill->TrueParticles[j];
            if (!truePart) continue;
            if (truePart->PDG == 310) { // K0
                k0Particles.push_back(truePart);
            }
        }

        // Find daughters of each K0 and check for exactly one π+ and one π-
        for (UInt_t k0Idx = 0; k0Idx < k0Particles.size(); k0Idx++) {
            AnaTrueParticleB* k0Particle = k0Particles[k0Idx];

            // Count daughters of this K0
            int piPlusCount = 0;
            int piMinusCount = 0;
            int totalDaughterCount = 0;

            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Check if this particle is a daughter of the K0 (using ParentID)
                if (truePart->ParentID == k0Particle->ID) {
                    totalDaughterCount++;
                    if (truePart->PDG == 211) {
                        piPlusCount++;
                        k0PiPlusDaughters.push_back(truePart);
                    } else if (truePart->PDG == -211) {
                        piMinusCount++;
                        k0PiMinusDaughters.push_back(truePart);
                    }
                }
            }

            // Only keep this K0 if it has exactly 2 daughters: 1 π+ and 1 π-
            if (totalDaughterCount == 2 && piPlusCount == 1 && piMinusCount == 1) {
                hasK0WithExactPiDaughters = true;
                // Print K0 decay reaction
                std::cout << "Event " << i << ": K0 decay found - K0(ID:" << k0Particle->ID << ") -> π+(ID:" << k0PiPlusDaughters[0]->ID << ") + π-(ID:" << k0PiMinusDaughters[0]->ID << ")" << std::endl;
                break; // We only need one qualifying K0 per event
            } else {
                // This K0 doesn't qualify - clear the daughters we found
                k0PiPlusDaughters.clear();
                k0PiMinusDaughters.clear();
            }
        }

        // Only display events with K0 having exactly one π+ and one π- daughter
        if (hasK0WithExactPiDaughters) {
            longPionEventCount++;
            std::cout << "=== Entry " << i << " (K0->π+π- Event #" << longPionEventCount << ") ===" << std::endl;
            AnaEventInfoB* eventInfo = spill->EventInfo;

            // We already have the K0 particle from the selection above
            AnaTrueParticleB* k0Particle = k0Particles[0]; // The qualifying K0

            // First, collect all long π+ and π- particles that are the specific K0 daughters
            std::vector<AnaParticleB*> longPiPlusParticles;
            std::vector<AnaParticleB*> longPiMinusParticles;

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

                    // Check if this is one of our specific K0 daughters and travels > 10cm
                    bool isK0Daughter = false;
                    if (truePart->ID == k0PiPlusDaughters[0]->ID) {
                        isK0Daughter = true;
                    } else if (truePart->ID == k0PiMinusDaughters[0]->ID) {
                        isK0Daughter = true;
                    }

                    if (isK0Daughter) {
                        // Calculate track length
                        float startX = particlePD->PositionStart[0];
                        float startY = particlePD->PositionStart[1];
                        float startZ = particlePD->PositionStart[2];
                        float endX = particlePD->PositionEnd[0];
                        float endY = particlePD->PositionEnd[1];
                        float endZ = particlePD->PositionEnd[2];

                        if (startX > -900 && startY > -900 && startZ > -900 &&
                            endX > -900 && endY > -900 && endZ > -900) {
                            float trackLength = sqrt((endX - startX)*(endX - startX) +
                                                   (endY - startY)*(endY - startY) +
                                                   (endZ - startZ)*(endZ - startZ));

                            if (trackLength > 10.0) { // 10cm threshold
                                if (truePart->PDG == 211) {
                                    longPiPlusParticles.push_back(particle);
                                } else if (truePart->PDG == -211) {
                                    longPiMinusParticles.push_back(particle);
                                }
                            }
                        }
                    }
                }
            }

            // Check if we have both π+ and π- daughters that travel > 10cm
            bool hasLongPiPlus = (longPiPlusParticles.size() > 0);
            bool hasLongPiMinus = (longPiMinusParticles.size() > 0);

            // Only proceed if we have both pions traveling > 10cm
            if (!hasLongPiPlus || !hasLongPiMinus) {
                continue; // Skip this event
            }

            // Create new graphs for this event only
            std::map<int, TGraph*> eventParticleGraphs;
            std::set<AnaParticleB*> selectedParticles; // Track which particles are selected for highlighting
            int particleCount = 0;
            int piPlusCount = 0;
            int piMinusCount = 0;

            // Add the K0 daughter pions to selected particles for highlighting
            for (UInt_t idx = 0; idx < longPiPlusParticles.size(); idx++) {
                selectedParticles.insert(longPiPlusParticles[idx]);
            }
            for (UInt_t idx = 0; idx < longPiMinusParticles.size(); idx++) {
                selectedParticles.insert(longPiMinusParticles[idx]);
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
                    int pointsAdded = 0;
                    for (int plane = 0; plane < 3; plane++) {
                        for (UInt_t h = 0; h < particlePD->Hits[plane].size(); h++) {
                            AnaHitPD hit = particlePD->Hits[plane][h];
                            if (hit.Position.X() > -900 && hit.Position.Y() > -900) {
                                int nPoints = eventParticleGraphs[particleType]->GetN();
                                eventParticleGraphs[particleType]->SetPoint(nPoints, hit.Position.X(), hit.Position.Y());
                                pointsAdded++;
                            }
                        }
                    }

                    // Count particles for verification
                    if (truePart) {
                        if (truePart->PDG == 211) piPlusCount++;
                        if (truePart->PDG == -211) piMinusCount++;
                    }
                }
            }

            // Create canvas title with K0 daughter information
            std::string canvasTitle = "Event " + std::to_string(eventInfo->Event) + " (K^{0}#rightarrow#pi^{+}#pi^{-} Event #" + std::to_string(longPionEventCount) + ")";
            canvasTitle += " - K0 Daughters: " + std::to_string(longPiPlusParticles.size()) + " #pi^{+}, " + std::to_string(longPiMinusParticles.size()) + " #pi^{-}";

            TCanvas *c1 = new TCanvas("c1", canvasTitle.c_str(), 1400, 700);
            c1->Divide(2, 1); // Two subplots side by side

            // Create graphs for both projections
            std::map<int, TGraph*> eventParticleGraphsXY;
            std::map<int, TGraph*> eventParticleGraphsXZ;

            // Special graphs for highlighted closest pair pions
            TGraph* highlightedPiPlusXY = nullptr;
            TGraph* highlightedPiMinusXY = nullptr;
            TGraph* highlightedPiPlusXZ = nullptr;
            TGraph* highlightedPiMinusXZ = nullptr;

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

            // Create highlighted graphs for K0 daughter pions
            if (longPiPlusParticles.size() > 0) {
                highlightedPiPlusXY = new TGraph();
                highlightedPiPlusXY->SetMarkerStyle(24); // Open circle
                highlightedPiPlusXY->SetMarkerSize(1.2);
                highlightedPiPlusXY->SetMarkerColor(kBlue);
                highlightedPiPlusXY->SetLineColor(kBlue);
                highlightedPiPlusXY->SetLineWidth(2);

                highlightedPiPlusXZ = new TGraph();
                highlightedPiPlusXZ->SetMarkerStyle(24);
                highlightedPiPlusXZ->SetMarkerSize(1.2);
                highlightedPiPlusXZ->SetMarkerColor(kBlue);
                highlightedPiPlusXZ->SetLineColor(kBlue);
                highlightedPiPlusXZ->SetLineWidth(2);
            }

            if (longPiMinusParticles.size() > 0) {
                highlightedPiMinusXY = new TGraph();
                highlightedPiMinusXY->SetMarkerStyle(25); // Open square
                highlightedPiMinusXY->SetMarkerSize(1.2);
                highlightedPiMinusXY->SetMarkerColor(kOrange);
                highlightedPiMinusXY->SetLineColor(kOrange);
                highlightedPiMinusXY->SetLineWidth(2);

                highlightedPiMinusXZ = new TGraph();
                highlightedPiMinusXZ->SetMarkerStyle(25);
                highlightedPiMinusXZ->SetMarkerSize(1.2);
                highlightedPiMinusXZ->SetMarkerColor(kOrange);
                highlightedPiMinusXZ->SetLineColor(kOrange);
                highlightedPiMinusXZ->SetLineWidth(2);
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

                                // Add to highlighted graphs if this is a K0 daughter pion
                                if (selectedParticles.find(particle) != selectedParticles.end() && truePart->PDG == 211 && highlightedPiPlusXY) {
                                    int nHighlighted = highlightedPiPlusXY->GetN();
                                    highlightedPiPlusXY->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Y());
                                } else if (selectedParticles.find(particle) != selectedParticles.end() && truePart->PDG == -211 && highlightedPiMinusXY) {
                                    int nHighlighted = highlightedPiMinusXY->GetN();
                                    highlightedPiMinusXY->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Y());
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

                                // Add to highlighted graphs if this is a K0 daughter pion
                                if (selectedParticles.find(particle) != selectedParticles.end() && truePart->PDG == 211 && highlightedPiPlusXZ) {
                                    int nHighlighted = highlightedPiPlusXZ->GetN();
                                    highlightedPiPlusXZ->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Z());
                                } else if (selectedParticles.find(particle) != selectedParticles.end() && truePart->PDG == -211 && highlightedPiMinusXZ) {
                                    int nHighlighted = highlightedPiMinusXZ->GetN();
                                    highlightedPiMinusXZ->SetPoint(nHighlighted, hit.Position.X(), hit.Position.Z());
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

            // Draw highlighted closest pair pions in XY projection
            if (highlightedPiPlusXY && highlightedPiPlusXY->GetN() > 0) {
                highlightedPiPlusXY->Draw("P SAME");
            }
            if (highlightedPiMinusXY && highlightedPiMinusXY->GetN() > 0) {
                highlightedPiMinusXY->Draw("P SAME");
            }

            // Draw K0 trajectory in XY projection
            if (k0Particle) {
                TLine* k0TrajectoryXY = new TLine(k0Particle->Position[0], k0Particle->Position[1],
                                                 k0Particle->PositionEnd[0], k0Particle->PositionEnd[1]);
                k0TrajectoryXY->SetLineColor(kRed);
                k0TrajectoryXY->SetLineWidth(3);
                k0TrajectoryXY->SetLineStyle(2);
                k0TrajectoryXY->Draw();
            }

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

            // Draw highlighted closest pair pions in XZ projection
            if (highlightedPiPlusXZ && highlightedPiPlusXZ->GetN() > 0) {
                highlightedPiPlusXZ->Draw("P SAME");
            }
            if (highlightedPiMinusXZ && highlightedPiMinusXZ->GetN() > 0) {
                highlightedPiMinusXZ->Draw("P SAME");
            }

            // Draw K0 trajectory in XZ projection
            if (k0Particle) {
                TLine* k0TrajectoryXZ = new TLine(k0Particle->Position[0], k0Particle->Position[2],
                                                 k0Particle->PositionEnd[0], k0Particle->PositionEnd[2]);
                k0TrajectoryXZ->SetLineColor(kRed);
                k0TrajectoryXZ->SetLineWidth(3);
                k0TrajectoryXZ->SetLineStyle(2);
                k0TrajectoryXZ->Draw();
            }

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

            // Add highlighted closest pair pions to legend
            if (highlightedPiPlusXY && highlightedPiPlusXY->GetN() > 0) {
                legend->AddEntry(highlightedPiPlusXY, "#pi^{+} (closest pair)", "P");
            }
            if (highlightedPiMinusXY && highlightedPiMinusXY->GetN() > 0) {
                legend->AddEntry(highlightedPiMinusXY, "#pi^{-} (closest pair)", "P");
            }

            // Add K0 trajectory to legend
            if (k0Particle) {
                TLine* k0LegendLine = new TLine(0, 0, 1, 1);
                k0LegendLine->SetLineColor(kRed);
                k0LegendLine->SetLineWidth(3);
                k0LegendLine->SetLineStyle(2);
                legend->AddEntry(k0LegendLine, "K^{0}", "l");
            }

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
            if (highlightedPiPlusXY) delete highlightedPiPlusXY;
            if (highlightedPiMinusXY) delete highlightedPiMinusXY;
            if (highlightedPiPlusXZ) delete highlightedPiPlusXZ;
            if (highlightedPiMinusXZ) delete highlightedPiMinusXZ;

            // Clean up K0 trajectory lines
            if (k0Particle) {
                // Note: TLine objects drawn on canvas are automatically cleaned up by ROOT
                // when the canvas is destroyed, so we don't need to manually delete them
            }

        }

        // Progress output every 1000 events
        if (i % 1000 == 0) {
            std::cout << "Processing entry " << i << "/" << nEntries << " (K0->π+π- events found: " << longPionEventCount << ")" << std::endl;
        }
    }

    // Print final statistics
    std::cout << "\n=== Final Statistics ===" << std::endl;
    std::cout << "Total K0->π+π- events found and displayed: " << longPionEventCount << std::endl;

    f->Close();
}
