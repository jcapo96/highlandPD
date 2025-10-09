// File: K0Statistics.C
// Analyzes all events to find K0 particles that are daughters of K+ beam particles
// Based on truth information
// Modified to process all files in DATA directory and produce combined results

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
// Data classes should be available through ROOT's dictionary system

void K0Statistics(const char* dataDir = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA") {

    // Get list of all .root files in the DATA directory
    TSystemDirectory dir(dataDir, dataDir);
    TList *files = dir.GetListOfFiles();

    if (!files) {
        std::cerr << "Error: Could not access directory " << dataDir << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;
    TIter next(files);
    TSystemFile *file;

    while ((file = (TSystemFile*)next())) {
        std::string filename = file->GetName();
        if (filename.length() > 5 && filename.substr(filename.length() - 5) == ".root") {
            // Exclude the reco2 file as requested
            if (filename.find("reco2") != std::string::npos) {
                std::cout << "Skipping excluded file: " << filename << std::endl;
                continue;
            }
            std::string fullPath = std::string(dataDir) + "/" + filename;
            rootFiles.push_back(fullPath);
        }
    }

    if (rootFiles.empty()) {
        std::cerr << "Error: No .root files found in " << dataDir << std::endl;
        return;
    }

    std::cout << "Found " << rootFiles.size() << " root files to process:" << std::endl;
    for (const auto& file : rootFiles) {
        std::cout << "  " << file << std::endl;
    }
    std::cout << std::endl;

    // Global statistics counters (accumulated across all files)
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

    // Statistics for beam particles (proton and pion beams)
    int eventsWithProtonBeam = 0;
    int eventsWithK0FromProtonBeam = 0;
    int totalK0FromProtonBeam = 0;
    int eventsWithPionBeam = 0;
    int eventsWithK0FromPionBeam = 0;
    int totalK0FromPionBeam = 0;

    // Statistics for K0 from protons (PDG = 2212)
    int eventsWithK0FromProton = 0;
    int totalK0FromProton = 0;
    std::map<int, int> k0DecayModesFromProton; // Decay modes for K0 from protons

    // Statistics for K0 from pions (PDG = 211, -211, 111)
    int eventsWithK0FromPion = 0;
    int totalK0FromPion = 0;
    std::map<int, int> k0DecayModesFromPion; // Decay modes for K0 from pions

    // Combined statistics for all parent types
    int eventsWithK0FromAnyParent = 0;
    int totalK0FromAnyParent = 0;
    std::map<int, int> k0DecayModesFromAnyParent; // Decay modes for K0 from any parent

    // Detailed statistics (accumulated across all files)
    std::map<int, int> k0DecayModes; // PDG -> count
    std::vector<int> k0DecayChains; // Store decay chain information
    std::map<std::string, int> otherDecayModes; // Track specific "other" decay patterns

    // Reconstruction statistics for K0 daughters
    int totalK0Daughters = 0;
    int reconstructedK0Daughters = 0;
    int totalPiPlusDaughters = 0;
    int reconstructedPiPlusDaughters = 0;
    int totalPiMinusDaughters = 0;
    int reconstructedPiMinusDaughters = 0;
    int totalPiZeroDaughters = 0;
    int reconstructedPiZeroDaughters = 0;

    int fileCounter = 0;
    Long64_t globalEventCounter = 0;

    // Process each file
    for (const auto& filename : rootFiles) {
        fileCounter++;
        std::cout << "Processing file " << fileCounter << "/" << rootFiles.size()
                  << ": " << filename << std::endl;

        // Open the file
        TFile *f = TFile::Open(filename.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file " << filename << std::endl;
            continue;
        }

        // Get the tree
        TTree *MiniTree = (TTree*)f->Get("MiniTree");
        if (!MiniTree) {
            std::cerr << "Could not find MiniTree in file " << filename << std::endl;
            f->Close();
            delete f;
            continue;
        }

        // Set up the object for the branch
        AnaSpill* spill = nullptr;
        MiniTree->SetBranchAddress("Spill", &spill);

        Long64_t nEntries = MiniTree->GetEntries();
        std::cout << "  Events in this file: " << nEntries << std::endl;

        // File-specific counters (will be added to global counters)
        int fileEventsWithKPlusBeam = 0;
        int fileEventsWithK0FromKPlusBeam = 0;
        int fileTotalK0FromKPlusBeam = 0;
        int fileEventsWithK0FromKPlusNonBeam = 0;
        int fileTotalK0FromKPlusNonBeam = 0;
        int fileTotalK0WithPiDaughters = 0;
        int fileEventsWithK0PiDaughters = 0;
        int fileTotalK0FromAnyKPlus = 0;
        int fileEventsWithK0FromAnyKPlus = 0;

        // File-specific counters for new beam types
        int fileEventsWithProtonBeam = 0;
        int fileEventsWithK0FromProtonBeam = 0;
        int fileTotalK0FromProtonBeam = 0;
        int fileEventsWithPionBeam = 0;
        int fileEventsWithK0FromPionBeam = 0;
        int fileTotalK0FromPionBeam = 0;

        // File-specific counters for new parent types
        int fileEventsWithK0FromProton = 0;
        int fileTotalK0FromProton = 0;
        int fileEventsWithK0FromPion = 0;
        int fileTotalK0FromPion = 0;
        int fileEventsWithK0FromAnyParent = 0;
        int fileTotalK0FromAnyParent = 0;

        // File-specific reconstruction statistics
        int fileTotalK0Daughters = 0;
        int fileReconstructedK0Daughters = 0;
        int fileTotalPiPlusDaughters = 0;
        int fileReconstructedPiPlusDaughters = 0;
        int fileTotalPiMinusDaughters = 0;
        int fileReconstructedPiMinusDaughters = 0;
        int fileTotalPiZeroDaughters = 0;
        int fileReconstructedPiZeroDaughters = 0;

        for (Long64_t i = 0; i < nEntries; i++) {
            MiniTree->GetEntry(i);
            if (!spill) continue;

            totalEvents++;
            globalEventCounter++;

            // Check if there's a beam particle
            AnaParticleB* beamParticle = nullptr;
            if (spill->Beam) {
                AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
                if (beam && beam->BeamParticle) {
                    beamParticle = beam->BeamParticle;
                }
            }

            // Check beam particle type (K+, proton, or pion)
            bool hasKPlusBeam = false;
            bool hasProtonBeam = false;
            bool hasPionBeam = false;
            AnaTrueParticleB* beamTrue = nullptr;
            if (beamParticle) {
                beamTrue = beamParticle->GetTrueParticle();
                if (beamTrue) {
                    if (beamTrue->PDG == 321) { // K+
                        hasKPlusBeam = true;
                        eventsWithKPlusBeam++;
                        fileEventsWithKPlusBeam++;
                    } else if (beamTrue->PDG == 2212) { // Proton
                        hasProtonBeam = true;
                        eventsWithProtonBeam++;
                        fileEventsWithProtonBeam++;
                    } else if (beamTrue->PDG == 211 || beamTrue->PDG == -211 || beamTrue->PDG == 111) { // Pions
                        hasPionBeam = true;
                        eventsWithPionBeam++;
                        fileEventsWithPionBeam++;
                    }
                }
            }

            UInt_t nTrueParticles = spill->TrueParticles.size();

            // Find all potential parent particles in the event
            std::vector<AnaTrueParticleB*> kPlusParticles;
            std::vector<AnaTrueParticleB*> protonParticles;
            std::vector<AnaTrueParticleB*> pionParticles;

            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                if (truePart->PDG == 321) { // K+
                    kPlusParticles.push_back(truePart);
                } else if (truePart->PDG == 2212) { // Proton
                    protonParticles.push_back(truePart);
                } else if (truePart->PDG == 211 || truePart->PDG == -211 || truePart->PDG == 111) { // Pions
                    pionParticles.push_back(truePart);
                }
            }

            // Find all K0 particles and categorize by parent type
            std::vector<AnaTrueParticleB*> k0FromKPlusBeam;
            std::vector<AnaTrueParticleB*> k0FromKPlusNonBeam;
            std::vector<AnaTrueParticleB*> k0FromAnyKPlus;
            std::vector<AnaTrueParticleB*> k0FromProtonBeam;
            std::vector<AnaTrueParticleB*> k0FromPionBeam;
            std::vector<AnaTrueParticleB*> k0FromProton;
            std::vector<AnaTrueParticleB*> k0FromPion;
            std::vector<AnaTrueParticleB*> k0FromAnyParent;

            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Check if this is a K0 (PDG = 310)
                if (truePart->PDG == 310) {
                    bool foundParent = false;

                    // Check if this K0 is a daughter of any K+ particle
                    for (UInt_t k = 0; k < kPlusParticles.size(); k++) {
                        AnaTrueParticleB* kPlus = kPlusParticles[k];
                        if (truePart->ParentID == kPlus->ID) {
                            k0FromAnyKPlus.push_back(truePart);
                            k0FromAnyParent.push_back(truePart);
                            totalK0FromAnyKPlus++;
                            fileTotalK0FromAnyKPlus++;
                            totalK0FromAnyParent++;
                            fileTotalK0FromAnyParent++;

                            // Check if this K+ is the beam particle
                            if (beamTrue && kPlus->ID == beamTrue->ID) {
                                k0FromKPlusBeam.push_back(truePart);
                                totalK0FromKPlusBeam++;
                                fileTotalK0FromKPlusBeam++;
                            } else {
                                k0FromKPlusNonBeam.push_back(truePart);
                                totalK0FromKPlusNonBeam++;
                                fileTotalK0FromKPlusNonBeam++;
                            }
                            foundParent = true;
                            break; // Found the parent, no need to check other K+ particles
                        }
                    }

                    // Check if this K0 is a daughter of any proton
                    if (!foundParent) {
                        for (UInt_t k = 0; k < protonParticles.size(); k++) {
                            AnaTrueParticleB* proton = protonParticles[k];
                            if (truePart->ParentID == proton->ID) {
                                k0FromProton.push_back(truePart);
                                k0FromAnyParent.push_back(truePart);
                                totalK0FromProton++;
                                fileTotalK0FromProton++;
                                totalK0FromAnyParent++;
                                fileTotalK0FromAnyParent++;

                                // Check if this proton is the beam particle
                                if (beamTrue && proton->ID == beamTrue->ID) {
                                    k0FromProtonBeam.push_back(truePart);
                                    totalK0FromProtonBeam++;
                                    fileTotalK0FromProtonBeam++;
                                }
                                foundParent = true;
                                break;
                            }
                        }
                    }

                    // Check if this K0 is a daughter of any pion
                    if (!foundParent) {
                        for (UInt_t k = 0; k < pionParticles.size(); k++) {
                            AnaTrueParticleB* pion = pionParticles[k];
                            if (truePart->ParentID == pion->ID) {
                                k0FromPion.push_back(truePart);
                                k0FromAnyParent.push_back(truePart);
                                totalK0FromPion++;
                                fileTotalK0FromPion++;
                                totalK0FromAnyParent++;
                                fileTotalK0FromAnyParent++;

                                // Check if this pion is the beam particle
                                if (beamTrue && pion->ID == beamTrue->ID) {
                                    k0FromPionBeam.push_back(truePart);
                                    totalK0FromPionBeam++;
                                    fileTotalK0FromPionBeam++;
                                }
                                foundParent = true;
                                break;
                            }
                        }
                    }
                }
            }

            // Update event counters
            if (k0FromKPlusBeam.size() > 0) {
                eventsWithK0FromKPlusBeam++;
                fileEventsWithK0FromKPlusBeam++;
            }
            if (k0FromKPlusNonBeam.size() > 0) {
                eventsWithK0FromKPlusNonBeam++;
                fileEventsWithK0FromKPlusNonBeam++;
            }
            if (k0FromAnyKPlus.size() > 0) {
                eventsWithK0FromAnyKPlus++;
                fileEventsWithK0FromAnyKPlus++;
            }
            if (k0FromProtonBeam.size() > 0) {
                eventsWithK0FromProtonBeam++;
                fileEventsWithK0FromProtonBeam++;
            }
            if (k0FromPionBeam.size() > 0) {
                eventsWithK0FromPionBeam++;
                fileEventsWithK0FromPionBeam++;
            }
            if (k0FromProton.size() > 0) {
                eventsWithK0FromProton++;
                fileEventsWithK0FromProton++;
            }
            if (k0FromPion.size() > 0) {
                eventsWithK0FromPion++;
                fileEventsWithK0FromPion++;
            }
            if (k0FromAnyParent.size() > 0) {
                eventsWithK0FromAnyParent++;
                fileEventsWithK0FromAnyParent++;
            }

            // Analyze decays of all K0 particles from any parent
            for (UInt_t k0Idx = 0; k0Idx < k0FromAnyParent.size(); k0Idx++) {
                AnaTrueParticleB* k0Particle = k0FromAnyParent[k0Idx];

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

                        // Count total daughters (reconstruction counting moved to separate loop)
                        totalK0Daughters++;
                        fileTotalK0Daughters++;

                        // Count by particle type
                        if (pdg == 211) {
                            piPlusCount++;
                            totalPiPlusDaughters++;
                            fileTotalPiPlusDaughters++;
                        }
                        else if (pdg == -211) {
                            piMinusCount++;
                            totalPiMinusDaughters++;
                            fileTotalPiMinusDaughters++;
                        }
                        else if (pdg == 111) {
                            piZeroCount++;
                            totalPiZeroDaughters++;
                            fileTotalPiZeroDaughters++;
                        }
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

                        // Create a string representation of the decay pattern
                        std::string decayPattern = "";
                        if (piPlusCount > 0) decayPattern += std::to_string(piPlusCount) + "π+";
                        if (piMinusCount > 0) {
                            if (!decayPattern.empty()) decayPattern += "+";
                            decayPattern += std::to_string(piMinusCount) + "π-";
                        }
                        if (piZeroCount > 0) {
                            if (!decayPattern.empty()) decayPattern += "+";
                            decayPattern += std::to_string(piZeroCount) + "π0";
                        }
                        if (electronCount > 0) {
                            if (!decayPattern.empty()) decayPattern += "+";
                            decayPattern += std::to_string(electronCount) + "e";
                        }
                        if (muonCount > 0) {
                            if (!decayPattern.empty()) decayPattern += "+";
                            decayPattern += std::to_string(muonCount) + "μ";
                        }
                        if (otherCount > 0) {
                            if (!decayPattern.empty()) decayPattern += "+";
                            decayPattern += std::to_string(otherCount) + "other";
                        }
                        otherDecayModes[decayPattern]++;
                    }

                    k0DecayModes[decayMode]++;
                    k0DecayModesFromAnyParent[decayMode]++;

                    // Track decay modes by parent type
                    // Check if this K0 comes from K+, proton, or pion
                    bool isFromKPlus = false;
                    bool isFromProton = false;
                    bool isFromPion = false;

                    // Check parent type
                    for (UInt_t k = 0; k < kPlusParticles.size(); k++) {
                        if (k0Particle->ParentID == kPlusParticles[k]->ID) {
                            isFromKPlus = true;
                            break;
                        }
                    }
                    if (!isFromKPlus) {
                        for (UInt_t k = 0; k < protonParticles.size(); k++) {
                            if (k0Particle->ParentID == protonParticles[k]->ID) {
                                isFromProton = true;
                                break;
                            }
                        }
                    }
                    if (!isFromKPlus && !isFromProton) {
                        for (UInt_t k = 0; k < pionParticles.size(); k++) {
                            if (k0Particle->ParentID == pionParticles[k]->ID) {
                                isFromPion = true;
                                break;
                            }
                        }
                    }

                    // Record decay mode by parent type
                    if (isFromKPlus) {
                        // Already handled by original K+ analysis
                    } else if (isFromProton) {
                        k0DecayModesFromProton[decayMode]++;
                    } else if (isFromPion) {
                        k0DecayModesFromPion[decayMode]++;
                    }

                    // Store decay chain information
                    if (decayMode == 1) { // π+π- decay
                        totalK0WithPiDaughters++;
                        fileTotalK0WithPiDaughters++;
                        if (eventsWithK0PiDaughters == 0 ||
                            std::find(k0DecayChains.begin(), k0DecayChains.end(), globalEventCounter) == k0DecayChains.end()) {
                            eventsWithK0PiDaughters++;
                            fileEventsWithK0PiDaughters++;
                            k0DecayChains.push_back(globalEventCounter);
                        }
                    }

                    // Print detailed decay information for first few events across all files
                    if (globalEventCounter <= 10) {
                        std::cout << "Global Event " << globalEventCounter << " (File " << fileCounter << ", Local " << i
                                  << "): K+ -> K0(ID:" << k0Particle->ID << ") -> ";
                        for (UInt_t dIdx = 0; dIdx < k0Daughters.size(); dIdx++) {
                            AnaTrueParticleB* daughter = k0Daughters[dIdx];
                            if (dIdx > 0) std::cout << " + ";
                            std::cout << "PDG:" << daughter->PDG;
                        }
                        std::cout << " (Mode: " << decayMode << ")" << std::endl;
                    }
                }
            }

            // Now analyze reconstructed particles to find K0 daughters
            // First, collect all unique K0 daughter true particle IDs that we found
            std::set<Int_t> k0DaughterTrueIDs;
            std::set<Int_t> allK0IDs;

            // Debug: Find all K0s first
            for (UInt_t k = 0; k < nTrueParticles; k++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[k];
                if (!truePart) continue;
                if (truePart->PDG == 310) {
                    allK0IDs.insert(truePart->ID);
                }
            }

            // Debug: Find K0 daughters
            for (UInt_t k = 0; k < nTrueParticles; k++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[k];
                if (!truePart) continue;

                // Check if this true particle is a daughter of a K0
                for (UInt_t l = 0; l < nTrueParticles; l++) {
                    AnaTrueParticleB* candidate = spill->TrueParticles[l];
                    if (!candidate) continue;

                    // Check if this candidate is a K0 and is the parent of our true particle
                    if (candidate->PDG == 310 && candidate->ID == truePart->ParentID) {
                        k0DaughterTrueIDs.insert(truePart->ID);
                        break;
                    }
                }
            }

            // Debug: Show detailed info for first few events
            if (i < 3) {
                std::cout << "Debug Event " << i << ": spill->Bunches.size() = " << spill->Bunches.size()
                          << ", Total K0s found: " << allK0IDs.size()
                          << ", K0 daughter true IDs found: " << k0DaughterTrueIDs.size() << std::endl;

                // Debug: Show some K0 daughter details
                if (k0DaughterTrueIDs.size() > 0) {
                    std::cout << "  K0 daughter true IDs: ";
                    int count = 0;
                    for (auto id : k0DaughterTrueIDs) {
                        if (count < 5) {
                            std::cout << id << " ";
                            count++;
                        }
                    }
                    std::cout << std::endl;
                }
            }

            // Now look for reconstructed particles that correspond to these K0 daughters
            std::set<Int_t> alreadyCountedReconParticles; // To avoid double counting across bunches
            std::map<Int_t, Int_t> trueIDToReconCount; // Track how many reconstructed particles point to each true particle

            for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
                AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
                UInt_t nReconParticles = bunch->Particles.size();

                for (UInt_t j = 0; j < nReconParticles; j++) {
                    AnaParticleB* reconPart = bunch->Particles[j];
                    if (!reconPart) continue;

                    // Check if we already counted this reconstructed particle
                    if (alreadyCountedReconParticles.find(reconPart->UniqueID) != alreadyCountedReconParticles.end()) {
                        continue; // Already counted this reconstructed particle
                    }

                    // Check if this reconstructed particle has an associated true particle
                    if (!reconPart->TrueObject) continue;

                    // Cast to AnaTrueParticleB to access the true particle info
                    AnaTrueParticleB* truePart = static_cast<AnaTrueParticleB*>(reconPart->TrueObject);
                    if (!truePart) continue;

                    // Get PDG for this true particle
                    int pdg = truePart->PDG;

                    // Debug: Track all reconstructed particles with true objects for first few events
                    if (i < 3 && reconstructedK0Daughters < 10) {
                        std::cout << "  Recon particle: UniqueID=" << reconPart->UniqueID
                                  << ", TrueID=" << truePart->ID << ", PDG=" << pdg
                                  << ", ParentID=" << truePart->ParentID << std::endl;
                    }

                    // Track potential particle splitting
                    trueIDToReconCount[truePart->ID]++;

                    // Check if this true particle is one of the K0 daughters we identified
                    if (k0DaughterTrueIDs.find(truePart->ID) != k0DaughterTrueIDs.end()) {
                        // This reconstructed particle corresponds to a K0 daughter

                        // Mark this reconstructed particle as counted
                        alreadyCountedReconParticles.insert(reconPart->UniqueID);

                        // Count total reconstructed K0 daughters
                        reconstructedK0Daughters++;
                        fileReconstructedK0Daughters++;

                        // Count by particle type
                        if (pdg == 211) {
                            reconstructedPiPlusDaughters++;
                            fileReconstructedPiPlusDaughters++;
                        }
                        else if (pdg == -211) {
                            reconstructedPiMinusDaughters++;
                            fileReconstructedPiMinusDaughters++;
                        }
                        else if (pdg == 111) {
                            reconstructedPiZeroDaughters++;
                            fileReconstructedPiZeroDaughters++;
                        }

                        // Debug output for first few reconstructed daughters
                        if (reconstructedK0Daughters <= 5) {
                            std::cout << "Reconstructed K0 daughter " << reconstructedK0Daughters
                                      << " PDG=" << pdg << " TrueID=" << truePart->ID
                                      << " ReconUniqueID=" << reconPart->UniqueID << std::endl;
                        }
                    }
                    // Debug: Check if we're finding reconstructed particles that should be K0 daughters but aren't in our set
                    else if (i < 3 && (pdg == 211 || pdg == -211 || pdg == 111)) {
                        // This is a pion/pion0 that might be a K0 daughter but we didn't identify it
                        bool foundParent = false;
                        for (UInt_t m = 0; m < nTrueParticles; m++) {
                            AnaTrueParticleB* parent = spill->TrueParticles[m];
                            if (!parent) continue;
                            if (parent->PDG == 310 && parent->ID == truePart->ParentID) {
                                foundParent = true;
                                break;
                            }
                        }
                        if (foundParent) {
                            std::cout << "  WARNING: Found reconstructed " << pdg
                                      << " with K0 parent (TrueID=" << truePart->ID
                                      << ", ParentID=" << truePart->ParentID
                                      << ") but not in k0DaughterTrueIDs set!" << std::endl;
                        }
                    }
                }
            }

            // Debug: Check for particle splitting in first few events
            if (i < 3) {
                std::cout << "  Particle splitting analysis:" << std::endl;
                for (auto& pair : trueIDToReconCount) {
                    if (pair.second > 1) {
                        std::cout << "    TrueID " << pair.first << " has " << pair.second
                                  << " reconstructed particles (SPLITTING DETECTED!)" << std::endl;
                    }
                }
                std::cout << "  Total reconstructed particles: " << trueIDToReconCount.size() << std::endl;
                std::cout << "  Total reconstructed K0 daughters found: " << reconstructedK0Daughters << std::endl;
            }

            // Progress indicator for large files
            if (i % 10000 == 0 && i > 0) {
                std::cout << "    File progress: " << i << "/" << nEntries << std::endl;
            }
        }

        // Print file summary
        std::cout << "  File " << fileCounter << " summary:" << std::endl;
        std::cout << "    Events with K+ beam: " << fileEventsWithKPlusBeam << std::endl;
        std::cout << "    Events with K0 from K+ beam: " << fileEventsWithK0FromKPlusBeam << std::endl;
        std::cout << "    Total K0 from K+ beam: " << fileTotalK0FromKPlusBeam << std::endl;
        std::cout << "    Events with K0 from K+ non-beam: " << fileEventsWithK0FromKPlusNonBeam << std::endl;
        std::cout << "    Total K0 from K+ non-beam: " << fileTotalK0FromKPlusNonBeam << std::endl;
        std::cout << "    Total K0 from any K+: " << fileTotalK0FromAnyKPlus << std::endl;
        std::cout << "    Events with proton beam: " << fileEventsWithProtonBeam << std::endl;
        std::cout << "    Events with K0 from proton beam: " << fileEventsWithK0FromProtonBeam << std::endl;
        std::cout << "    Total K0 from proton beam: " << fileTotalK0FromProtonBeam << std::endl;
        std::cout << "    Events with pion beam: " << fileEventsWithPionBeam << std::endl;
        std::cout << "    Events with K0 from pion beam: " << fileEventsWithK0FromPionBeam << std::endl;
        std::cout << "    Total K0 from pion beam: " << fileTotalK0FromPionBeam << std::endl;
        std::cout << "    Events with K0 from protons: " << fileEventsWithK0FromProton << std::endl;
        std::cout << "    Total K0 from protons: " << fileTotalK0FromProton << std::endl;
        std::cout << "    Events with K0 from pions: " << fileEventsWithK0FromPion << std::endl;
        std::cout << "    Total K0 from pions: " << fileTotalK0FromPion << std::endl;
        std::cout << "    Events with K0 from any parent: " << fileEventsWithK0FromAnyParent << std::endl;
        std::cout << "    Total K0 from any parent: " << fileTotalK0FromAnyParent << std::endl;
        std::cout << "    K0 with π+π- daughters: " << fileTotalK0WithPiDaughters << std::endl;
        std::cout << "    Events with K0→π+π- decays: " << fileEventsWithK0PiDaughters << std::endl;
        std::cout << "    Total K0 daughters: " << fileTotalK0Daughters << std::endl;
        std::cout << "    Reconstructed K0 daughters: " << fileReconstructedK0Daughters << std::endl;
        std::cout << "    K0 daughter reconstruction efficiency: " << (fileTotalK0Daughters > 0 ? (double)fileReconstructedK0Daughters / fileTotalK0Daughters * 100.0 : 0.0) << "%" << std::endl;
        std::cout << "    Total π+ daughters: " << fileTotalPiPlusDaughters << std::endl;
        std::cout << "    Reconstructed π+ daughters: " << fileReconstructedPiPlusDaughters << std::endl;
        std::cout << "    π+ reconstruction efficiency: " << (fileTotalPiPlusDaughters > 0 ? (double)fileReconstructedPiPlusDaughters / fileTotalPiPlusDaughters * 100.0 : 0.0) << "%" << std::endl;
        std::cout << "    Total π- daughters: " << fileTotalPiMinusDaughters << std::endl;
        std::cout << "    Reconstructed π- daughters: " << fileReconstructedPiMinusDaughters << std::endl;
        std::cout << "    π- reconstruction efficiency: " << (fileTotalPiMinusDaughters > 0 ? (double)fileReconstructedPiMinusDaughters / fileTotalPiMinusDaughters * 100.0 : 0.0) << "%" << std::endl;
        std::cout << "    Total π0 daughters: " << fileTotalPiZeroDaughters << std::endl;
        std::cout << "    Reconstructed π0 daughters: " << fileReconstructedPiZeroDaughters << std::endl;
        std::cout << "    π0 reconstruction efficiency: " << (fileTotalPiZeroDaughters > 0 ? (double)fileReconstructedPiZeroDaughters / fileTotalPiZeroDaughters * 100.0 : 0.0) << "%" << std::endl;
        std::cout << std::endl;

        // Close the file
        f->Close();
        delete f;
    }

    // Print final combined statistics
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "COMBINED K0 STATISTICS ANALYSIS (ALL FILES)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Files processed: " << rootFiles.size() << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << "Events with K+ beam: " << eventsWithKPlusBeam << std::endl;
    std::cout << "Events with K0 daughters of K+ beam: " << eventsWithK0FromKPlusBeam << std::endl;
    std::cout << "Events with K0 daughters of K+ non-beam: " << eventsWithK0FromKPlusNonBeam << std::endl;
    std::cout << "Events with K0 daughters of any K+: " << eventsWithK0FromAnyKPlus << std::endl;
    std::cout << "Events with proton beam: " << eventsWithProtonBeam << std::endl;
    std::cout << "Events with K0 daughters of proton beam: " << eventsWithK0FromProtonBeam << std::endl;
    std::cout << "Events with pion beam: " << eventsWithPionBeam << std::endl;
    std::cout << "Events with K0 daughters of pion beam: " << eventsWithK0FromPionBeam << std::endl;
    std::cout << "Events with K0 daughters of protons: " << eventsWithK0FromProton << std::endl;
    std::cout << "Events with K0 daughters of pions: " << eventsWithK0FromPion << std::endl;
    std::cout << "Events with K0 daughters of any parent: " << eventsWithK0FromAnyParent << std::endl;
    std::cout << "Total K0 particles from K+ beam: " << totalK0FromKPlusBeam << std::endl;
    std::cout << "Total K0 particles from K+ non-beam: " << totalK0FromKPlusNonBeam << std::endl;
    std::cout << "Total K0 particles from any K+: " << totalK0FromAnyKPlus << std::endl;
    std::cout << "Total K0 particles from proton beam: " << totalK0FromProtonBeam << std::endl;
    std::cout << "Total K0 particles from pion beam: " << totalK0FromPionBeam << std::endl;
    std::cout << "Total K0 particles from protons: " << totalK0FromProton << std::endl;
    std::cout << "Total K0 particles from pions: " << totalK0FromPion << std::endl;
    std::cout << "Total K0 particles from any parent: " << totalK0FromAnyParent << std::endl;
    std::cout << "K0 particles with π+π- daughters: " << totalK0WithPiDaughters << std::endl;
    std::cout << "Events with K0→π+π- decays: " << eventsWithK0PiDaughters << std::endl;

    std::cout << "\n=== K0 Decay Modes (Combined) ===" << std::endl;
    std::cout << "Mode 1 (π+π-): " << k0DecayModes[1] << std::endl;
    std::cout << "Mode 2 (π0π0): " << k0DecayModes[2] << std::endl;
    std::cout << "Mode 3 (π0e+νe): " << k0DecayModes[3] << std::endl;
    std::cout << "Mode 4 (π0μ+νμ): " << k0DecayModes[4] << std::endl;
    std::cout << "Mode 5 (e+e-): " << k0DecayModes[5] << std::endl;
    std::cout << "Mode 6 (μ+μ-): " << k0DecayModes[6] << std::endl;
    std::cout << "Mode 99 (Other): " << k0DecayModes[99] << std::endl;

    // Print decay modes by parent type
    if (totalK0FromProton > 0) {
        std::cout << "\n=== K0 Decay Modes from Protons ===" << std::endl;
        std::cout << "Mode 1 (π+π-): " << k0DecayModesFromProton[1] << std::endl;
        std::cout << "Mode 2 (π0π0): " << k0DecayModesFromProton[2] << std::endl;
        std::cout << "Mode 3 (π0e+νe): " << k0DecayModesFromProton[3] << std::endl;
        std::cout << "Mode 4 (π0μ+νμ): " << k0DecayModesFromProton[4] << std::endl;
        std::cout << "Mode 5 (e+e-): " << k0DecayModesFromProton[5] << std::endl;
        std::cout << "Mode 6 (μ+μ-): " << k0DecayModesFromProton[6] << std::endl;
        std::cout << "Mode 99 (Other): " << k0DecayModesFromProton[99] << std::endl;
    }

    if (totalK0FromPion > 0) {
        std::cout << "\n=== K0 Decay Modes from Pions ===" << std::endl;
        std::cout << "Mode 1 (π+π-): " << k0DecayModesFromPion[1] << std::endl;
        std::cout << "Mode 2 (π0π0): " << k0DecayModesFromPion[2] << std::endl;
        std::cout << "Mode 3 (π0e+νe): " << k0DecayModesFromPion[3] << std::endl;
        std::cout << "Mode 4 (π0μ+νμ): " << k0DecayModesFromPion[4] << std::endl;
        std::cout << "Mode 5 (e+e-): " << k0DecayModesFromPion[5] << std::endl;
        std::cout << "Mode 6 (μ+μ-): " << k0DecayModesFromPion[6] << std::endl;
        std::cout << "Mode 99 (Other): " << k0DecayModesFromPion[99] << std::endl;
    }

    // Print detailed breakdown of "Other" decay modes
    if (k0DecayModes[99] > 0) {
        std::cout << "\n=== Detailed 'Other' Decay Modes ===" << std::endl;
        std::cout << "Total 'Other' decays: " << k0DecayModes[99] << std::endl;

        // Sort by count (most frequent first)
        std::vector<std::pair<std::string, int>> sortedOtherModes;
        for (const auto& pair : otherDecayModes) {
            sortedOtherModes.push_back(pair);
        }
        std::sort(sortedOtherModes.begin(), sortedOtherModes.end(),
                 [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
                     return a.second > b.second;
                 });

        std::cout << "Top 'Other' decay patterns:" << std::endl;
        int count = 0;
        for (const auto& pair : sortedOtherModes) {
            if (count >= 20) { // Limit to top 20
                std::cout << "... and " << (sortedOtherModes.size() - 20) << " more patterns" << std::endl;
                break;
            }
            std::cout << "  " << pair.first << ": " << pair.second << " decays ("
                      << (double)pair.second / k0DecayModes[99] * 100.0 << "%)" << std::endl;
            count++;
        }
    }

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

    if (totalK0FromProton > 0) {
        std::cout << "\n=== Decay Fractions (All K0 from Proton parents) ===" << std::endl;
        std::cout << "π+π- fraction: " << (double)k0DecayModesFromProton[1] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "π0π0 fraction: " << (double)k0DecayModesFromProton[2] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "π0e+νe fraction: " << (double)k0DecayModesFromProton[3] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "π0μ+νμ fraction: " << (double)k0DecayModesFromProton[4] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "e+e- fraction: " << (double)k0DecayModesFromProton[5] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "μ+μ- fraction: " << (double)k0DecayModesFromProton[6] / totalK0FromProton * 100.0 << "%" << std::endl;
        std::cout << "Other fraction: " << (double)k0DecayModesFromProton[99] / totalK0FromProton * 100.0 << "%" << std::endl;
    }

    if (totalK0FromPion > 0) {
        std::cout << "\n=== Decay Fractions (All K0 from Pion parents) ===" << std::endl;
        std::cout << "π+π- fraction: " << (double)k0DecayModesFromPion[1] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "π0π0 fraction: " << (double)k0DecayModesFromPion[2] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "π0e+νe fraction: " << (double)k0DecayModesFromPion[3] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "π0μ+νμ fraction: " << (double)k0DecayModesFromPion[4] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "e+e- fraction: " << (double)k0DecayModesFromPion[5] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "μ+μ- fraction: " << (double)k0DecayModesFromPion[6] / totalK0FromPion * 100.0 << "%" << std::endl;
        std::cout << "Other fraction: " << (double)k0DecayModesFromPion[99] / totalK0FromPion * 100.0 << "%" << std::endl;
    }

    if (totalK0FromAnyParent > 0) {
        std::cout << "\n=== Decay Fractions (All K0 from Any Parent) ===" << std::endl;
        std::cout << "π+π- fraction: " << (double)k0DecayModesFromAnyParent[1] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "π0π0 fraction: " << (double)k0DecayModesFromAnyParent[2] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "π0e+νe fraction: " << (double)k0DecayModesFromAnyParent[3] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "π0μ+νμ fraction: " << (double)k0DecayModesFromAnyParent[4] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "e+e- fraction: " << (double)k0DecayModesFromAnyParent[5] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "μ+μ- fraction: " << (double)k0DecayModesFromAnyParent[6] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
        std::cout << "Other fraction: " << (double)k0DecayModesFromAnyParent[99] / totalK0FromAnyParent * 100.0 << "%" << std::endl;
    }

    std::cout << "\n=== Beam Interaction Fractions ===" << std::endl;
    if (eventsWithKPlusBeam > 0) {
        std::cout << "K+ beam events with K0 production: " << (double)eventsWithK0FromKPlusBeam / eventsWithKPlusBeam * 100.0 << "%" << std::endl;
        std::cout << "K+ beam events with K0→π+π-: " << (double)eventsWithK0PiDaughters / eventsWithKPlusBeam * 100.0 << "%" << std::endl;
    }
    if (eventsWithProtonBeam > 0) {
        std::cout << "Proton beam events with K0 production: " << (double)eventsWithK0FromProtonBeam / eventsWithProtonBeam * 100.0 << "%" << std::endl;
    }
    if (eventsWithPionBeam > 0) {
        std::cout << "Pion beam events with K0 production: " << (double)eventsWithK0FromPionBeam / eventsWithPionBeam * 100.0 << "%" << std::endl;
    }

    std::cout << "\n=== Final Summary ===" << std::endl;
    std::cout << "Processed " << rootFiles.size() << " files with " << totalEvents << " total events" << std::endl;
    std::cout << "Beam particle statistics:" << std::endl;
    std::cout << "  K+ beam events: " << eventsWithKPlusBeam << std::endl;
    std::cout << "  Proton beam events: " << eventsWithProtonBeam << std::endl;
    std::cout << "  Pion beam events: " << eventsWithPionBeam << std::endl;
    std::cout << "K0 production statistics:" << std::endl;
    std::cout << "  Found " << totalK0FromKPlusBeam << " K0 particles from K+ beam interactions" << std::endl;
    std::cout << "  Found " << totalK0FromProtonBeam << " K0 particles from proton beam interactions" << std::endl;
    std::cout << "  Found " << totalK0FromPionBeam << " K0 particles from pion beam interactions" << std::endl;
    std::cout << "  Found " << totalK0FromKPlusNonBeam << " K0 particles from K+ non-beam interactions" << std::endl;
    std::cout << "  Found " << totalK0FromProton << " K0 particles from proton interactions" << std::endl;
    std::cout << "  Found " << totalK0FromPion << " K0 particles from pion interactions" << std::endl;
    std::cout << "  Total K0 particles from any K+ parent: " << totalK0FromAnyKPlus << std::endl;
    std::cout << "  Total K0 particles from any parent: " << totalK0FromAnyParent << std::endl;
    std::cout << "Decay statistics:" << std::endl;
    std::cout << "  Of all K0 from K+ parents, " << totalK0WithPiDaughters << " decayed to π+π- ("
              << (totalK0FromAnyKPlus > 0 ? (double)totalK0WithPiDaughters / totalK0FromAnyKPlus * 100.0 : 0.0) << "%)" << std::endl;
    std::cout << "  Of all K0 from any parent, " << k0DecayModesFromAnyParent[1] << " decayed to π+π- ("
              << (totalK0FromAnyParent > 0 ? (double)k0DecayModesFromAnyParent[1] / totalK0FromAnyParent * 100.0 : 0.0) << "%)" << std::endl;

    std::cout << "\nReconstruction statistics for K0 daughters:" << std::endl;
    std::cout << "  Total K0 daughters: " << totalK0Daughters << std::endl;
    std::cout << "  Reconstructed K0 daughters: " << reconstructedK0Daughters << std::endl;
    std::cout << "  Overall K0 daughter reconstruction efficiency: " << (totalK0Daughters > 0 ? (double)reconstructedK0Daughters / totalK0Daughters * 100.0 : 0.0) << "%" << std::endl;
    std::cout << "  Total π+ daughters: " << totalPiPlusDaughters << std::endl;
    std::cout << "  Reconstructed π+ daughters: " << reconstructedPiPlusDaughters << std::endl;
    std::cout << "  π+ reconstruction efficiency: " << (totalPiPlusDaughters > 0 ? (double)reconstructedPiPlusDaughters / totalPiPlusDaughters * 100.0 : 0.0) << "%" << std::endl;
    std::cout << "  Total π- daughters: " << totalPiMinusDaughters << std::endl;
    std::cout << "  Reconstructed π- daughters: " << reconstructedPiMinusDaughters << std::endl;
    std::cout << "  π- reconstruction efficiency: " << (totalPiMinusDaughters > 0 ? (double)reconstructedPiMinusDaughters / totalPiMinusDaughters * 100.0 : 0.0) << "%" << std::endl;
    std::cout << "  Total π0 daughters: " << totalPiZeroDaughters << std::endl;
    std::cout << "  Reconstructed π0 daughters: " << reconstructedPiZeroDaughters << std::endl;
    std::cout << "  π0 reconstruction efficiency: " << (totalPiZeroDaughters > 0 ? (double)reconstructedPiZeroDaughters / totalPiZeroDaughters * 100.0 : 0.0) << "%" << std::endl;

    std::cout << "\nAnalysis complete!" << std::endl;
}