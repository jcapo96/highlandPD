// File: K0ReconstructionAnalysis.C
// Analyzes K0 particles and their reconstruction efficiency
// Focuses on K0 from K+ (beam and non-beam) that decay into π+π-
// Counts how many π+ and π- have associated reconstructed objects

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

// Helper function to find reconstructed particle associated with a true particle
AnaParticlePD* FindReconstructedParticle(AnaSpill* spill, Int_t trueParticleID) {
    for (UInt_t bunchIdx = 0; bunchIdx < spill->Bunches.size(); bunchIdx++) {
        AnaBunchB* bunch = static_cast<AnaBunchB*>(spill->Bunches[bunchIdx]);
        if (!bunch) continue;

        for (UInt_t j = 0; j < bunch->Particles.size(); j++) {
            AnaParticleB* particle = bunch->Particles[j];
            if (!particle) continue;

            // Cast to AnaParticlePD to access ProtoDUNE-specific data
            AnaParticlePD* reconPart = dynamic_cast<AnaParticlePD*>(particle);
            if (!reconPart) continue;

            // Check if this reconstructed particle is associated with the true particle
            AnaTrueParticleB* associatedTrue = reconPart->GetTrueParticle();
            if (associatedTrue && associatedTrue->ID == trueParticleID) {
                return reconPart;
            }
        }
    }
    return nullptr;
}

void K0ReconstructionAnalysis(const char* dataDir = "/Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA") {

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

    // Global statistics counters
    int totalEvents = 0;

    // K0 statistics by parent type
    int totalK0FromKPlusBeam = 0;
    int totalK0FromKPlusNonBeam = 0;
    int eventsWithK0FromKPlusBeam = 0;
    int eventsWithK0FromKPlusNonBeam = 0;

    // K0 -> π+π- decay statistics
    int k0FromKPlusBeam_PiPiDecay = 0;
    int k0FromKPlusNonBeam_PiPiDecay = 0;
    int eventsWithK0FromKPlusBeam_PiPiDecay = 0;
    int eventsWithK0FromKPlusNonBeam_PiPiDecay = 0;

    // Reconstruction statistics for π+π- daughters
    int piPlusFromKPlusBeamK0_Reconstructed = 0;
    int piMinusFromKPlusBeamK0_Reconstructed = 0;
    int piPlusFromKPlusNonBeamK0_Reconstructed = 0;
    int piMinusFromKPlusNonBeamK0_Reconstructed = 0;

    // Total π+π- daughters (for efficiency calculation)
    int totalPiPlusFromKPlusBeamK0 = 0;
    int totalPiMinusFromKPlusBeamK0 = 0;
    int totalPiPlusFromKPlusNonBeamK0 = 0;
    int totalPiMinusFromKPlusNonBeamK0 = 0;

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

        // File-specific counters
        int fileK0FromKPlusBeam = 0;
        int fileK0FromKPlusNonBeam = 0;
        int fileEventsWithK0FromKPlusBeam = 0;
        int fileEventsWithK0FromKPlusNonBeam = 0;
        int fileK0FromKPlusBeam_PiPiDecay = 0;
        int fileK0FromKPlusNonBeam_PiPiDecay = 0;
        int fileEventsWithK0FromKPlusBeam_PiPiDecay = 0;
        int fileEventsWithK0FromKPlusNonBeam_PiPiDecay = 0;
        int filePiPlusFromKPlusBeamK0_Reconstructed = 0;
        int filePiMinusFromKPlusBeamK0_Reconstructed = 0;
        int filePiPlusFromKPlusNonBeamK0_Reconstructed = 0;
        int filePiMinusFromKPlusNonBeamK0_Reconstructed = 0;
        int fileTotalPiPlusFromKPlusBeamK0 = 0;
        int fileTotalPiMinusFromKPlusBeamK0 = 0;
        int fileTotalPiPlusFromKPlusNonBeamK0 = 0;
        int fileTotalPiMinusFromKPlusNonBeamK0 = 0;

        // Process each event
        for (Long64_t i = 0; i < nEntries; i++) {
            globalEventCounter++;
            totalEvents++;

            MiniTree->GetEntry(i);

            if (!spill) continue;

            // Get beam particle information
            AnaParticleB* beamParticle = nullptr;
            if (spill->Beam) {
                AnaBeam* beam = static_cast<AnaBeam*>(spill->Beam);
                if (beam && beam->BeamParticle) {
                    beamParticle = beam->BeamParticle;
                }
            }

            AnaTrueParticleB* beamTrue = nullptr;
            if (beamParticle) {
                beamTrue = beamParticle->GetTrueParticle();
            }

            UInt_t nTrueParticles = spill->TrueParticles.size();

            // Find all K+ particles
            std::vector<AnaTrueParticleB*> kPlusParticles;
            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;
                if (truePart->PDG == 321) { // K+
                    kPlusParticles.push_back(truePart);
                }
            }

            // Find all K0 particles and categorize by parent type
            std::vector<AnaTrueParticleB*> k0FromKPlusBeam;
            std::vector<AnaTrueParticleB*> k0FromKPlusNonBeam;

            for (UInt_t j = 0; j < nTrueParticles; j++) {
                AnaTrueParticleB* truePart = spill->TrueParticles[j];
                if (!truePart) continue;

                // Check if this is a K0 (PDG = 310)
                if (truePart->PDG == 310) {
                    // Check if this K0 is a daughter of any K+ particle
                    for (UInt_t k = 0; k < kPlusParticles.size(); k++) {
                        AnaTrueParticleB* kPlus = kPlusParticles[k];
                        if (truePart->ParentID == kPlus->ID) {
                            // Check if this K+ is the beam particle
                            if (beamTrue && kPlus->ID == beamTrue->ID) {
                                k0FromKPlusBeam.push_back(truePart);
                                totalK0FromKPlusBeam++;
                                fileK0FromKPlusBeam++;
                            } else {
                                k0FromKPlusNonBeam.push_back(truePart);
                                totalK0FromKPlusNonBeam++;
                                fileK0FromKPlusNonBeam++;
                            }
                            break; // Found the parent, no need to check other K+ particles
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

            // Analyze K0 decays for beam K+ daughters
            for (UInt_t k0Idx = 0; k0Idx < k0FromKPlusBeam.size(); k0Idx++) {
                AnaTrueParticleB* k0Particle = k0FromKPlusBeam[k0Idx];

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

                // Check if K0 decays into exactly π+π-
                if (k0Daughters.size() == 2) {
                    int piPlusCount = 0;
                    int piMinusCount = 0;
                    AnaTrueParticleB* piPlusDaughter = nullptr;
                    AnaTrueParticleB* piMinusDaughter = nullptr;

                    for (UInt_t dIdx = 0; dIdx < k0Daughters.size(); dIdx++) {
                        AnaTrueParticleB* daughter = k0Daughters[dIdx];
                        int pdg = daughter->PDG;

                        if (pdg == 211) {
                            piPlusCount++;
                            piPlusDaughter = daughter;
                        } else if (pdg == -211) {
                            piMinusCount++;
                            piMinusDaughter = daughter;
                        }
                    }

                    // Check if this is exactly π+π- decay
                    if (piPlusCount == 1 && piMinusCount == 1) {
                        k0FromKPlusBeam_PiPiDecay++;
                        fileK0FromKPlusBeam_PiPiDecay++;

                        // Check reconstruction for π+
                        if (piPlusDaughter) {
                            totalPiPlusFromKPlusBeamK0++;
                            fileTotalPiPlusFromKPlusBeamK0++;

                            // Look for reconstructed particle associated with this π+
                            AnaParticlePD* reconPart = FindReconstructedParticle(spill, piPlusDaughter->ID);
                            if (reconPart) {
                                piPlusFromKPlusBeamK0_Reconstructed++;
                                filePiPlusFromKPlusBeamK0_Reconstructed++;
                            }
                        }

                        // Check reconstruction for π-
                        if (piMinusDaughter) {
                            totalPiMinusFromKPlusBeamK0++;
                            fileTotalPiMinusFromKPlusBeamK0++;

                            // Look for reconstructed particle associated with this π-
                            AnaParticlePD* reconPart = FindReconstructedParticle(spill, piMinusDaughter->ID);
                            if (reconPart) {
                                piMinusFromKPlusBeamK0_Reconstructed++;
                                filePiMinusFromKPlusBeamK0_Reconstructed++;
                            }
                        }
                    }
                }
            }

            // Analyze K0 decays for non-beam K+ daughters
            for (UInt_t k0Idx = 0; k0Idx < k0FromKPlusNonBeam.size(); k0Idx++) {
                AnaTrueParticleB* k0Particle = k0FromKPlusNonBeam[k0Idx];

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

                // Check if K0 decays into exactly π+π-
                if (k0Daughters.size() == 2) {
                    int piPlusCount = 0;
                    int piMinusCount = 0;
                    AnaTrueParticleB* piPlusDaughter = nullptr;
                    AnaTrueParticleB* piMinusDaughter = nullptr;

                    for (UInt_t dIdx = 0; dIdx < k0Daughters.size(); dIdx++) {
                        AnaTrueParticleB* daughter = k0Daughters[dIdx];
                        int pdg = daughter->PDG;

                        if (pdg == 211) {
                            piPlusCount++;
                            piPlusDaughter = daughter;
                        } else if (pdg == -211) {
                            piMinusCount++;
                            piMinusDaughter = daughter;
                        }
                    }

                    // Check if this is exactly π+π- decay
                    if (piPlusCount == 1 && piMinusCount == 1) {
                        k0FromKPlusNonBeam_PiPiDecay++;
                        fileK0FromKPlusNonBeam_PiPiDecay++;

                        // Check reconstruction for π+
                        if (piPlusDaughter) {
                            totalPiPlusFromKPlusNonBeamK0++;
                            fileTotalPiPlusFromKPlusNonBeamK0++;

                            // Look for reconstructed particle associated with this π+
                            AnaParticlePD* reconPart = FindReconstructedParticle(spill, piPlusDaughter->ID);
                            if (reconPart) {
                                piPlusFromKPlusNonBeamK0_Reconstructed++;
                                filePiPlusFromKPlusNonBeamK0_Reconstructed++;
                            }
                        }

                        // Check reconstruction for π-
                        if (piMinusDaughter) {
                            totalPiMinusFromKPlusNonBeamK0++;
                            fileTotalPiMinusFromKPlusNonBeamK0++;

                            // Look for reconstructed particle associated with this π-
                            AnaParticlePD* reconPart = FindReconstructedParticle(spill, piMinusDaughter->ID);
                            if (reconPart) {
                                piMinusFromKPlusNonBeamK0_Reconstructed++;
                                filePiMinusFromKPlusNonBeamK0_Reconstructed++;
                            }
                        }
                    }
                }
            }

            // Update event counters for π+π- decays
            if (fileK0FromKPlusBeam_PiPiDecay > 0) {
                eventsWithK0FromKPlusBeam_PiPiDecay++;
                fileEventsWithK0FromKPlusBeam_PiPiDecay++;
            }
            if (fileK0FromKPlusNonBeam_PiPiDecay > 0) {
                eventsWithK0FromKPlusNonBeam_PiPiDecay++;
                fileEventsWithK0FromKPlusNonBeam_PiPiDecay++;
            }
        }

        // Print file summary
        std::cout << "  File summary:" << std::endl;
        std::cout << "    K0 from K+ beam: " << fileK0FromKPlusBeam << std::endl;
        std::cout << "    K0 from K+ non-beam: " << fileK0FromKPlusNonBeam << std::endl;
        std::cout << "    K0 from K+ beam -> π+π-: " << fileK0FromKPlusBeam_PiPiDecay << std::endl;
        std::cout << "    K0 from K+ non-beam -> π+π-: " << fileK0FromKPlusNonBeam_PiPiDecay << std::endl;
        std::cout << "    π+ from K+ beam K0 reconstructed: " << filePiPlusFromKPlusBeamK0_Reconstructed
                  << "/" << fileTotalPiPlusFromKPlusBeamK0 << std::endl;
        std::cout << "    π- from K+ beam K0 reconstructed: " << filePiMinusFromKPlusBeamK0_Reconstructed
                  << "/" << fileTotalPiMinusFromKPlusBeamK0 << std::endl;
        std::cout << "    π+ from K+ non-beam K0 reconstructed: " << filePiPlusFromKPlusNonBeamK0_Reconstructed
                  << "/" << fileTotalPiPlusFromKPlusNonBeamK0 << std::endl;
        std::cout << "    π- from K+ non-beam K0 reconstructed: " << filePiMinusFromKPlusNonBeamK0_Reconstructed
                  << "/" << fileTotalPiMinusFromKPlusNonBeamK0 << std::endl;
        std::cout << std::endl;

        f->Close();
        delete f;
    }

    // Print final results
    std::cout << "==========================================" << std::endl;
    std::cout << "K0 RECONSTRUCTION ANALYSIS SUMMARY" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << std::endl;

    std::cout << "K0 FROM K+ BEAM PARTICLES:" << std::endl;
    std::cout << "  Total K0: " << totalK0FromKPlusBeam << std::endl;
    std::cout << "  Events with K0: " << eventsWithK0FromKPlusBeam << std::endl;
    std::cout << "  K0 -> π+π- decays: " << k0FromKPlusBeam_PiPiDecay << std::endl;
    std::cout << "  Events with K0 -> π+π-: " << eventsWithK0FromKPlusBeam_PiPiDecay << std::endl;
    std::cout << "  π+ reconstructed: " << piPlusFromKPlusBeamK0_Reconstructed
              << "/" << totalPiPlusFromKPlusBeamK0
              << " (" << (totalPiPlusFromKPlusBeamK0 > 0 ? 100.0 * piPlusFromKPlusBeamK0_Reconstructed / totalPiPlusFromKPlusBeamK0 : 0.0) << "%)" << std::endl;
    std::cout << "  π- reconstructed: " << piMinusFromKPlusBeamK0_Reconstructed
              << "/" << totalPiMinusFromKPlusBeamK0
              << " (" << (totalPiMinusFromKPlusBeamK0 > 0 ? 100.0 * piMinusFromKPlusBeamK0_Reconstructed / totalPiMinusFromKPlusBeamK0 : 0.0) << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << "K0 FROM K+ NON-BEAM PARTICLES:" << std::endl;
    std::cout << "  Total K0: " << totalK0FromKPlusNonBeam << std::endl;
    std::cout << "  Events with K0: " << eventsWithK0FromKPlusNonBeam << std::endl;
    std::cout << "  K0 -> π+π- decays: " << k0FromKPlusNonBeam_PiPiDecay << std::endl;
    std::cout << "  Events with K0 -> π+π-: " << eventsWithK0FromKPlusNonBeam_PiPiDecay << std::endl;
    std::cout << "  π+ reconstructed: " << piPlusFromKPlusNonBeamK0_Reconstructed
              << "/" << totalPiPlusFromKPlusNonBeamK0
              << " (" << (totalPiPlusFromKPlusNonBeamK0 > 0 ? 100.0 * piPlusFromKPlusNonBeamK0_Reconstructed / totalPiPlusFromKPlusNonBeamK0 : 0.0) << "%)" << std::endl;
    std::cout << "  π- reconstructed: " << piMinusFromKPlusNonBeamK0_Reconstructed
              << "/" << totalPiMinusFromKPlusNonBeamK0
              << " (" << (totalPiMinusFromKPlusNonBeamK0 > 0 ? 100.0 * piMinusFromKPlusNonBeamK0_Reconstructed / totalPiMinusFromKPlusNonBeamK0 : 0.0) << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << "COMBINED K0 FROM ALL K+:" << std::endl;
    std::cout << "  Total K0: " << (totalK0FromKPlusBeam + totalK0FromKPlusNonBeam) << std::endl;
    std::cout << "  Events with K0: " << (eventsWithK0FromKPlusBeam + eventsWithK0FromKPlusNonBeam) << std::endl;
    std::cout << "  K0 -> π+π- decays: " << (k0FromKPlusBeam_PiPiDecay + k0FromKPlusNonBeam_PiPiDecay) << std::endl;
    std::cout << "  Events with K0 -> π+π-: " << (eventsWithK0FromKPlusBeam_PiPiDecay + eventsWithK0FromKPlusNonBeam_PiPiDecay) << std::endl;
    std::cout << "  π+ reconstructed: " << (piPlusFromKPlusBeamK0_Reconstructed + piPlusFromKPlusNonBeamK0_Reconstructed)
              << "/" << (totalPiPlusFromKPlusBeamK0 + totalPiPlusFromKPlusNonBeamK0)
              << " (" << ((totalPiPlusFromKPlusBeamK0 + totalPiPlusFromKPlusNonBeamK0) > 0 ?
                  100.0 * (piPlusFromKPlusBeamK0_Reconstructed + piPlusFromKPlusNonBeamK0_Reconstructed) /
                  (totalPiPlusFromKPlusBeamK0 + totalPiPlusFromKPlusNonBeamK0) : 0.0) << "%)" << std::endl;
    std::cout << "  π- reconstructed: " << (piMinusFromKPlusBeamK0_Reconstructed + piMinusFromKPlusNonBeamK0_Reconstructed)
              << "/" << (totalPiMinusFromKPlusBeamK0 + totalPiMinusFromKPlusNonBeamK0)
              << " (" << ((totalPiMinusFromKPlusBeamK0 + totalPiMinusFromKPlusNonBeamK0) > 0 ?
                  100.0 * (piMinusFromKPlusBeamK0_Reconstructed + piMinusFromKPlusNonBeamK0_Reconstructed) /
                  (totalPiMinusFromKPlusBeamK0 + totalPiMinusFromKPlusNonBeamK0) : 0.0) << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << "Analysis completed successfully!" << std::endl;
}
