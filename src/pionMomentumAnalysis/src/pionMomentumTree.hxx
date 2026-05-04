#ifndef pionMomentumTree_h
#define pionMomentumTree_h

#include "OutputManager.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include <vector>

namespace pionMomentumTree {

  // Registers micro-tree branches for main-track MCS/TLEFit and beam true-info (see enumPionMomentumMicroTrees).
  void AddPionMomentumVariables_BeamParticleReco(OutputManager& output);
  // Fills those branches from reco particles, MCS vectors, and optional beam true particle.
  void FillPionMomentumVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beampart,
                                                  const std::vector<double>& mainmcsdeltathetaPair,
                                                  const std::vector<double>& mainmcsdeltathetaxzPair,
                                                  const std::vector<double>& mainmcsdeltathetayzPair,
                                                  const std::vector<double>& mainmcssegmentlengthSingle,
                                                  const std::vector<double>& mainmcssegmentfitchi2ndfSingle,
                                                  const std::vector<double>& mainmcssegmentrrtoendSingle,
                                                  double mainmcsmomentumGeV,
                                                  double maintlefitmomentumGeV,
                                                  double mainbraggchi2pibb,
                                                  Int_t mainbraggdedxnhits,
                                                  Int_t mainntruetrajpoints,
                                                  const std::vector<double>& maintruemcsdeltathetaPair,
                                                  const std::vector<double>& maintruemcsdeltathetaxzPair,
                                                  const std::vector<double>& maintruemcsdeltathetayzPair,
                                                  const std::vector<double>& maintruemcssegmentlengthSingle,
                                                  const std::vector<double>& maintruemcssegmentfitchi2ndfSingle,
                                                  const std::vector<double>& maintruemcssegmentrrtoendSingle,
                                                  double maintruemcsmomentumGeV);

  // Indices for pion-momentum micro-tree branches (see AddPionMomentumVariables_BeamParticleReco).
  enum enumPionMomentumMicroTrees {
    mainispandora = standardPDTree::enumStandardMicroTreesLast_standardPDTree + 1,  // main track: Pandora reconstruction flag
    maintruepdg,                                                                  // main track: true PDG
    maintruestartmomentum,                                                        // main track: true |p| at start [GeV/c]
    maintrueendmomentum,                                                          // main track: true |p| at end [GeV/c]
    maintruelength,                                                               // main track: true path length [cm]
    maintrueid,                                                                   // main track: true particle ID
    beamtrueid,                                                                   // beam particle: true particle ID
    mainlength,                                                                   // main track: reconstructed length [cm]
    mainmcsnsegments,                                                             // main track: MCS segment count (vector size)
    mainmcsnpairs,                                                                // main track: MCS deflection-pair count (vector size)
    mainmcsdeltatheta,                                                            // main track: MCS sqrt(theta_xz^2+theta_yz^2) [rad], per pair
    mainmcsdeltathetaxz,                                                          // main track: MCS projected angle theta_xz [rad], per pair
    mainmcsdeltathetayz,                                                          // main track: MCS projected angle theta_yz [rad], per pair
    mainmcssegmentlength,                                                         // main track: MCS mean segment length [cm], per segment
    mainmcssegmentfitchi2ndf,                                                   // main track: MCS segment straight-line fit quality [cm^2]
    mainmcssegmentrrtoend,                                                      // main track: MCS segment midpoint RR distance to track end [cm]
    mainmcsmomentum,                                                              // main track: MCS momentum estimate [GeV/c]
    maintlefitmomentum,                                                           // main track: TLEFit momentum estimate [GeV/c]
    mainnhitscollection,                                                          // main track: number of collection-plane hits
    mainntrjpoints,                                                               // main track: number of trajectory points
    mainnvalidrrhits,                                                             // main track: collection hits with finite residual range
    mainnvalidxyzhits,                                                            // main track: collection hits with valid XYZ
    mainnvalidxyztrj,                                                             // main track: trajectory points with valid XYZ
    mainntruetrajpoints,  // true traj point count used for true MCS (beam when IDs match and main has no points)
    maintruemcsnsegments,                                                         // main track: true-traj MCS segment count (vector size)
    maintruemcsnpairs,                                                            // main track: true-traj MCS pair count (vector size)
    maintruemcsdeltatheta,                                                        // main track: true-traj MCS combined angle [rad], per pair
    maintruemcsdeltathetaxz,                                                      // main track: true-traj MCS theta_xz [rad], per pair
    maintruemcsdeltathetayz,                                                      // main track: true-traj MCS theta_yz [rad], per pair
    maintruemcssegmentlength,                                                     // main track: true-traj MCS mean segment length [cm]
    maintruemcssegmentfitchi2ndf,                                               // main track: true-traj MCS segment fit quality [cm^2]
    maintruemcssegmentrrtoend,                                                    // main track: true-traj MCS segment RR to end [cm]
    maintruemcsmomentum,                                                          // main track: true-traj MCS momentum [GeV/c]
    beamtruelaststeplen,                                                          // beam true: last step length [cm]
    beamtruenelastic,                                                             // beam true: number of elastic scatters
    beamtrueidetotaldep,                                                          // beam true: total IDE energy deposit
    beamtruenhits,                                                                // beam true: hit count
    beamtrueelasticn,                                                             // beam true: elastic-scatter list size (vector size)
    beamtrueelasticcostheta,                                                      // beam true: elastic scatter cos(theta)
    beamtrueelasticx,                                                             // beam true: elastic scatter X [cm]
    beamtrueelasticy,                                                             // beam true: elastic scatter Y [cm]
    beamtrueelasticz,                                                             // beam true: elastic scatter Z [cm]
    beamtrueelasticdeltae,                                                      // beam true: elastic scatter delta E
    beamtrueelasticideedep,                                                       // beam true: elastic scatter step IDE edep
    beamtrueprocessesn,                                                           // beam true: process-name list size (vector size)
    beamtrueprocesses,                                                            // beam true: process names along track (strings)
    beamtrueincidentn,                                                            // beam true: thin-slice incident-energy list size
    beamtrueincidentenergies,                                                     // beam true: thin-slice incident energies
    mainbraggchi2pibb,  // main track: χ²_π± mean per arXiv:2409.18288 Eq. 6.1 (BB ⟨dE/dx⟩ vs meas, RR<param, σ param); smaller≈stopping
    mainbraggdedxnhits,  // main track: collection hits used for Bragg-window χ² (Eq. 6.1)
  };

}  // namespace pionMomentumTree

#endif
