#ifndef pionMomentumTree_h
#define pionMomentumTree_h

#include "OutputManager.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include <vector>

namespace pionMomentumTree {

  // Registers micro-tree branches for main-track MCS/TLEFit and beam true-info (see enumPionMomentumMicroTrees).
  void AddPionMomentumVariables_BeamParticleReco(OutputManager& output);
  /// Per true_beam daughter (same ID order as beam true \c Daughters): true PDG, true |p| at start, Bragg-window π template χ²/hit on best-matched reco track (if any).
  void AddPionMomentumVariables_BeamTrueDaughterBragg(OutputManager& output);
  /// Flat rows: beam-daughter TLE truncation scan (optional; see RunBeamDaughterTLETruncationScan).
  void AddPionMomentumVariables_BeamDaughterTleTruncScan(OutputManager& output);
  /// Flat rows: beam-daughter MCS truncation scan (optional; see RunBeamDaughterMCSTruncationScan).
  void AddPionMomentumVariables_BeamDaughterMcsTruncScan(OutputManager& output);
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

  void FillPionMomentumVariables_BeamTrueDaughterBragg(OutputManager& output, const std::vector<Int_t>& trueIds,
                                                        const std::vector<Int_t>& truePdgs,
                                                        const std::vector<Float_t>& trueStartMomGeV,
                                                        const std::vector<Float_t>& trueEndMomGeV,
                                                        const std::vector<Float_t>& trueStartEkinGeV,
                                                        const std::vector<Float_t>& trueEndEkinGeV,
                                                        const std::vector<Float_t>& csdaMomentumGeV,
                                                        const std::vector<Float_t>& csdaKineticEnergyGeV,
                                                        const std::vector<Int_t>& trueEndProcess,
                                                        const std::vector<Float_t>& trueEndToRecoStartDistCm,
                                                        const std::vector<Float_t>& trueEndToRecoEndDistCm,
                                                        const std::vector<double>& braggChi2MeanPerHit,
                                                        const std::vector<Int_t>& braggNhits);

  void FillPionMomentumVariables_BeamDaughterTleTruncScan(
      OutputManager& output, const std::vector<Int_t>& truncDauIdx, const std::vector<Int_t>& truncK,
      const std::vector<Int_t>& truncNhitsInt, const std::vector<Float_t>& truncTrueEkin0GeV,
      const std::vector<Float_t>& truncPtleGeV, const std::vector<Float_t>& truncFracRes);
  void FillPionMomentumVariables_BeamDaughterMcsTruncScan(
      OutputManager& output, const std::vector<Int_t>& truncDauIdx, const std::vector<Int_t>& truncK,
      const std::vector<Int_t>& truncNsegments, const std::vector<Float_t>& truncTrueEkin0GeV,
      const std::vector<Float_t>& truncPmcsGeV, const std::vector<Float_t>& truncFracRes);

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
    maintrueendtorecoendcm,                                                       // main track: distance true-end to reco-end [cm]
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
    mainbraggchi2pibb,  // main track: mean χ²/hit like Chi2PID(211) (dedx_range_pi template) in 0<RR<Rmax window
    mainbraggdedxnhits,  // main track: collection hits used for Bragg-window template χ²
    beamdaughtern,           // beam true first-level daughters: vector size (same order as true_beam daughter IDs)
    beamdaughtertrueid,     // per daughter: true particle ID
    beamdaughtertruepdg,    // per daughter: true PDG
    beamdaughtertruestartmom,  // per daughter: true |p| at start [GeV/c]
    beamdaughtertrueendmom,    // per daughter: true |p| at end [GeV/c]
    beamdaughtertruestartekin,  // per daughter: true kinetic energy at start [GeV]
    beamdaughtertrueendekin,    // per daughter: true kinetic energy at end [GeV]
    beamdaughtercsdamomentum,   // per daughter: CSDA momentum estimate from reco range [GeV/c]
    beamdaughtercsdakineticenergy,  // per daughter: CSDA kinetic-energy estimate from reco range [GeV]
    beamdaughtertrueendprocess,  // per daughter: true end process enum (AnaTrueParticleB::ProcessEnum as int)
    beamdaughtertrueendtorecostartcm,  // per daughter: distance true-end to reco-start [cm]
    beamdaughtertrueendtorecoendcm,    // per daughter: distance true-end to reco-end [cm]
    beamdaughterbraggchi2pibb,  // per daughter: Bragg-window χ²_π± Eq. 6.1 on reco daughter (-999 if invalid/fail)
    beamdaughterbraggdedxnhits,  // per daughter: hits used for χ² (-999 if unmatched / fail)
    beamdaughtertletruncn,       // TLE trunc scan: row count (vector size for following branches)
    beamdaughtertletruncdauidx,  // TLE trunc: daughter index in beamdaughter* arrays [0, beamdaughtern)
    beamdaughtertletrunck,       // TLE trunc: extra end skip k (effective skipHitsLast = base + k)
    beamdaughtertletruncnhitsint,  // TLE trunc: interior collection hits after skips + dE/dx window
    beamdaughtertletrunctrueekin0,  // TLE trunc: true start Ekin [GeV] (reference)
    beamdaughtertletruncptle,       // TLE trunc: TLEFit momentum [GeV/c]
    beamdaughtertletruncfracres,    // TLE trunc: (E_K^TLE - E_K^true) / E_K^true, E from momentum & mass
    beamdaughtermcstruncn,       // MCS trunc scan: row count (vector size for following branches)
    beamdaughtermcstruncdauidx,  // MCS trunc: daughter index in beamdaughter* arrays [0, beamdaughtern)
    beamdaughtermcstrunck,       // MCS trunc: extra low-RR segment drop k (dropLast = 3 + k)
    beamdaughtermcstruncnsegments,  // MCS trunc: surviving scattering segments after drops
    beamdaughtermcstrunctrueekin0,  // MCS trunc: true start Ekin [GeV] (reference)
    beamdaughtermcstruncpmcs,       // MCS trunc: MCS momentum [GeV/c]
    beamdaughtermcstruncfracres,    // MCS trunc: (E_K^MCS - E_K^true) / E_K^true
  };

}  // namespace pionMomentumTree

#endif
