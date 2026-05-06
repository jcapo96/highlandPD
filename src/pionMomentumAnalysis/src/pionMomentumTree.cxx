#include "pionMomentumTree.hxx"
#include <algorithm>

namespace pionMomentumTree {

//********************************************************************
void AddPionMomentumVariables_BeamParticleReco(OutputManager& output) {
//********************************************************************
  const UInt_t nmaxmcs = 4096;
  AddVarI(output, mainispandora, "Main track Pandora flag");
  AddVarI(output, maintruepdg, "Main track true PDG");
  AddVarF(output, maintruestartmomentum, "Main track true start momentum [GeV/c]");
  AddVarF(output, maintrueendmomentum, "Main track true end momentum [GeV/c]");
  AddVarF(output, maintruelength, "Main track true length [cm]");
  AddVarI(output, maintrueid, "Main track true ID");
  AddVarI(output, beamtrueid, "Beam particle true ID");
  AddVarF(output, mainlength, "Main track reconstructed length [cm]");
  AddVarF(output, maintrueendtorecoendcm, "Main track distance true-end to reco-end [cm] (-999 if invalid)");
  AddVarMaxSizeVF(output, mainmcsdeltatheta, "Main track MCS projected-angle magnitude sqrt(theta_xz^2+theta_yz^2) [rad]", mainmcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcsdeltathetaxz, "Main track MCS projected angle theta_xz [rad]", mainmcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcsdeltathetayz, "Main track MCS projected angle theta_yz [rad]", mainmcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcssegmentlength, "Main track MCS mean segment length [cm]", mainmcsnsegments, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcssegmentfitchi2ndf, "Main track MCS segment linear-fit quality (chi2/ndf-like) [cm^{2}]", mainmcsnsegments, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcssegmentrrtoend, "Main track MCS segment distance to track end (RR midpoint) [cm]", mainmcsnsegments, nmaxmcs);
  AddVarF(output, mainmcsmomentum, "Main track MCS momentum estimate [GeV/c]");
  AddVarF(output, maintlefitmomentum, "Main track TLEFit momentum estimate [GeV/c]");
  AddVarI(output, mainnhitscollection, "Main track number of collection-plane hits");
  AddVarI(output, mainntrjpoints, "Main track number of trajectory points");
  AddVarI(output, mainnvalidrrhits, "Main track number of collection hits with valid residual range");
  AddVarI(output, mainnvalidxyzhits, "Main track collection hits with valid XYZ");
  AddVarI(output, mainnvalidxyztrj, "Main track trajectory points with valid XYZ");
  AddVarI(output, mainntruetrajpoints,
          "True particle TrjPoints count used for true MCS (beam TrueObject when same ID as main and main has no points)");
  AddVarMaxSizeVF(output, maintruemcsdeltatheta,
                   "Main track true-traj MCS projected-angle magnitude sqrt(theta_xz^2+theta_yz^2) [rad]",
                   maintruemcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, maintruemcsdeltathetaxz, "Main track true-traj MCS projected angle theta_xz [rad]",
                  maintruemcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, maintruemcsdeltathetayz, "Main track true-traj MCS projected angle theta_yz [rad]",
                  maintruemcsnpairs, nmaxmcs);
  AddVarMaxSizeVF(output, maintruemcssegmentlength, "Main track true-traj MCS mean segment length [cm]",
                  maintruemcsnsegments, nmaxmcs);
  AddVarMaxSizeVF(output, maintruemcssegmentfitchi2ndf,
                  "Main track true-traj MCS segment linear-fit quality (chi2/ndf-like) [cm^{2}]", maintruemcsnsegments,
                  nmaxmcs);
  AddVarMaxSizeVF(output, maintruemcssegmentrrtoend,
                  "Main track true-traj MCS segment distance to track end (RR midpoint) [cm]", maintruemcsnsegments,
                  nmaxmcs);
  AddVarF(output, maintruemcsmomentum, "Main track true-traj MCS momentum estimate [GeV/c]");

  const UInt_t nmaxelastic = 512;
  const UInt_t nmaxproc = 512;
  const UInt_t nmaxincident = 4096;
  AddVarD(output, beamtruelaststeplen, "Beam true particle last step length (true_beam_last_len) [cm]");
  AddVarI(output, beamtruenelastic, "Beam true particle number of elastic scatters");
  AddVarD(output, beamtrueidetotaldep, "Beam true particle total IDE energy deposit (true_beam_IDE_totalDep)");
  AddVarI(output, beamtruenhits, "Beam true particle number of hits (true_beam_nHits)");
  AddVarMaxSizeVD(output, beamtrueelasticcostheta, "Beam elastic scatter cos(theta)", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVD(output, beamtrueelasticx, "Beam elastic scatter X [cm]", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVD(output, beamtrueelasticy, "Beam elastic scatter Y [cm]", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVD(output, beamtrueelasticz, "Beam elastic scatter Z [cm]", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVD(output, beamtrueelasticdeltae, "Beam elastic scatter delta E", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVD(output, beamtrueelasticideedep, "Beam elastic scatter step IDE edep", beamtrueelasticn, nmaxelastic);
  AddVarMaxSizeVC(output, beamtrueprocesses, "Beam true processes along track (strings)", beamtrueprocessesn, nmaxproc);
  AddVarMaxSizeVD(output, beamtrueincidentenergies, "Beam thin-slice incident energies", beamtrueincidentn, nmaxincident);
  AddVarF(output, mainbraggchi2pibb,
          "Main track mean chi2/hit like Chi2PID(211): dedx_range_pi template + errors, hits with 0<RR<Rmax only");
  AddVarI(output, mainbraggdedxnhits, "Main track collection hits used for Bragg-window Chi2PID-style chi2_pi");
}

//********************************************************************
void AddPionMomentumVariables_BeamTrueDaughterBragg(OutputManager& output) {
//********************************************************************
  constexpr UInt_t nmaxdau = 64;
  AddVarMaxSizeVI(output, beamdaughtertrueid, "Beam true daughter: true particle ID", beamdaughtern, nmaxdau);
  AddVarMaxSizeVI(output, beamdaughtertruepdg, "Beam true daughter: true PDG", beamdaughtern, nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertruestartmom, "Beam true daughter: true |p| at start [GeV/c]", beamdaughtern,
                  nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertrueendmom, "Beam true daughter: true |p| at end [GeV/c]", beamdaughtern,
                  nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertruestartekin, "Beam true daughter: true kinetic energy at start [GeV]",
                  beamdaughtern, nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertrueendekin, "Beam true daughter: true kinetic energy at end [GeV]",
                  beamdaughtern, nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtercsdamomentum,
                  "Beam daughter reco track: pion CSDA momentum estimate from range [GeV/c]", beamdaughtern, nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtercsdakineticenergy,
                  "Beam daughter reco track: pion CSDA kinetic energy estimate from range [GeV]", beamdaughtern, nmaxdau);
  AddVarMaxSizeVI(output, beamdaughtertrueendprocess,
                  "Beam true daughter: true end process enum (AnaTrueParticleB::ProcessEnum as int)", beamdaughtern,
                  nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertrueendtorecostartcm,
                  "Beam true daughter: distance true-end to reco-start [cm] (-999 if invalid)", beamdaughtern,
                  nmaxdau);
  AddVarMaxSizeVF(output, beamdaughtertrueendtorecoendcm,
                  "Beam true daughter: distance true-end to reco-end [cm] (-999 if invalid)", beamdaughtern, nmaxdau);
  AddVarMaxSizeVF(output, beamdaughterbraggchi2pibb,
                  "Beam true daughter: Bragg-window chi2_pi Eq. 6.1 on reco daughter track (-999 if invalid/fail)",
                  beamdaughtern, nmaxdau);
  AddVarMaxSizeVI(output, beamdaughterbraggdedxnhits,
                  "Beam true daughter: collection hits used in Bragg chi2 (-999 if none/fail)", beamdaughtern, nmaxdau);
}

//********************************************************************
void AddPionMomentumVariables_BeamDaughterTleTruncScan(OutputManager& output) {
//********************************************************************
  constexpr UInt_t nmaxrows = 4096;
  AddVarMaxSizeVI(output, beamdaughtertletruncdauidx,
                  "Beam daughter TLE trunc scan: daughter index in beamdaughter* arrays [0, beamdaughtern)",
                  beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVI(output, beamdaughtertletrunck,
                  "Beam daughter TLE trunc scan: extra end-veto k (skipHitsLast = FreeRange base + k)",
                  beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVI(output, beamdaughtertletruncnhitsint,
                  "Beam daughter TLE trunc scan: interior hits after skips + dE/dx window", beamdaughtertletruncn,
                  nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtertletrunctrueekin0,
                  "Beam daughter TLE trunc scan: true start Ekin [GeV] (reference)", beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtertletruncekincsdafull,
                  "Beam daughter TLE trunc scan: CSDA Ekin from full reco track length [GeV] (reference for data)",
                  beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtertletruncptle, "Beam daughter TLE trunc scan: TLEFit momentum [GeV/c]",
                  beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtertletruncfracres,
                  "Beam daughter TLE trunc scan: (E_K^TLE - E_K^true) / E_K^true", beamdaughtertletruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtertletruncfracrescsda,
                  "Beam daughter TLE trunc scan: (E_K^TLE - E_K^CSDA,full) / E_K^CSDA,full", beamdaughtertletruncn,
                  nmaxrows);
}

//********************************************************************
void AddPionMomentumVariables_BeamDaughterMcsTruncScan(OutputManager& output) {
//********************************************************************
  constexpr UInt_t nmaxrows = 4096;
  AddVarMaxSizeVI(output, beamdaughtermcstruncdauidx,
                  "Beam daughter MCS trunc scan: daughter index in beamdaughter* arrays [0, beamdaughtern)",
                  beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVI(output, beamdaughtermcstrunck,
                  "Beam daughter MCS trunc scan: extra low-RR segment drop k (dropLast = 3 + k)",
                  beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVI(output, beamdaughtermcstruncnsegments,
                  "Beam daughter MCS trunc scan: surviving scattering segments after drops",
                  beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtermcstrunctrueekin0,
                  "Beam daughter MCS trunc scan: true start Ekin [GeV] (reference)", beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtermcstruncekincsdafull,
                  "Beam daughter MCS trunc scan: CSDA Ekin from full reco track length [GeV] (reference for data)",
                  beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtermcstruncpmcs, "Beam daughter MCS trunc scan: MCS momentum [GeV/c]",
                  beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtermcstruncfracres,
                  "Beam daughter MCS trunc scan: (E_K^MCS - E_K^true) / E_K^true", beamdaughtermcstruncn, nmaxrows);
  AddVarMaxSizeVF(output, beamdaughtermcstruncfracrescsda,
                  "Beam daughter MCS trunc scan: (E_K^MCS - E_K^CSDA,full) / E_K^CSDA,full", beamdaughtermcstruncn,
                  nmaxrows);
}

//********************************************************************
void FillPionMomentumVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beampart,
                                                const std::vector<double>& mainmcsdeltathetaPair,
                                                const std::vector<double>& mainmcsdeltathetaxzPair,
                                                const std::vector<double>& mainmcsdeltathetayzPair,
                                                const std::vector<double>& mainmcssegmentlengthSingle,
                                                const std::vector<double>& mainmcssegmentfitchi2ndfSingle,
                                                const std::vector<double>& mainmcssegmentrrtoendSingle,
                                                double mainmcsmomentumGeV,
                                                double maintlefitmomentumGeV,
                                                double mainbraggchi2pibbIn,
                                                Int_t mainbraggdedxnhitsIn,
                                                Int_t mainntruetrajpointsIn,
                                                const std::vector<double>& maintruemcsdeltathetaPair,
                                                const std::vector<double>& maintruemcsdeltathetaxzPair,
                                                const std::vector<double>& maintruemcsdeltathetayzPair,
                                                const std::vector<double>& maintruemcssegmentlengthSingle,
                                                const std::vector<double>& maintruemcssegmentfitchi2ndfSingle,
                                                const std::vector<double>& maintruemcssegmentrrtoendSingle,
                                                double maintruemcsmomentumGeV) {
//********************************************************************
  const size_t nmaxmcs = 4096;
  Int_t isPandora = -999;
  Int_t truePdg = -999;
  Float_t pStartTrue = -999.f;
  Float_t pEndTrue = -999.f;
  Float_t trueLength = -999.f;
  Float_t recoLength = -999.f;
  Float_t trueEndToRecoEndCm = -999.f;
  Int_t trueIdBeamParticle = -999;
  Int_t trueIdMainTrack = -999;
  Int_t nHitsCollection = 0;
  Int_t nTrjPoints = 0;
  Int_t nValidRRHits = 0;
  Int_t nValidXYZHits = 0;
  Int_t nValidXYZTrj = 0;
  if (part) isPandora = part->isPandora;
  if (part && part->TrueObject) {
    const AnaTrueParticlePD* tpart = static_cast<const AnaTrueParticlePD*>(part->TrueObject);
    if (tpart) {
      truePdg = tpart->PDG;
      pStartTrue = tpart->Momentum;
      pEndTrue = tpart->MomentumEnd;
      trueLength = tpart->Length;
      trueIdMainTrack = tpart->ID;
      const TVector3 trueEnd(static_cast<double>(tpart->PositionEnd[0]), static_cast<double>(tpart->PositionEnd[1]),
                             static_cast<double>(tpart->PositionEnd[2]));
      const TVector3 recoEnd(static_cast<double>(part->PositionEnd[0]), static_cast<double>(part->PositionEnd[1]),
                             static_cast<double>(part->PositionEnd[2]));
      const auto validXYZ = [](const TVector3& p) {
        return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) && p.X() > -900.0 &&
               p.Y() > -900.0 && p.Z() > -900.0;
      };
      if (validXYZ(trueEnd) && validXYZ(recoEnd)) trueEndToRecoEndCm = static_cast<Float_t>((trueEnd - recoEnd).Mag());
    }
  }
  if (part) {
    recoLength = part->Length;
    nHitsCollection = static_cast<Int_t>(part->Hits[2].size());
    nTrjPoints = static_cast<Int_t>(part->TrjPoints.size());
    for (size_t ihit = 0; ihit < part->Hits[2].size(); ++ihit) {
      const TVector3& p = part->Hits[2][ihit].Position;
      if (std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) &&
          p.X() > -900.0 && p.Y() > -900.0 && p.Z() > -900.0) {
        ++nValidXYZHits;
      }
      const double rr = static_cast<double>(part->Hits[2][ihit].ResidualRange);
      if (std::isfinite(rr) && rr > -900.0) ++nValidRRHits;
    }
    for (size_t itp = 0; itp < part->TrjPoints.size(); ++itp) {
      const TVector3& p = part->TrjPoints[itp].Position;
      if (std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) &&
          p.X() > -900.0 && p.Y() > -900.0 && p.Z() > -900.0) {
        ++nValidXYZTrj;
      }
    }
  }
  const AnaTrueParticlePD* beamTrue = nullptr;
  if (beampart && beampart->TrueObject) {
    const AnaTrueParticlePD* tpart = static_cast<const AnaTrueParticlePD*>(beampart->TrueObject);
    if (tpart) {
      trueIdBeamParticle = tpart->ID;
      beamTrue = tpart;
    }
  }
  output.FillVar(mainispandora, isPandora);
  output.FillVar(maintruepdg, truePdg);
  output.FillVar(maintruestartmomentum, pStartTrue);
  output.FillVar(maintrueendmomentum, pEndTrue);
  output.FillVar(maintruelength, trueLength);
  output.FillVar(maintrueid, trueIdMainTrack);
  output.FillVar(beamtrueid, trueIdBeamParticle);
  output.FillVar(mainlength, recoLength);
  output.FillVar(maintrueendtorecoendcm, trueEndToRecoEndCm);
  output.FillVar(mainnhitscollection, nHitsCollection);
  output.FillVar(mainntrjpoints, nTrjPoints);
  output.FillVar(mainnvalidrrhits, nValidRRHits);
  output.FillVar(mainnvalidxyzhits, nValidXYZHits);
  output.FillVar(mainnvalidxyztrj, nValidXYZTrj);
  output.FillVar(mainntruetrajpoints, mainntruetrajpointsIn);
  output.FillVar(mainmcsmomentum, std::isfinite(mainmcsmomentumGeV) ? static_cast<Float_t>(mainmcsmomentumGeV) : -999.f);
  output.FillVar(maintruemcsmomentum,
                 std::isfinite(maintruemcsmomentumGeV) ? static_cast<Float_t>(maintruemcsmomentumGeV) : -999.f);
  output.FillVar(maintlefitmomentum,
                 std::isfinite(maintlefitmomentumGeV) ? static_cast<Float_t>(maintlefitmomentumGeV) : -999.f);

  const size_t npairs =
      std::min(mainmcsdeltathetaPair.size(),
               std::min(mainmcsdeltathetaxzPair.size(), mainmcsdeltathetayzPair.size()));
  const size_t nfillPairs = std::min(npairs, nmaxmcs);
  for (size_t i = 0; i < nfillPairs; ++i) {
    output.FillVectorVar(mainmcsdeltatheta, static_cast<Float_t>(mainmcsdeltathetaPair[i]));
    output.FillVectorVar(mainmcsdeltathetaxz, static_cast<Float_t>(mainmcsdeltathetaxzPair[i]));
    output.FillVectorVar(mainmcsdeltathetayz, static_cast<Float_t>(mainmcsdeltathetayzPair[i]));
    output.IncrementCounterForVar(mainmcsdeltatheta);
  }

  const size_t nsegments =
      std::min(mainmcssegmentlengthSingle.size(),
               std::min(mainmcssegmentfitchi2ndfSingle.size(), mainmcssegmentrrtoendSingle.size()));
  const size_t nfillSeg = std::min(nsegments, nmaxmcs);
  for (size_t i = 0; i < nfillSeg; ++i) {
    output.FillVectorVar(mainmcssegmentlength, static_cast<Float_t>(mainmcssegmentlengthSingle[i]));
    output.FillVectorVar(mainmcssegmentfitchi2ndf, static_cast<Float_t>(mainmcssegmentfitchi2ndfSingle[i]));
    output.FillVectorVar(mainmcssegmentrrtoend, static_cast<Float_t>(mainmcssegmentrrtoendSingle[i]));
    output.IncrementCounterForVar(mainmcssegmentlength);
  }

  const size_t npairsTrue =
      std::min(maintruemcsdeltathetaPair.size(),
               std::min(maintruemcsdeltathetaxzPair.size(), maintruemcsdeltathetayzPair.size()));
  const size_t nfillPairsTrue = std::min(npairsTrue, nmaxmcs);
  for (size_t i = 0; i < nfillPairsTrue; ++i) {
    output.FillVectorVar(maintruemcsdeltatheta, static_cast<Float_t>(maintruemcsdeltathetaPair[i]));
    output.FillVectorVar(maintruemcsdeltathetaxz, static_cast<Float_t>(maintruemcsdeltathetaxzPair[i]));
    output.FillVectorVar(maintruemcsdeltathetayz, static_cast<Float_t>(maintruemcsdeltathetayzPair[i]));
    output.IncrementCounterForVar(maintruemcsdeltatheta);
  }

  const size_t nsegmentsTrue =
      std::min(maintruemcssegmentlengthSingle.size(),
               std::min(maintruemcssegmentfitchi2ndfSingle.size(), maintruemcssegmentrrtoendSingle.size()));
  const size_t nfillSegTrue = std::min(nsegmentsTrue, nmaxmcs);
  for (size_t i = 0; i < nfillSegTrue; ++i) {
    output.FillVectorVar(maintruemcssegmentlength, static_cast<Float_t>(maintruemcssegmentlengthSingle[i]));
    output.FillVectorVar(maintruemcssegmentfitchi2ndf, static_cast<Float_t>(maintruemcssegmentfitchi2ndfSingle[i]));
    output.FillVectorVar(maintruemcssegmentrrtoend, static_cast<Float_t>(maintruemcssegmentrrtoendSingle[i]));
    output.IncrementCounterForVar(maintruemcssegmentlength);
  }

  const UInt_t nmaxelastic = 512;
  const UInt_t nmaxproc = 512;
  const UInt_t nmaxincident = 4096;
  if (beamTrue) {
    output.FillVar(beamtruelaststeplen, beamTrue->TrueBeamLastStepLen);
    output.FillVar(beamtruenelastic, beamTrue->TrueBeamNElasticScatters);
    output.FillVar(beamtrueidetotaldep, beamTrue->TrueBeamIDEtotalDep);
    output.FillVar(beamtruenhits, beamTrue->TrueBeamNHits);
    const size_t nc = std::min(
        {beamTrue->TrueBeamElasticCosTheta.size(), beamTrue->TrueBeamElasticX.size(),
         beamTrue->TrueBeamElasticY.size(), beamTrue->TrueBeamElasticZ.size(),
         beamTrue->TrueBeamElasticDeltaE.size(), beamTrue->TrueBeamElasticIDEedep.size()});
    const size_t nfillE = std::min(nc, static_cast<size_t>(nmaxelastic));
    for (size_t i = 0; i < nfillE; ++i) {
      output.FillVectorVar(beamtrueelasticcostheta, beamTrue->TrueBeamElasticCosTheta[i]);
      output.FillVectorVar(beamtrueelasticx, beamTrue->TrueBeamElasticX[i]);
      output.FillVectorVar(beamtrueelasticy, beamTrue->TrueBeamElasticY[i]);
      output.FillVectorVar(beamtrueelasticz, beamTrue->TrueBeamElasticZ[i]);
      output.FillVectorVar(beamtrueelasticdeltae, beamTrue->TrueBeamElasticDeltaE[i]);
      output.FillVectorVar(beamtrueelasticideedep, beamTrue->TrueBeamElasticIDEedep[i]);
      output.IncrementCounterForVar(beamtrueelasticcostheta);
    }
    const size_t np = std::min(beamTrue->TrueBeamProcesses.size(), static_cast<size_t>(nmaxproc));
    for (size_t i = 0; i < np; ++i) {
      output.FillVectorVar(beamtrueprocesses, beamTrue->TrueBeamProcesses[i]);
      output.IncrementCounterForVar(beamtrueprocesses);
    }
    const size_t ni =
        std::min(beamTrue->TrueBeamIncidentEnergies.size(), static_cast<size_t>(nmaxincident));
    for (size_t i = 0; i < ni; ++i) {
      output.FillVectorVar(beamtrueincidentenergies, beamTrue->TrueBeamIncidentEnergies[i]);
      output.IncrementCounterForVar(beamtrueincidentenergies);
    }
  } else {
    output.FillVar(beamtruelaststeplen, -999.);
    output.FillVar(beamtruenelastic, -999);
    output.FillVar(beamtrueidetotaldep, -999.);
    output.FillVar(beamtruenhits, -999);
  }

  output.FillVar(mainbraggchi2pibb,
                 std::isfinite(mainbraggchi2pibbIn) ? static_cast<Float_t>(mainbraggchi2pibbIn) : -999.f);
  output.FillVar(mainbraggdedxnhits, mainbraggdedxnhitsIn);
}

//********************************************************************
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
                                                    const std::vector<Int_t>& braggNhits) {
//********************************************************************
  constexpr size_t nmaxdau = 64;
  const size_t n = std::min(
      {trueIds.size(), truePdgs.size(), trueStartMomGeV.size(), trueEndMomGeV.size(), trueStartEkinGeV.size(),
       trueEndEkinGeV.size(), csdaMomentumGeV.size(), csdaKineticEnergyGeV.size(), trueEndProcess.size(),
       trueEndToRecoStartDistCm.size(), trueEndToRecoEndDistCm.size(), braggChi2MeanPerHit.size(), braggNhits.size(),
       nmaxdau});
  for (size_t i = 0; i < n; ++i) {
    output.FillVectorVar(beamdaughtertrueid, trueIds[i]);
    output.FillVectorVar(beamdaughtertruepdg, truePdgs[i]);
    output.FillVectorVar(beamdaughtertruestartmom, trueStartMomGeV[i]);
    output.FillVectorVar(beamdaughtertrueendmom, trueEndMomGeV[i]);
    output.FillVectorVar(beamdaughtertruestartekin, trueStartEkinGeV[i]);
    output.FillVectorVar(beamdaughtertrueendekin, trueEndEkinGeV[i]);
    output.FillVectorVar(beamdaughtercsdamomentum, csdaMomentumGeV[i]);
    output.FillVectorVar(beamdaughtercsdakineticenergy, csdaKineticEnergyGeV[i]);
    output.FillVectorVar(beamdaughtertrueendprocess, trueEndProcess[i]);
    output.FillVectorVar(beamdaughtertrueendtorecostartcm, trueEndToRecoStartDistCm[i]);
    output.FillVectorVar(beamdaughtertrueendtorecoendcm, trueEndToRecoEndDistCm[i]);
    output.FillVectorVar(beamdaughterbraggchi2pibb, static_cast<Float_t>(braggChi2MeanPerHit[i]));
    output.FillVectorVar(beamdaughterbraggdedxnhits, braggNhits[i]);
    output.IncrementCounterForVar(beamdaughtertrueid);
  }
}

//********************************************************************
void FillPionMomentumVariables_BeamDaughterTleTruncScan(
    OutputManager& output, const std::vector<Int_t>& truncDauIdx, const std::vector<Int_t>& truncK,
    const std::vector<Int_t>& truncNhitsInt, const std::vector<Float_t>& truncTrueEkin0GeV,
    const std::vector<Float_t>& truncEkinCsdaFullGeV, const std::vector<Float_t>& truncPtleGeV,
    const std::vector<Float_t>& truncFracResTrueRef, const std::vector<Float_t>& truncFracResCsdaRef) {
//********************************************************************
  constexpr size_t nmaxrows = 4096;
  const size_t n = std::min({truncDauIdx.size(), truncK.size(), truncNhitsInt.size(), truncTrueEkin0GeV.size(),
                             truncEkinCsdaFullGeV.size(), truncPtleGeV.size(), truncFracResTrueRef.size(),
                             truncFracResCsdaRef.size(), nmaxrows});
  for (size_t i = 0; i < n; ++i) {
    output.FillVectorVar(beamdaughtertletruncdauidx, truncDauIdx[i]);
    output.FillVectorVar(beamdaughtertletrunck, truncK[i]);
    output.FillVectorVar(beamdaughtertletruncnhitsint, truncNhitsInt[i]);
    output.FillVectorVar(beamdaughtertletrunctrueekin0, truncTrueEkin0GeV[i]);
    output.FillVectorVar(beamdaughtertletruncekincsdafull, truncEkinCsdaFullGeV[i]);
    output.FillVectorVar(beamdaughtertletruncptle, truncPtleGeV[i]);
    output.FillVectorVar(beamdaughtertletruncfracres, truncFracResTrueRef[i]);
    output.FillVectorVar(beamdaughtertletruncfracrescsda, truncFracResCsdaRef[i]);
    output.IncrementCounterForVar(beamdaughtertletruncdauidx);
  }
}

//********************************************************************
void FillPionMomentumVariables_BeamDaughterMcsTruncScan(
    OutputManager& output, const std::vector<Int_t>& truncDauIdx, const std::vector<Int_t>& truncK,
    const std::vector<Int_t>& truncNsegments, const std::vector<Float_t>& truncTrueEkin0GeV,
    const std::vector<Float_t>& truncEkinCsdaFullGeV, const std::vector<Float_t>& truncPmcsGeV,
    const std::vector<Float_t>& truncFracResTrueRef, const std::vector<Float_t>& truncFracResCsdaRef) {
//********************************************************************
  constexpr size_t nmaxrows = 4096;
  const size_t n = std::min({truncDauIdx.size(), truncK.size(), truncNsegments.size(), truncTrueEkin0GeV.size(),
                             truncEkinCsdaFullGeV.size(), truncPmcsGeV.size(), truncFracResTrueRef.size(),
                             truncFracResCsdaRef.size(), nmaxrows});
  for (size_t i = 0; i < n; ++i) {
    output.FillVectorVar(beamdaughtermcstruncdauidx, truncDauIdx[i]);
    output.FillVectorVar(beamdaughtermcstrunck, truncK[i]);
    output.FillVectorVar(beamdaughtermcstruncnsegments, truncNsegments[i]);
    output.FillVectorVar(beamdaughtermcstrunctrueekin0, truncTrueEkin0GeV[i]);
    output.FillVectorVar(beamdaughtermcstruncekincsdafull, truncEkinCsdaFullGeV[i]);
    output.FillVectorVar(beamdaughtermcstruncpmcs, truncPmcsGeV[i]);
    output.FillVectorVar(beamdaughtermcstruncfracres, truncFracResTrueRef[i]);
    output.FillVectorVar(beamdaughtermcstruncfracrescsda, truncFracResCsdaRef[i]);
    output.IncrementCounterForVar(beamdaughtermcstruncdauidx);
  }
}

}  // namespace pionMomentumTree
