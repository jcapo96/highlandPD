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
          "Main track chi2_pi mean (arXiv:2409.18288 Eq. 6.1): mean [(dEdx-dEdx_BB)^2/sigma^2] for 0<RR<Rmax; smaller≈stopping");
  AddVarI(output, mainbraggdedxnhits, "Main track collection hits used for Bragg-window chi2_pi (Eq. 6.1)");
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

}  // namespace pionMomentumTree
