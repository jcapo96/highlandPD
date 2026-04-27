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
  AddVarMaxSizeVF(output, mainmcsdeltatheta, "Main track MCS delta theta [rad]", mainmcsnsegments, nmaxmcs);
  AddVarMaxSizeVF(output, mainmcssegmentlength, "Main track MCS mean segment length [cm]", mainmcsnsegments, nmaxmcs);
  AddVarF(output, mainmcsmomentum, "Main track MCS momentum estimate [GeV/c]");
  AddVarF(output, maintlefitmomentum, "Main track TLEFit momentum estimate [GeV/c]");
  AddVarI(output, mainnhitscollection, "Main track number of collection-plane hits");
  AddVarI(output, mainntrjpoints, "Main track number of trajectory points");
  AddVarI(output, mainnvalidrrhits, "Main track number of collection hits with valid residual range");
  AddVarI(output, mainnvalidxyzhits, "Main track collection hits with valid XYZ");
  AddVarI(output, mainnvalidxyztrj, "Main track trajectory points with valid XYZ");
}

//********************************************************************
void FillPionMomentumVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beampart,
                                                const std::vector<double>& mainthetascatter,
                                                const std::vector<double>& mainsegmentlength,
                                                double mainmcsmomentumGeV,
                                                double maintlefitmomentumGeV) {
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
  if (beampart && beampart->TrueObject) {
    const AnaTrueParticlePD* tpart = static_cast<const AnaTrueParticlePD*>(beampart->TrueObject);
    if (tpart) trueIdBeamParticle = tpart->ID;
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
  output.FillVar(mainmcsmomentum, std::isfinite(mainmcsmomentumGeV) ? static_cast<Float_t>(mainmcsmomentumGeV) : -999.f);
  output.FillVar(maintlefitmomentum,
                 std::isfinite(maintlefitmomentumGeV) ? static_cast<Float_t>(maintlefitmomentumGeV) : -999.f);

  const size_t npaired = std::min(mainthetascatter.size(), mainsegmentlength.size());
  const size_t nfill = std::min(npaired, nmaxmcs);
  for (size_t i = 0; i < nfill; ++i) {
    output.FillVectorVar(mainmcsdeltatheta, static_cast<Float_t>(mainthetascatter[i]));
    output.FillVectorVar(mainmcssegmentlength, static_cast<Float_t>(mainsegmentlength[i]));
    output.IncrementCounterForVar(mainmcsdeltatheta);
  }
}

}  // namespace pionMomentumTree
