#include "ParticlePositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include <cassert>


//********************************************************************
ParticlePositionSCECorrection::ParticlePositionSCECorrection(){
//********************************************************************

  _sce = new SpaceCharge();
  _sce->Initialize();
  _sce0 = NULL;
  // SCE systematic (k0s_sce build): scale the calibration displacement map by SCEScale (HighLAND SCE.dat sigma 0.08).
  _scale = ND::params().GetParameterD("neutralKaonAnalysis.SCEScale");
  if(_scale != 1.0){
    std::cout << "ParticlePositionSCECorrection: SCE map scaled by " << _scale << std::endl;
    _sce->ApplyGlobalVariation(_scale);
    _sce0 = new SpaceCharge();
    _sce0->Initialize();
  }
}

//********************************************************************
void ParticlePositionSCECorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);
  const double xBiasCorrectionCm = ND::params().GetParameterD("neutralKaonAnalysis.VertexXBiasCorrectionCm");

  // Loop over particles
  for(UInt_t ipart = 0; ipart < bunch->Particles.size(); ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);

    if (!original) continue; //?

    int ntps = part->TrjPoints.size();
    //if no trajectory points correct only position start end
    if(ntps<1){
      _sce->ApplyParticlePositionCorrection(part);
    }
    else{
      //correct trajectory points
      _sce->ApplyTrjPointPositionCorrection(part);
      _sce->ApplyTrjPointDirectionCorrection(part);

      //correct length
      part->Length = pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(part);

      //correct end and start position/direction
      pdAnaUtils::ComputeParticlePositionAndDirection(part);
    }

    // SCE systematic: hits are not moved by the correction, so shift them by (scale-1) times the nominal displacement.
    if(_sce0){
      for (int plane = 0; plane < 3; ++plane) {
        for (size_t ihit = 0; ihit < part->Hits[plane].size(); ++ihit) {
          AnaHitPD& hit = part->Hits[plane][ihit];
          if (hit.Position.X() < -900 || hit.Position_NoSCE.X() < -900) continue;
          TVector3 d = _sce0->GetCalPosOffsets(hit.Position_NoSCE, hit.TPCid);
          hit.Position.SetX(hit.Position.X() - (_scale - 1.0) * d.X());
          hit.Position.SetY(hit.Position.Y() + (_scale - 1.0) * d.Y());
          hit.Position.SetZ(hit.Position.Z() + (_scale - 1.0) * d.Z());
        }
      }
    }

    // Apply global X-bias shift consistently to all reconstructed particle positions.
    if (xBiasCorrectionCm != 0.0) {
      part->PositionStart[0] -= xBiasCorrectionCm;
      part->PositionEnd[0] -= xBiasCorrectionCm;
      part->PositionStartSCE[0] -= xBiasCorrectionCm;
      part->PositionEndSCE[0] -= xBiasCorrectionCm;

      for (int plane = 0; plane < 3; ++plane) {
        for (size_t ihit = 0; ihit < part->Hits[plane].size(); ++ihit) {
          AnaHitPD& hit = part->Hits[plane][ihit];
          if (hit.Position.X() > -900) {
            hit.Position.SetX(hit.Position.X() - xBiasCorrectionCm);
          }
          if (hit.Position_NoSCE.X() > -900) {
            hit.Position_NoSCE.SetX(hit.Position_NoSCE.X() - xBiasCorrectionCm);
          }
        }
      }

      for (size_t itp = 0; itp < part->TrjPoints.size(); ++itp) {
        AnaTrajectoryPointPD& trj = part->TrjPoints[itp];
        if (trj.Position.X() > -900) {
          trj.Position.SetX(trj.Position.X() - xBiasCorrectionCm);
        }
        if (trj.Position_NoSCE.X() > -900) {
          trj.Position_NoSCE.SetX(trj.Position_NoSCE.X() - xBiasCorrectionCm);
        }
      }
    }
  }

  for (UInt_t it = 0; it < spill.TrueParticles.size(); ++it) {
    AnaTrueParticlePD* truePart = dynamic_cast<AnaTrueParticlePD*>(spill.TrueParticles[it]);
    if (!truePart || truePart->TrjPoints.empty()) continue;
    _sce->ApplyTrjPointPositionCorrection(truePart);
    _sce->ApplyTrjPointDirectionCorrection(truePart);

    if (xBiasCorrectionCm != 0.0) {
      for (size_t itp = 0; itp < truePart->TrjPoints.size(); ++itp) {
        AnaTrueTrajectoryPointPD& trj = truePart->TrjPoints[itp];
        if (trj.Position.X() > -900) {
          trj.Position.SetX(trj.Position.X() - xBiasCorrectionCm);
        }
        if (trj.Position_NoSCE.X() > -900) {
          trj.Position_NoSCE.SetX(trj.Position_NoSCE.X() - xBiasCorrectionCm);
        }
      }
    }
  }
}

