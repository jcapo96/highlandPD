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
}

