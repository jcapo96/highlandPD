#include "DaughterFinderCorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>

int _n_reco = 0;
int _n_abs = 0;
int _n_dis_below_2 = 0;


//********************************************************************
DaughterFinderCorrection::DaughterFinderCorrection(){
//********************************************************************

}

//********************************************************************
void DaughterFinderCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  // Loop over particles
  for(UInt_t ipart = 0; ipart < bunch->Particles.size(); ipart++){
    //get beam track
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    if(!part->isPandora)continue;
    //loop over daughters
    for(int idau = 0; idau < (int)part->Daughters.size(); idau++){
      //is it a kaon candidate?
      AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[idau]);
      std::pair<double,int> pid = pdAnaUtils::Chi2PID(*dau,321);
      double chi = pid.first/pid.second;
      double trunc = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dau->Hits[2]);
      if(!(abs(chi-2.5)<2.5 && abs(trunc-3)<1))continue;
      
      //has it zero daughters?
      if(dau->Daughters.size()!=0)continue;
      //get true particle
      AnaTrueParticlePD* truedau = static_cast<AnaTrueParticlePD*>(dau->TrueObject);
      
      //is it a kaon?
      if(truedau->PDG!=321)continue;
      //is there a muon among its daughters?
      for(int itdau = 0; itdau < (int)truedau->Daughters.size(); itdau++){
	AnaTrueParticlePD* truedau_truedau = pdAnaUtils::GetTrueParticle(spill.TrueParticles,truedau->Daughters[itdau]);
	if(!truedau_truedau)continue;
	if(abs(truedau_truedau->PDG)!=13)continue;
	AnaParticlePD* muon = pdAnaUtils::GetRecoParticleWithAssociatedTrueID(bunch->Particles,truedau_truedau->ID);
	if(!muon)_n_abs++;
	else{     
	  _n_reco++;
	  double dis = 0;
	  for(int i = 0; i < 3; i++)dis += pow(muon->PositionStart[i]-dau->PositionEnd[i],2);
	  dis = sqrt(dis);
	  if(dis<2)_n_dis_below_2++;
	}
	std::cout << _n_reco << " (" << _n_dis_below_2 << "), " << _n_abs << std::endl;
	
      }
    }
    break; //only one beam particle
  }
}

