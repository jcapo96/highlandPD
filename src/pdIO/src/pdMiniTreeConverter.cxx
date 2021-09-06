#include "pdMiniTreeConverter.hxx"



//*****************************************************************************
Int_t pdMiniTreeConverter::GetSpill(Long64_t& entry, AnaSpillC*& spill){
//*****************************************************************************

  Int_t status = MiniTreeConverter::GetSpill(entry, spill);
  
  
  for (size_t b=0;b<spill->Bunches.size();b++){
    AnaBunchPD* bunch = static_cast<AnaBunchPD*>(spill->Bunches[b]);
    for (UInt_t i=0;i<bunch->Particles.size();i++){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[i]);
      
      for (UInt_t j=0;j<part->Hits[2].size();j++){      
        AnaHitPD& hit = part->Hits[2][j];
        
        hit.Signal.clear();
        for (UInt_t w =0;w<bunch->CNNwires.size();w++){
          if ((UInt_t)(bunch->CNNwires[w].wire) == hit.Channel && (UInt_t)(bunch->CNNwires[w].time) == hit.StartTick ){
            for (UInt_t t=hit.StartTick;t<=hit.EndTick;t++){
              hit.Signal.push_back(bunch->CNNwires[w].adcs[t-bunch->CNNwires[w].time]);            
            }
          }
        }      
      }
      
    }
    
  }
  return status;  
}
