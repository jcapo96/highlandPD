#include "CNNUtils.hxx"

//*******************************************************
void CNNUtils::produce(std::vector<AnaHitPD>& hits){
//*******************************************************  

  std::vector< std::pair<unsigned int, float> > points;
  for (auto const & hit : hits){                  
    points.emplace_back(hit.WireID.Wire, hit.PeakTime);
  }
  auto batch_out = fPointIdAlg.predictIdVectors(points);
  if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: #points = "<<  points.size() << " batch_out.size() = "  << batch_out.size() << std::endl;
  if (points.size() != batch_out.size()){
      //                        throw cet::exception("EmTrackMichelId") << "hits processing failed" << std::endl;
      std::cout << "hits processing failed: " << points.size() << " " << batch_out.size()  << std::endl;
  }
}
  
//*******************************************************
void CNNUtils::FillHits(const AnaEventPD& evt, std::vector<AnaHitPD*>& hitPtrList){
//*******************************************************

  
  if (debug_level>=0) std::cout << "EmTrackMichelId.cxx: #parts "<< evt.nParticles << std::endl;  
  for (Int_t i=0;i<evt.nParticles/2+1;i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(evt.Particles[i]);

    std::cout << "ANSELMO: particle " << i << " --> #hits = " << part->Hits[2].size() << std::endl;       
    
    if (debug_level>=0){
      std::cout << "EmTrackMichelId.cxx: hits list for this particle: " << std::endl;       
      for (UInt_t j=0;j<part->Hits[2].size();j++){

        AnaHitPD& hit = part->Hits[2][j];
        
        std::cout << " - " << j << ": wire, time, ampl, integral =  "
                  << hit.WireID.Wire << " "
                  << hit.PeakTime << " "
                  << hit.PeakAmplitude << " "
                  << hit.Integral << " "
                  << hit.Channel << " "
                  << hit.StartTick << "-"
                  << hit.EndTick << " "
          //                  << evt.ADC[hit.Channel][hit.StartTick] << " "
                  << hit.Signal.size() << " " 
                  << std::endl;       
      }
    }
    
    for (UInt_t j=0;j<part->Hits[2].size();j++){

      AnaHitPD& hit = part->Hits[2][j];

      // Fill hits with some "meaningful"  adc values
      /*
      for (size_t k=hit.fPeakTime-5;k<hit.fPeakTime+5;k++){
        hit.fSignal.push_back(hit.fPeakAmplitude*10);
      }
      */
      fPointIdAlg.setWireDriftDataFromHit(hit);
      //      hitListHandle.push_back(&(hit));
      hitPtrList.push_back(&(hit));
    }
  }
}

