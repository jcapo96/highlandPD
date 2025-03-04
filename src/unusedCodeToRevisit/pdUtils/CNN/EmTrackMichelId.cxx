#include "EmTrackMichelId.hxx"

/*nnet::EmTrackMichelId::EmTrackMichelId(EmTrackMichelId::Parameters const& config) :
        EDProducer{config},
        fBatchSize(config().BatchSize()),
        fPointIdAlg(config().PointIdAlg()),
        fMVAWriter(producesCollector(), "emtrkmichel"),
        fWireProducerLabel(config().WireLabel()),
        fHitModuleLabel(config().HitModuleLabel()),
        fClusterModuleLabel(config().ClusterModuleLabel()),
        fTrackModuleLabel(config().TrackModuleLabel()),
        fViews(config().Views()),

        fNewClustersTag(
            config.get_PSet().get<std::string>("module_label"), "",
            art::ServiceHandle<art::TriggerNamesService const>()->getProcessName())
{
    fMVAWriter.produces_using< my_recob::Hit >();

    if (!fClusterModuleLabel.label().empty())
    {
        produces< std::vector<my_recob::Cluster> >();
            produces< art::Assns<my_recob::Cluster, my_recob::Hit> >();

        fMVAWriter.produces_using< my_recob::Cluster >();
        fDoClusters = true;
    }
    else { fDoClusters = false; }

    if (!fTrackModuleLabel.label().empty())
    {
        fMVAWriter.produces_using< my_recob::Track >();
        fDoTracks = true;
    }
    else { fDoTracks = false; }
    }*/

nnet::EmTrackMichelId::EmTrackMichelId():
//        EDProducer{config},
        fBatchSize(1),
        //        fPointIdAlg(config().PointIdAlg()),
        //        fMVAWriter(producesCollector(), "emtrkmichel"),
        fWireProducerLabel("wclsdatasp:gaus"),
        fHitModuleLabel("hitpdune"),
        fClusterModuleLabel("pandora"),
        fTrackModuleLabel("dummy")  // anselmo
//        fViews(config().Views())  // anselmo
//        fNewClustersTag()
{


  fViews.clear();  //anselmo
  
  if (!fClusterModuleLabel.empty())fDoClusters = true;
  else fDoClusters = false;

  if (!fTrackModuleLabel.empty())fDoTracks = true;
  else fDoTracks = false;

  //fNewClustersTag;? // input tag for the clusters produced by this module

}
// ------------------------------------------------------

void nnet::EmTrackMichelId::produce(AnaEventPD & evt)
{

  /*
      - loop over Cryostats
        - loop over TPCs in that cryostat
          - loop over planes/views in that TPC
            - loop over batches of hits in that plane, with size fBatchSize
              - loop over hits in that batch
                - fill the points vector: points.emplace_back(hit.WireID().Wire, hit.PeakTime());
              - call TF for that batch of hits: fPointIdAlg.predictIdVectors(points);
                - loop over points (= pairs of wire-PeakTime)
                  - call PointIdAlg::bufferPatch
                    - fill a 2D vector with adc values in a window of wires and time around the nominal point : patch[wire][drifttime]. 
                      The function doing that is DataProviderAlg::patchFromOriginalView(wire, drift, fPatchSizeW, fPatchSizeD, patch);
                      The window for Original view is W = fPatchSizeW=44, D = fPatchSizeD*fDriftWindow = 48*6=288
                - Fill the inps 3D vector:  inps[hit index in batch][wire index][drifttime index]
                - call fNNet->Run(inps);
                  - Fills the input map for TF: input_map(s, r, c, 0) = row[c];  where s=hit index in batch, r=wire index, c=drifttime index, row[c]=adc value
                  - Call TF:  g->run(_x);

                  
      That means that the input tree should only contain adc values for a window of wires and times around the hit, for all hits in the tracks we are interested in. 
      The pion selection uses CNN for beam particle daughters (1.3 in average), which usually don't have many hits (avg=80)
      Data Reduction: 
       - For each hit 44 adc vectors (one for each wire) of 288 time samples are needed. However those 44 vectors are common to many hits, since for a particle those hits are contiguous.  
       - Not all 288 time samples have contents, thus one can probaly save a vector of pairs (time,adc) for the meaningful samples 

    */


  
  std::cout << "next event: " << evt.EventInfo.Run << " / " << evt.EventInfo.Event << std::endl;
  /*
  mf::LogVerbatim("EmTrackMichelId") << "next event: " << evt.run() << " / " << evt.id().event();

  auto wireHandle = evt.getValidHandle< std::vector<my_recob::Wire> >(fWireProducerLabel);
  */

  //  detinfo::DetectorClocksData clockData;
  //  detinfo::DetectorPropertiesData detProp;

  //  std::vector<my_recob::Wire> wireHandle;
  //  FillWires(evt,wireHandle);

  
  unsigned int cryo, tpc, view;
  
  // ******************* get and sort hits ********************

  //  std::vector<AnaHitPD*> hitListHandle;   // anselmo: not used
  std::vector<AnaHitPD*> hitPtrList;
  //  auto hitListHandle = evt.getValidHandle< std::vector<my_recob::Hit> >(fHitModuleLabel);
  //  art::fill_ptr_vector(hitPtrList, hitListHandle);

  FillHits(evt,hitPtrList);
     
  if (debug_level>=0) std::cout << "EmTrackMichelId.cxx: total #hits "<< hitPtrList.size() << std::endl; 

  EmTrackMichelId::cryo_tpc_view_keymap hitMap;
  for (auto const& h : hitPtrList)
    {
      view = h->WireID.Plane;
      if (!isViewSelected(view)) continue;
      
      cryo = h->WireID.Cryostat;
      tpc = h->WireID.TPC;
      
      //      hitMap[cryo][tpc][view].push_back(h.key());
      //      hitMap[cryo][tpc][view].push_back((size_t)h);
      // anselmo

      //      std::cout << h->fChannel << " " << cryo << " " << tpc << " " << view << " " << h->WireID().Wire << std::endl; 
      hitMap[cryo][tpc][view].push_back(h);
    }

    // ********************* classify hits **********************
  //    auto hitID = fMVAWriter.initOutputs<my_recob::Hit>(fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());  // anselmo. Only used for writter

    std::vector< char > hitInFA(hitPtrList.size(), 0); // tag hits in fid. area as 1, use 0 for hits close to the projectrion edges

    for (auto const & pcryo : hitMap)
    {
        cryo = pcryo.first;
        for (auto const & ptpc : pcryo.second)
        {
            tpc = ptpc.first;
            for (auto const & pview : ptpc.second)
            {
                view = pview.first;
                if (!isViewSelected(view)) continue; // should not happen, hits were selected

                if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: cryo, tpc, view = "<< cryo << " " << tpc << " " << view  << " --> #hits = " << pview.second.size() << std::endl;
                
                //                fPointIdAlg.setWireDriftData(clockData, detProp, wireHandle, view, tpc, cryo);  // anselmo


                // (1) do all hits in this plane ------------------------------------------------
                for (size_t idx = 0; idx < pview.second.size(); idx += fBatchSize)
                {
                  if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: processing hit with idx = "<< idx << std::endl;
                    std::vector< std::pair<unsigned int, float> > points;
                    std::vector< size_t > keys;
                    for (size_t k = 0; k < fBatchSize; ++k)
                    {

                        if (idx + k >= pview.second.size()) { break; } // careful about the tail
                        // anselmo:  use directly the pointer, not using the key
                        //                        size_t h = pview.second[idx+k]; // h is the Ptr< my_recob::Hit >::key()
                        size_t h = (size_t)(pview.second[idx+k]);
                        //                        const AnaHitPD & hit = *(hitPtrList[h]);
                        const AnaHitPD & hit = *(pview.second[idx+k]);
                        if (debug_level>=2) std::cout << "  EmTrackMichelId.cxx: idx+k  = "<< idx+k << ", hit (wireID, peak time) = " << hit.WireID.Wire << " " << hit.PeakTime << std::endl;
                  
                        points.emplace_back(hit.WireID.Wire, hit.PeakTime);
                        //                        keys.push_back(h);
                        keys.push_back((size_t)h);
                    }
                    auto batch_out = fPointIdAlg.predictIdVectors(points);
                    if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: #points = "<<  points.size() << " batch_out.size() = "  << batch_out.size() << std::endl;
                    if (points.size() != batch_out.size())
                    {
                      //                        throw cet::exception("EmTrackMichelId") << "hits processing failed" << std::endl;
                      std::cout << "hits processing failed: " << points.size() << " " << batch_out.size()  << std::endl;
                    }
                    for (size_t k = 0; k < points.size(); ++k)
                    {
                      //                        size_t h = keys[k];  // anselmo
                        //                        fMVAWriter.setOutput(hitID, h, batch_out[k]);   // anselmo
                        if (debug_level>=2) std::cout << "  EmTrackMichelId.cxx: points[k].first, points[k].second = "<<  points[k].first << " " << points[k].second << std::endl;
                        if (fPointIdAlg.isInsideFiducialRegion(points[k].first, points[k].second))
                          { /*hitInFA[h] = 1;     */}   // anselmo: must change h per hit index
                    }
                } // hits done ------------------------------------------------------------------
            }
        }
    }

  /*

  // (2) do clusters when hits are ready in all planes ----------------------------------------
    if (fDoClusters)
    {
        // **************** prepare for new clusters ****************
            auto clusters = std::make_unique< std::vector< my_recob::Cluster > >();
            auto clu2hit = std::make_unique< art::Assns< my_recob::Cluster, my_recob::Hit > >();

        // ************** get and sort input clusters ***************
        auto cluListHandle = evt.getValidHandle< std::vector<my_recob::Cluster> >(fClusterModuleLabel);
            std::vector< art::Ptr<my_recob::Cluster> > cluPtrList;
            art::fill_ptr_vector(cluPtrList, cluListHandle);

        EmTrackMichelId::cryo_tpc_view_keymap cluMap;
            for (auto const& c : cluPtrList)
            {
                view = c->Plane().Plane;
                if (!isViewSelected(view)) continue;

                cryo = c->Plane().Cryostat;
                tpc = c->Plane().TPC;

                cluMap[cryo][tpc][view].push_back(c.key());
            }

        auto cluID = fMVAWriter.initOutputs<my_recob::Cluster>(fNewClustersTag, fPointIdAlg.outputLabels());

        unsigned int cidx = 0; // new clusters index
        art::FindManyP< my_recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
        std::vector< bool > hitUsed(hitPtrList.size(), false); // tag hits used in clusters
        for (auto const & pcryo : cluMap)
        {
            cryo = pcryo.first;
            for (auto const & ptpc : pcryo.second)
            {
                tpc = ptpc.first;
                for (auto const & pview : ptpc.second)
                {
                    view = pview.first;
                    if (!isViewSelected(view)) continue; // should not happen, clusters were pre-selected

                    for (size_t c : pview.second) // c is the Ptr< my_recob::Cluster >::key()
                    {
                                auto v = hitsFromClusters.at(c);
                                if (v.empty()) continue;

                        for (auto const & hit : v)
                        {
                            if (hitUsed[hit.key()]) { mf::LogWarning("EmTrackMichelId") << "hit already used in another cluster"; }
                            hitUsed[hit.key()] = true;
                        }

                        auto vout = fMVAWriter.getOutput<my_recob::Hit>(v,
                            [&](art::Ptr<my_recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });

                            float pvalue = vout[0] / (vout[0] + vout[1]);
                            mf::LogVerbatim("EmTrackMichelId") << "cluster in tpc:" << tpc << " view:" << view
                            << " size:" << v.size() << " p:" << pvalue;

                                clusters->emplace_back(
                                    my_recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
                                            v.size(), 0.0F, 0.0F, cidx, (geo::View_t)view, v.front()->WireID().planeID()));
                            util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
                            cidx++;

                            fMVAWriter.addOutput(cluID, vout); // add copy of the input cluster
                    }

                    // (2b) make single-hit clusters --------------------------------------------
                    for (size_t h : hitMap[cryo][tpc][view]) // h is the Ptr< my_recob::Hit >::key()
                    {
                        if (hitUsed[h]) continue;

                        auto vout = fMVAWriter.getOutput<my_recob::Hit>(h);
                                float pvalue = vout[0] / (vout[0] + vout[1]);

                                mf::LogVerbatim("EmTrackMichelId") << "single hit in tpc:" << tpc << " view:" << view
                                        << " wire:" << hitPtrList[h]->WireID().Wire << " drift:" << hitPtrList[h]->PeakTime() << " p:" << pvalue;

                                art::PtrVector< my_recob::Hit > cluster_hits;
                                cluster_hits.push_back(hitPtrList[h]);
                                clusters->emplace_back(
                                        my_recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
                                            1, 0.0F, 0.0F, cidx, (geo::View_t)view, hitPtrList[h]->WireID().planeID()));
                                util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
                                cidx++;

                                fMVAWriter.addOutput(cluID, vout); // add single-hit cluster tagging unclutered hit
                    }
                    mf::LogVerbatim("EmTrackMichelId") << "...produced " << cidx - pview.second.size() << " single-hit clusters.";
                }
            }
        }

        evt.put(std::move(clusters));
            evt.put(std::move(clu2hit));
        } // all clusters done ----------------------------------------------------------------------
  */
    // (3) do tracks when all hits in all cryo/tpc/plane are done -------------------------------
    if (fDoTracks)
    {
      /*
        auto trkListHandle = evt.getValidHandle< std::vector<my_recob::Track> >(fTrackModuleLabel);
        art::FindManyP< my_recob::Hit > hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
        std::vector< std::vector< art::Ptr<my_recob::Hit> > > trkHitPtrList(trkListHandle->size());
      */

      // Anselmo
      std::vector<my_recob::Track>  trkListHandle;
      std::vector<std::vector<AnaHitPD> > hitsFromTracks;
      std::vector<std::vector<AnaHitPD> > trkHitPtrList;

      for (Int_t i=0;i<evt.nParticles;i++){
        AnaParticlePD* part = static_cast<AnaParticlePD*>(evt.Particles[i]);
        my_recob::Track track(*part);
        trkListHandle.push_back(track);
        hitsFromTracks.push_back(part->Hits[2]);        
      }
      trkHitPtrList.resize(trkListHandle.size());

      if (debug_level>=0) std::cout << "EmTrackMichelId.cxx: #tracks "<< trkListHandle.size() << std::endl; 
      
        for (size_t t = 0; t < trkListHandle.size(); ++t)
        {
            auto v = hitsFromTracks.at(t);
            size_t nh[3] = { 0, 0, 0 };
            if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: Track " << t << " --> #hits = "<< v.size() << std::endl; 
            for (auto const & hptr : v) { ++nh[hptr.View]; }  // count hits in eac view
            size_t best_view = 2; // collection
            if ((nh[0] >= nh[1]) && (nh[0] > 2 * nh[2])) best_view = 0; // ind1
            if ((nh[1] >= nh[0]) && (nh[1] > 2 * nh[2])) best_view = 1; // ind2
            size_t k = 0;

            if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: hits per view = " << nh[0] << " " << nh[1] << " " << nh[2] << std::endl; 
            if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: best view = " << best_view << std::endl; 
            while (!isViewSelected(best_view))
            {
                best_view = (best_view + 1) % 3;
                //                if (++k > 3) { throw cet::exception("EmTrackMichelId") << "No views selected at all?" << std::endl; }
                if (++k > 3) { std::cout << "No views selected at all?" << std::endl; }
            }

            if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: best view (new) = " << best_view << std::endl; 
            for (auto const & hptr : v)
            {
	      if ((UInt_t)hptr.View == best_view) trkHitPtrList[t].emplace_back(hptr);
            }

            if (debug_level>=1) std::cout << " EmTrackMichelId.cxx: trkHitPtrList[t].size() = " << trkHitPtrList[t].size() << std::endl; 
            
        }
        /*
        auto trkID = fMVAWriter.initOutputs<my_recob::Track>(fTrackModuleLabel, trkHitPtrList.size(), fPointIdAlg.outputLabels());
        for (size_t t = 0; t < trkHitPtrList.size(); ++t) // t is the Ptr< my_recob::Track >::key()
        {
            auto vout = fMVAWriter.getOutput<my_recob::Hit>(trkHitPtrList[t],
                [&](art::Ptr<my_recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });
            fMVAWriter.setOutput(trkID, t, vout);
        }
        */
    }
    // tracks done ------------------------------------------------------------------------------

    //        fMVAWriter.saveOutputs(evt);
}
// ------------------------------------------------------

bool nnet::EmTrackMichelId::isViewSelected(int view) const
{
  if (debug_level>=10) std::cout << "  EmTrackMichelId.cxx:isViewSelected. View = " << view << std::endl;   
  if (fViews.empty()) return true;
  else
    {
      bool selected = false;
      for (auto k : fViews) if (k == view) { selected = true; break; }
      return selected;
    }
}
// ------------------------------------------------------


//*******************************************************
void nnet::EmTrackMichelId::FillWires(const AnaEventPD& evt, std::vector<my_recob::Wire>& wires){
//*******************************************************
  
  // anselmo: fill the wire vector 15400 wires
  //  for (size_t i=0;i<15400;i++){
  /*
  for (size_t i=0;i<1;i++){
    my_recob::Wire wire;
    wire.fChannel=149;
    for (size_t j=0;j<6000;j++){
      if (j>4000 && j<5000)
        wire.fSignal.push_back(149.33);
      else
        wire.fSignal.push_back(0);
    }

    wireHandle.push_back(wire);
  }
  */

}



  
//*******************************************************
void nnet::EmTrackMichelId::FillHits(const AnaEventPD& evt, std::vector<AnaHitPD*>& hitPtrList){
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

//*******************************************************
void nnet::EmTrackMichelId::DumpADCs(const AnaEventPD& evt){
//*******************************************************
  
  std::cout << "ANSELMO: List of ADC values" << std::endl;
  for (UInt_t j=0;j<evt.CNNwires.size();j++){
    for (UInt_t k=0;k<evt.CNNwires[j].adcs.size();k++){
      if (evt.CNNwires[j].adcs[k]>0) 
        std::cout << j << " " << k << " " << evt.CNNwires[j].adcs[k] << std::endl;
    }
  }
}

