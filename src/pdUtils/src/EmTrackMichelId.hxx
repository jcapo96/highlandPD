/////////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackMichelId
// Module Type: producer
// File:        EmTrackMichelId_module.cc
// Authors:     D.Stefan (dorota.stefan@cern.ch), from DUNE, CERN/NCBJ
//              P.Plonski (pplonski86@gmail.com), from DUNE, WUT
//              R.Sulej (robert.sulej@cern.ch),   from DUNE, FNAL/NCBJ
//
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to distinguish EM-like activity from track-like objects. In addition the activity
// of Michel electrons is recognized. New clusters of hits are produced to include
// also unclustered hits and tag all the activity as EM/track in the same way.
// Note: Michel electrons are best tagged on the level of single hits, since
// clustering may have introduced mistakes (hits of muon and electron clustered
// together).
//
/////////////////////////////////////////////////////////////////////////////////////

#ifndef EmTrackMichelId_h
#define EmTrackMichelId_h

/*
  #include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"
#include "lardata/ArtDataHelper/MVAWriter.h"

#include <memory>
*/

#include "PointIdAlg.hxx"

namespace nnet {
  
  class EmTrackMichelId{// : public art::EDProducer {
  public:
    
    // these types to be replaced with use of feature proposed in redmine #12602
    typedef std::unordered_map< unsigned int, std::vector< size_t > > view_keymap;
    typedef std::unordered_map< unsigned int, view_keymap > tpc_view_keymap;
    typedef std::unordered_map< unsigned int, tpc_view_keymap > cryo_tpc_view_keymap;
    
    /*struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg {
      Name("PointIdAlg")
      };
                fhicl::Atom<size_t> BatchSize {
		Name("BatchSize"), Comment("number of samples processed in one batch")
                };
		
                fhicl::Atom<art::InputTag> WireLabel {
		Name("WireLabel"), Comment("tag of deconvoluted ADC on wires (recob::Wire)")
                };
		
                fhicl::Atom<art::InputTag> HitModuleLabel {
		Name("HitModuleLabel"), Comment("tag of hits to be EM/track / Michel tagged")
                };
		
                fhicl::Atom<art::InputTag> ClusterModuleLabel {
		Name("ClusterModuleLabel"),
		Comment("tag of clusters to be used as a source of EM/track / Michel tagged new clusters (incl. single-hit clusters ) using accumulated results from hits")
                };
		
                fhicl::Atom<art::InputTag> TrackModuleLabel {
		Name("TrackModuleLabel"),
		Comment("tag of 3D tracks to be EM/track / Michel tagged using accumulated results from hits in the best 2D projection")
                };
		
                fhicl::Sequence<int> Views {
		Name("Views"),
		Comment("tag clusters in selected views only, or in all views if empty list")
                };
		};*/
    //using Parameters = art::EDProducer::Table<Config>;
    //      explicit EmTrackMichelId(Parameters const & p);
    
    //      EmTrackMichelId(EmTrackMichelId const &) = delete;
    //      EmTrackMichelId(EmTrackMichelId &&) = delete;
    //      EmTrackMichelId & operator = (EmTrackMichelId const &) = delete;
    //      EmTrackMichelId & operator = (EmTrackMichelId &&) = delete;
    
    EmTrackMichelId();
    
  private:
    void produce(AnaEvent & e); //override;
    
    bool isViewSelected(int view) const;
    
    size_t fBatchSize;
    PointIdAlg fPointIdAlg;
    //anab::MVAWriter<4> fMVAWriter; // <-------------- using 4-output CNN model
    /*
      art::InputTag fWireProducerLabel;
      art::InputTag fHitModuleLabel;
      art::InputTag fClusterModuleLabel;
      art::InputTag fTrackModuleLabel;
    */
    
    std::string fWireProducerLabel;
    std::string fHitModuleLabel;
    std::string fClusterModuleLabel;
    std::string fTrackModuleLabel;
    bool fDoClusters, fDoTracks;
    
    std::vector< int > fViews{};
    
    std::string fNewClustersTag; // input tag for the clusters produced by this module
  };
}
 
  
#endif
