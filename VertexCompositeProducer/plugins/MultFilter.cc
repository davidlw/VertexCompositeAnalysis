// -*- C++ -*-
//
// Package:    VertexCompositeAnalysis/VertexCompositeProducer
// Class:      MultFilter
// 
/**\class MultFilter MultFilter.cc VertexCompositeAnalysis/VertexCompositeProducer/plugins/MultFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yousen Zhang
//         Created:  Thu, 12 Nov 2020 05:06:49 GMT
//
//


// system include files
#include <memory>
#include <limits>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

//
// class declaration
//

class MultFilter : public edm::stream::EDFilter<> {
   public:
      explicit MultFilter(const edm::ParameterSet&);
      ~MultFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
      const edm::EDGetTokenT<edm::ValueMap<int> > tok_nTracksVMap_;
      edm::ValueMap<int> nTracksVMap_;
      int nMultMin_, nMultMax_;

      bool useCent_, vtxSortByTrkSize_ ;
      const edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
      const edm::EDGetTokenT<int> tok_centBin_;
      int centBinMin_, centBinMax_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MultFilter::MultFilter(const edm::ParameterSet& iConfig) :
  tok_offlinePV_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  tok_nTracksVMap_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("nTracksVMap"))),
  nMultMin_(iConfig.getUntrackedParameter<int>("nMultMin", 0)),
  nMultMax_(iConfig.getUntrackedParameter<int>("nMultMax", std::numeric_limits<int>::max())),
  useCent_(iConfig.getUntrackedParameter<bool>("useCent", false)),
  vtxSortByTrkSize_(iConfig.getUntrackedParameter<bool>("vtxSortByTrkSize", true)),
  tok_centSrc_(consumes<reco::Centrality>(iConfig.getUntrackedParameter<edm::InputTag>("centrality"))),
  tok_centBin_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBin"))),
  centBinMin_(iConfig.getUntrackedParameter<int>("centBinMin", 0)),
  centBinMax_(iConfig.getUntrackedParameter<int>("centBinMax", std::numeric_limits<int>::max()))
{
   //now do what ever initialization is needed
}


MultFilter::~MultFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MultFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(tok_offlinePV_, vertexHandle);

  if (!vtxSortByTrkSize_ && useCent_)
    throw(cms::Exception("Config")
        << "hicentrality producer make use of the vertex with the highest multiplicity, not the one with leading pt2!");

  if (!useCent_) {
    // do not use centrality
    edm::Handle<edm::ValueMap<int> > nTracksVMap;
    iEvent.getByToken(tok_nTracksVMap_, nTracksVMap);

    if (nTracksVMap.isValid() ) { nTracksVMap_ = *nTracksVMap; }
    else {
      LogDebug("MultFilter") <<  "nTracksVMap is invalid" << std::endl;;
      return false;
    }

    if (vertexHandle.isValid() && vertexHandle->size()) {
      std::vector<reco::VertexRef> vertexRefs;
      vertexRefs.reserve( vertexHandle->size() );
      for (size_t ivtx=0; ivtx<vertexHandle->size(); ivtx++) {
        vertexRefs.push_back( reco::VertexRef( vertexHandle, ivtx ) );
      }
      if (vtxSortByTrkSize_) {
        auto byTracksSize = [] (const reco::VertexRef& v1, const reco::VertexRef& v2) -> bool { return v1->tracksSize() > v2->tracksSize(); };
        std::sort(vertexRefs.begin(), vertexRefs.end(), byTracksSize);
      }
      const auto& ntrk = nTracksVMap_[ vertexRefs.at(0) ];
      if (ntrk < nMultMin_ || ntrk >= nMultMax_) {
        LogDebug("MultFilter") << std::boolalpha
          << "UseCent=False, vtxSortByTrkSize="
          << vtxSortByTrkSize_ << ", Ntrk: " << ntrk
          << ", nMultMin: "<< nMultMin_ << ", nMultMax: " << nMultMax_
          << ", Failed"<< std::endl;
        return false;
      }
      LogDebug("MultFilter") << std::boolalpha
        << "UseCent=False, vtxSortByTrkSize="
        << vtxSortByTrkSize_ << ", Ntrk: " << ntrk
        << ", nMultMin: "<< nMultMin_ << ", nMultMax: " << nMultMax_
        << ", pass"<< std::endl;
    } else {
      return false;
    }
  }
  else {
    // use centrality
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);
    edm::Handle<int> centBin;
    iEvent.getByToken(tok_centBin_, centBin);
    const auto ntrk = cent->Ntracks();
    if (ntrk < nMultMin_ || ntrk >= nMultMax_) {
      LogDebug("MultFilter") << "UseCent=True, Ntrk: " << ntrk
        << ", nMultMin: "<< nMultMin_ << ", nMultMax: " << nMultMax_
        << ", Failed"<< std::endl;
      return false;
    }
    if (*centBin < centBinMin_ || *centBin >= centBinMax_) {
      LogDebug("MultFilter") << "UseCent=True, centBin: " << *centBin
        << ", centBinMin: "<< centBinMin_ << ", centBinMax: " << centBinMax_
        << ", Failed." << std::endl;
      return false;
    }
    LogDebug("MultFilter") << "UseCent=True, Ntrk: " << ntrk
      << ", nMultMin: "<< nMultMin_ << ", nMultMax: " << nMultMax_
      << ", centBin: " << *centBin
      << ", centBinMin: "<< centBinMin_ << ", centBinMax: " << centBinMax_
      << ", Passed." << std::endl;
  }

  return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MultFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MultFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MultFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MultFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MultFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MultFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MultFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"))->setComment("primary vertex collection");
  desc.add<edm::InputTag>("nTracksVMap", edm::InputTag("generalParticles", "nTracks"))->setComment("Ntrk ValueMap");
  desc.addUntracked<int>("nMultMin", 0)->setComment(">=");
  desc.addUntracked<int>("nMultMax", std::numeric_limits<int>::max())->setComment("<");
  desc.addUntracked<bool>("useCent", false);
  desc.addUntracked<bool>("vtxSortByTrkSize", true);
  desc.addUntracked<edm::InputTag>("centrality", edm::InputTag("hiCentrality"));
  desc.addUntracked<edm::InputTag>("centralityBin", edm::InputTag("centralityBin","HFtowers"));
  desc.addUntracked<int>("centBinMin", 0)->setComment(">=");
  desc.addUntracked<int>("centBinMax", std::numeric_limits<int>::max())->setComment("<");
  //descriptions.addDefault(desc);
  descriptions.add("ntrkFilter", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MultFilter);
