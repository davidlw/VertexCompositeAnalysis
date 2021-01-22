// -*- C++ -*-
//
// Package:    NTrackVertexMapper
// Class:      NTrackVertexMapper
// 
/**\class NTrackVertexMapper NTrackVertexMapper.cc VertexCompositeAnalysis/VertexCompositeProducer/plugin/NTrackVertexMapper.cc
//
*/
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"


class NTrackVertexMapper : public edm::stream::EDProducer<>
{
 public:
  explicit NTrackVertexMapper(const edm::ParameterSet & iConfig):
    token_vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    token_tracks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
  {
    produces<edm::ValueMap<int> >();
  };
  ~NTrackVertexMapper() override {};
  void produce(edm::Event&, const edm::EventSetup&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

 private:
  const edm::EDGetTokenT<reco::VertexCollection> token_vertices_;
  const edm::EDGetTokenT<reco::TrackCollection> token_tracks_;
};


void NTrackVertexMapper::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // extract the input collections
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(token_vertices_, vertexHandle);
  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByToken(token_tracks_, trackHandle);

  // initialise container
  std::vector<int> vtxNTrk((vertexHandle.isValid() ? vertexHandle->size() : 0), 0);

  // loop over general tracks
  if (vertexHandle.isValid() && trackHandle.isValid() && vertexHandle->size()) {
    for (const auto& trk : *trackHandle) {
      if (!trk.quality(reco::TrackBase::highPurity)) continue;
      if (trk.pt() <= 0.4 || std::abs(trk.eta()) >= 2.4) continue;
      if (trk.ptError()/trk.pt() >= 0.1) continue;
      // loop over primary vertices
      for (size_t i=0; i<vertexHandle->size(); i++) {
        const auto& pv = vertexHandle->at(i);
        const auto dz = trk.dz(pv.position());
        const auto dzErr2 = trk.dzError()*trk.dzError() + pv.zError()*pv.zError();
        if (dz*dz >= 9.0*dzErr2) continue;
        const auto dxy = trk.dxy(pv.position());
        const auto dxyErr2 = trk.dxyError()*trk.dxyError() + pv.xError()*pv.yError();
        if (dxy*dxy >= 9.0*dxyErr2) continue;
        vtxNTrk[i] += 1;
      }
    }
  }

  // store value map
  auto valueMap = std::make_unique<edm::ValueMap<int> >();
  edm::ValueMap<int>::Filler filler(*valueMap);
  filler.insert(vertexHandle, vtxNTrk.begin(), vtxNTrk.end());
  filler.fill();
  iEvent.put(std::move(valueMap));
};


void NTrackVertexMapper::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"))->setComment("primary vertex collection");
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))->setComment("track collection");
  descriptions.add("nTracks", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(NTrackVertexMapper);
