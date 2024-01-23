// -*- C++ -*-
//
// Package:    DeDxMuonProducer
// Class:      DeDxMuonProducer
// 
/**\class DeDxMuonProducer DeDxMuonProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/plugin/DeDxMuonProducer.cc
//
*/
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


class DeDxMuonProducer : public edm::stream::EDProducer<>
{
 public:
  explicit DeDxMuonProducer(const edm::ParameterSet & iConfig):
    token_dedx_(consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("dedx"))),
    token_muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
  {
    produces<edm::ValueMap<reco::DeDxData> >();
  };
  ~DeDxMuonProducer() override {};
  void produce(edm::Event&, const edm::EventSetup&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

 private:
  const edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx_;
  const edm::EDGetTokenT<pat::MuonCollection> token_muons_;
};


void DeDxMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // extract the input collections
  edm::Handle<edm::ValueMap<reco::DeDxData> > dedxHandle;
  iEvent.getByToken(token_dedx_, dedxHandle);
  edm::Handle<pat::MuonCollection> muonHandle;
  iEvent.getByToken(token_muons_, muonHandle);

  // initialise container
  std::vector<reco::DeDxData> dedxVec(muonHandle.isValid() ? muonHandle->size() : 0);

  // loop over muons
  if (dedxHandle.isValid() && muonHandle.isValid()) {
    for (size_t i=0; i<muonHandle->size(); i++) {
      const auto& muon = muonHandle->at(i);
      const auto& track = dynamic_cast<const reco::Muon*>(muon.originalObject())->track();
      if (track.isNull() || !dedxHandle->contains(track.id())) continue;
      dedxVec[i] = (*dedxHandle)[track];
    }
  }

  // store value map
  auto valueMap = std::make_unique<edm::ValueMap<reco::DeDxData> >();
  edm::ValueMap<reco::DeDxData>::Filler filler(*valueMap);
  filler.insert(muonHandle, dedxVec.begin(), dedxVec.end());
  filler.fill();
  iEvent.put(std::move(valueMap));
};


void DeDxMuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("dedx", edm::InputTag("dedxHarmonic2"))->setComment("dedx collection");
  desc.add<edm::InputTag>("muons", edm::InputTag("patMuonsWithoutTrigger"))->setComment("muon collection");
  descriptions.add("dedxMuon", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DeDxMuonProducer);
