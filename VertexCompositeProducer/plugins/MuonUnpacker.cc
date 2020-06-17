#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace pat {

  class MuonUnpacker : public edm::stream::EDProducer<> {

    public:

      explicit MuonUnpacker(const edm::ParameterSet & iConfig):
          muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
          triggerResultToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults")))
      {
        produces<pat::MuonCollection>();
      }
      ~MuonUnpacker() override {};

      void produce(edm::Event&, const edm::EventSetup&) override;

      static void fillDescriptions(edm::ConfigurationDescriptions&);

    private:

      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultToken_;

  };

}

void pat::MuonUnpacker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<pat::MuonCollection> muons;
  edm::Handle<edm::TriggerResults> triggerResults;

  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(triggerResultToken_, triggerResults);

  std::unique_ptr<pat::MuonCollection> output(new pat::MuonCollection(*muons));
  
  // Unpack trigger data for MiniAOD
  for (auto& muon : *output) {
    for(auto& o : muon.triggerObjectMatches()) {
      const_cast<pat::TriggerObjectStandAlone*>(&o)->unpackNamesAndLabels(iEvent, *triggerResults);
    }
  }

  iEvent.put(std::move(output));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void pat::MuonUnpacker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"))->setComment("muon input collection");
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults::HLT"))->setComment("trigger result collection");
  descriptions.add("unpackedMuons", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(MuonUnpacker);
