/* \class PATGenJetConstituentSelector
 */

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class PATGenJetConstituentSelector : public edm::stream::EDProducer<> {
public:

  using ConstituentsOutput = std::vector<reco::GenParticle>;

  PATGenJetConstituentSelector(edm::ParameterSet const& params) :
    srcToken_{consumes<edm::View<reco::GenJet>>(params.getParameter<edm::InputTag>("src"))},
    selector_{params.getParameter<std::string>("cut")}
  {
    produces<ConstituentsOutput>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions)
  {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("InputTag used for retrieving jets in event.");
    desc.add<std::string>("cut")->setComment("Cut used by which to select jets.  For example:\n"
                                             "  \"pt > 100.0 && abs(rapidity()) < 2.4\".");

    // addDefault must be used here instead of add unless this function is specialized
    // for different sets of template parameter types. Each specialization would need
    // a different module label. Otherwise the generated cfi filenames will conflict
    // for the different plugins.
    descriptions.addDefault(desc);
  }

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override
  {
    auto candsOut = std::make_unique<ConstituentsOutput>();

    edm::Handle<edm::View<reco::GenJet>> h_jets;
    iEvent.getByToken(srcToken_, h_jets);

    // Now set the Ptrs with the orphan handles.
    for (auto const& jet : *h_jets) {
      // Check the selection
      if (selector_(jet)) {
        for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
          const auto& cand = jet.daughterPtr(i);
          if (cand.isNonnull()) {
            auto new_cand = reco::GenParticle(cand->charge(), cand->p4(), math::XYZPoint(0,0,0), cand->pdgId(), 1, true);
            candsOut->push_back(new_cand);
          }
        }
      }
    }
    iEvent.put(std::move(candsOut));
  }

private:
  edm::EDGetTokenT<edm::View<reco::GenJet>> const srcToken_;
  StringCutObjectSelector<reco::GenJet> const selector_;
};

DEFINE_FWK_MODULE(PATGenJetConstituentSelector);
