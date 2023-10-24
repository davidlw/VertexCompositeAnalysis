#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GenParticleSuperChicFixer : public edm::stream::EDProducer<>
{

    public:

      explicit GenParticleSuperChicFixer(const edm::ParameterSet & iConfig):
          genParToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPars")))
      {
        produces<reco::GenParticleCollection>();
      }
      ~GenParticleSuperChicFixer() override {};

      void produce(edm::Event&, const edm::EventSetup&) override;

      static void fillDescriptions(edm::ConfigurationDescriptions&);

    private:

      const edm::EDGetTokenT<reco::GenParticleCollection> genParToken_;

};

void GenParticleSuperChicFixer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const auto& genPars = iEvent.get(genParToken_);
  const auto& ref = iEvent.getRefBeforePut<reco::GenParticleCollection>();
  std::unique_ptr<reco::GenParticleCollection> output(new reco::GenParticleCollection(genPars));

  // Add generated particles
  for(auto& par : *output) {
    const auto daus = par.daughterRefVector();
    par.resetDaughters(ref.id());
    for (const auto& dau : daus) {
      par.addDaughter(reco::GenParticleRef(ref, dau.key()));
    }
    const auto moms = par.motherRefVector();
    par.resetMothers(ref.id());
    for (const auto& mom : moms) {
      par.addMother(reco::GenParticleRef(ref, mom.key()));
    }
  }

  // Set status 21 to 1
  for(auto& par : *output) {
    if (par.status() == 21) {
      par.setStatus(1);
    }
  }

  iEvent.put(std::move(output));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenParticleSuperChicFixer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genPars", edm::InputTag("genParticles::SIM"))->setComment("generated particle input collection");
  descriptions.add("genSuperChicParticles", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleSuperChicFixer);
