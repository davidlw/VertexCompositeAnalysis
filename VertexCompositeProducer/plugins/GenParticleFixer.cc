#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GenParticleFixer : public edm::stream::EDProducer<>
{

    public:

      explicit GenParticleFixer(const edm::ParameterSet & iConfig):
          genParToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPars"))),
          momPdgId_(iConfig.getParameter<int>("momPdgId")),
          dauPdgId_(iConfig.getParameter<int>("dauPdgId"))
      {
        produces<reco::GenParticleCollection>();
      }
      ~GenParticleFixer() override {};

      void produce(edm::Event&, const edm::EventSetup&) override;

      static void fillDescriptions(edm::ConfigurationDescriptions&);

    private:

      edm::EDGetTokenT<reco::GenParticleCollection> genParToken_;
      int momPdgId_, dauPdgId_;

};

void GenParticleFixer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection> genPars;
  iEvent.getByToken(genParToken_, genPars);

  const reco::GenParticleRefProd ref = iEvent.getRefBeforePut<reco::GenParticleCollection>();
  std::unique_ptr<reco::GenParticleCollection> output(new reco::GenParticleCollection(*genPars));

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
  
  // Extract final state particles without mothers
  std::vector<size_t> idxV;
  for (size_t i=0; i<output->size(); i++) {
    const auto& par = output->at(i);
    if (par.status()==1 && par.numberOfMothers()==0 && par.numberOfDaughters()==0 && fabs(par.pdgId())==dauPdgId_ && par.vertex().Rho()==0.) {
      idxV.push_back(i);
    }
  }

  // Create mother of extracted final state particles
  if (idxV.size() == 2) {
    reco::GenParticle mom;
    int charge = 0;
    reco::Candidate::LorentzVector p4(0,0,0,0);
    for (const auto& i : idxV) {
      auto& par = output->at(i);
      charge += par.charge();
      p4 += par.p4();
    }
    const auto& vtx = output->at(idxV[0]).vertex();
    reco::GenParticle cand(charge, p4, vtx, momPdgId_, 2, true);
    cand.setCollisionId(0);
    for (const auto& i : {reco::GenStatusFlags::kIsPrompt, reco::GenStatusFlags::kIsDecayedLeptonHadron,
                          reco::GenStatusFlags::kIsFirstCopy, reco::GenStatusFlags::kIsLastCopy}) {
      cand.statusFlags().flags_[i] = 1;
    }

    // add mother
    output->push_back(cand);
    const auto& momIdx = output->size()-1;
    reco::GenParticleRef momRef(ref, momIdx);
    output->at(momIdx).resetDaughters(ref.id());
    for (const auto& dauIdx : idxV) {
      output->at(momIdx).addDaughter(reco::GenParticleRef(ref, dauIdx));
      output->at(dauIdx).addMother(momRef);
    }
  }
  else if (!idxV.empty()) {
    throw std::logic_error(Form("[ERROR] Found %lu final state particles without mothers", idxV.size()));
  }

  iEvent.put(std::move(output));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenParticleFixer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genPars", edm::InputTag("genParticles::SIM"))->setComment("generated particle input collection");
  desc.add<int>("momPdgId", 0)->setComment("mother pdg ID");
  desc.add<int>("dauPdgId", 13)->setComment("daughter pdg ID");
  descriptions.add("genParticles", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleFixer);
