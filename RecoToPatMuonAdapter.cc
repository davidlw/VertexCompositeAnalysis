#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/UserData.h"

class RecoToPatMuonAdapter : public edm::stream::EDProducer<> {
public:
  explicit RecoToPatMuonAdapter(const edm::ParameterSet&);
  ~RecoToPatMuonAdapter() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Convert reco::Muon to pat::Muon
  pat::Muon convertRecoToPatMuon(const reco::Muon& recoMuon) const;

  // Input tokens
  edm::EDGetTokenT<reco::MuonCollection> recoMuonToken_;
  
  // Configuration
  bool preserveUserData_;
  bool addQualityFlags_;
};

RecoToPatMuonAdapter::RecoToPatMuonAdapter(const edm::ParameterSet& iConfig)
    : recoMuonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      preserveUserData_(iConfig.getParameter<bool>("preserveUserData")),
      addQualityFlags_(iConfig.getParameter<bool>("addQualityFlags")) {
  
  produces<pat::MuonCollection>();
}

void RecoToPatMuonAdapter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input reco muon collection
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(recoMuonToken_, recoMuons);
  
  // Create output PAT muon collection
  auto patMuons = std::make_unique<pat::MuonCollection>();
  
  // Convert each reco muon to PAT muon
  for (const auto& recoMuon : *recoMuons) {
    pat::Muon patMuon = convertRecoToPatMuon(recoMuon);
    patMuons->push_back(patMuon);
  }
  
  // Store output collection
  iEvent.put(std::move(patMuons));
}

pat::Muon RecoToPatMuonAdapter::convertRecoToPatMuon(const reco::Muon& recoMuon) const {
  // Create PAT muon from reco muon
  pat::Muon patMuon(recoMuon);
  
  if (addQualityFlags_) {
    // Add muon quality information as user data
    patMuon.addUserInt("isGlobalMuon", recoMuon.isGlobalMuon() ? 1 : 0);
    patMuon.addUserInt("isTrackerMuon", recoMuon.isTrackerMuon() ? 1 : 0);
    patMuon.addUserInt("isStandAloneMuon", recoMuon.isStandAloneMuon() ? 1 : 0);
    patMuon.addUserInt("isCaloMuon", recoMuon.isCaloMuon() ? 1 : 0);
    patMuon.addUserInt("isPFMuon", recoMuon.isPFMuon() ? 1 : 0);
    
    // Add track quality information if available
    if (recoMuon.innerTrack().isNonnull()) {
      patMuon.addUserFloat("innerTrackChi2", recoMuon.innerTrack()->normalizedChi2());
      patMuon.addUserInt("innerTrackHits", recoMuon.innerTrack()->hitPattern().numberOfValidTrackerHits());
      patMuon.addUserInt("pixelLayersWithMeasurement", recoMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement());
      patMuon.addUserInt("trackerLayersWithMeasurement", recoMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    }
    
    if (recoMuon.globalTrack().isNonnull()) {
      patMuon.addUserFloat("globalTrackChi2", recoMuon.globalTrack()->normalizedChi2());
      patMuon.addUserInt("globalTrackMuonHits", recoMuon.globalTrack()->hitPattern().numberOfValidMuonHits());
    }
    
    // Add number of matched stations
    patMuon.addUserInt("numberOfMatchedStations", recoMuon.numberOfMatchedStations());
    
    // Add calo compatibility if it's a calo muon
    if (recoMuon.isCaloMuon()) {
      patMuon.addUserFloat("caloCompatibility", recoMuon.caloCompatibility());
    }
  }
  
  return patMuon;
}

void RecoToPatMuonAdapter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("muons"));
  desc.add<bool>("preserveUserData", true);
  desc.add<bool>("addQualityFlags", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoToPatMuonAdapter);