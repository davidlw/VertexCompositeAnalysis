#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>

class PatCaloMuonMerger : public edm::stream::EDProducer<> {
public:
  explicit PatCaloMuonMerger(const edm::ParameterSet&);
  ~PatCaloMuonMerger() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Input tokens
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::CaloMuonCollection> caloMuonToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  
  // Configuration parameters
  double minCaloCompatibility_;
  double deltaR_;
  bool embedCaloMuon_;
};

PatCaloMuonMerger::PatCaloMuonMerger(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      caloMuonToken_(consumes<reco::CaloMuonCollection>(iConfig.getParameter<edm::InputTag>("caloMuons"))),
      trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      minCaloCompatibility_(iConfig.getParameter<double>("minCaloCompatibility")),
      deltaR_(iConfig.getParameter<double>("deltaR")),
      embedCaloMuon_(iConfig.getParameter<bool>("embedCaloMuon")) {
  
  produces<pat::MuonCollection>();
}

PatCaloMuonMerger::~PatCaloMuonMerger() = default;

void PatCaloMuonMerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input collections
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  
  edm::Handle<reco::CaloMuonCollection> caloMuons;
  iEvent.getByToken(caloMuonToken_, caloMuons);
  
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  
  // Create output collection
  auto outputMuons = std::make_unique<pat::MuonCollection>();
  
  // Copy all input PAT muons to output
  for (const auto& muon : *muons) {
    outputMuons->push_back(muon);
  }
  
  // Process calo muons and merge with PAT muons
  for (const auto& caloMuon : *caloMuons) {
    // Check calo compatibility
    if (caloMuon.caloCompatibility() < minCaloCompatibility_)
      continue;
    
    // Check if this calo muon is already represented in PAT muons
    bool alreadyExists = false;
    size_t bestMatchIndex = 0;
    double bestDeltaR = 999.0;
    
    for (size_t i = 0; i < outputMuons->size(); ++i) {
      const auto& patMuon = outputMuons->at(i);
      
      // Calculate deltaR between calo muon and PAT muon
      double dEta = caloMuon.eta() - patMuon.eta();
      double dPhi = reco::deltaPhi(caloMuon.phi(), patMuon.phi());
      double dR = sqrt(dEta*dEta + dPhi*dPhi);
      
      if (dR < deltaR_ && dR < bestDeltaR) {
        bestDeltaR = dR;
        bestMatchIndex = i;
        alreadyExists = true;
      }
    }
    
    if (!alreadyExists) {
      // Create new PAT muon from calo muon
      pat::Muon newPatMuon;
      
      // Set basic kinematics from calo muon track
      if (caloMuon.trackRef().isNonnull()) {
        const reco::Track& track = *(caloMuon.trackRef());
        newPatMuon.setP4(reco::Particle::PolarLorentzVector(
          track.pt(), track.eta(), track.phi(), 0.1057)); // muon mass
        newPatMuon.setCharge(track.charge());
        newPatMuon.setVertex(reco::Particle::Point(track.vx(), track.vy(), track.vz()));
      }
      
      // Set muon type flags
      newPatMuon.setIsCaloMuon(true);
      
      // Embed calo muon information if requested
      if (embedCaloMuon_) {
        // Add calo muon as user data
        newPatMuon.addUserData("caloMuon", caloMuon);
        newPatMuon.addUserFloat("caloCompatibility", caloMuon.caloCompatibility());
      }
      
      outputMuons->push_back(newPatMuon);
    } else {
      // Enhance existing PAT muon with calo information
      if (embedCaloMuon_) {
        pat::Muon& existingMuon = outputMuons->at(bestMatchIndex);
        existingMuon.addUserData("caloMuon", caloMuon);
        existingMuon.addUserFloat("caloCompatibility", caloMuon.caloCompatibility());
        existingMuon.setIsCaloMuon(true);
      }
    }
  }
  
  // Store output collection
  iEvent.put(std::move(outputMuons));
}

void PatCaloMuonMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("caloMuons", edm::InputTag("calomuons"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<double>("minCaloCompatibility", 0.6);
  desc.add<double>("deltaR", 0.1);
  desc.add<bool>("embedCaloMuon", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PatCaloMuonMerger);