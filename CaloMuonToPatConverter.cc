#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/CaloMuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>

class CaloMuonToPatConverter : public edm::stream::EDProducer<> {
public:
  explicit CaloMuonToPatConverter(const edm::ParameterSet&);
  ~CaloMuonToPatConverter() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  
  // Convert reco::CaloMuon to pat::Muon
  pat::Muon convertCaloMuonToPat(const reco::CaloMuon& caloMuon) const;

  // Input tokens
  edm::EDGetTokenT<reco::CaloMuonCollection> caloMuonToken_;
  
  // Configuration
  double minCaloCompatibility_;
  bool addCaloInfo_;
};

CaloMuonToPatConverter::CaloMuonToPatConverter(const edm::ParameterSet& iConfig)
    : caloMuonToken_(consumes<reco::CaloMuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      minCaloCompatibility_(iConfig.getParameter<double>("minCaloCompatibility")),
      addCaloInfo_(iConfig.getParameter<bool>("addCaloInfo")) {
  
  produces<pat::MuonCollection>();
}

void CaloMuonToPatConverter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input calo muon collection
  edm::Handle<reco::CaloMuonCollection> caloMuons;
  iEvent.getByToken(caloMuonToken_, caloMuons);
  
  // Create output PAT muon collection
  auto patMuons = std::make_unique<pat::MuonCollection>();
  
  // Convert each calo muon to PAT muon
  for (const auto& caloMuon : *caloMuons) {
    // Apply calo compatibility cut
    if (caloMuon.caloCompatibility() < minCaloCompatibility_)
      continue;
      
    pat::Muon patMuon = convertCaloMuonToPat(caloMuon);
    patMuons->push_back(patMuon);
  }
  
  // Store output collection
  iEvent.put(std::move(patMuons));
}

pat::Muon CaloMuonToPatConverter::convertCaloMuonToPat(const reco::CaloMuon& caloMuon) const {
  // Create PAT muon
  pat::Muon patMuon;
  
  // Set basic kinematics from the associated track
  if (caloMuon.trackRef().isNonnull()) {
    const reco::Track& track = *(caloMuon.trackRef());
    
    // Set 4-momentum (assuming muon mass)
    patMuon.setP4(reco::Particle::PolarLorentzVector(
      track.pt(), track.eta(), track.phi(), 0.1057)); // muon mass in GeV
    
    // Set charge and vertex
    patMuon.setCharge(track.charge());
    patMuon.setVertex(reco::Particle::Point(track.vx(), track.vy(), track.vz()));
    
    // Embed the track reference
    patMuon.setInnerTrack(caloMuon.trackRef());
  }
  
  // Set muon type flags
  patMuon.setIsCaloMuon(true);
  
  if (addCaloInfo_) {
    // Add calo muon specific information as user data
    patMuon.addUserFloat("caloCompatibility", caloMuon.caloCompatibility());
    
    // Add calorimeter energy deposits
    patMuon.addUserFloat("caloEnergyEm", caloMuon.caloEnergyEm());
    patMuon.addUserFloat("caloEnergyHad", caloMuon.caloEnergyHad());
    patMuon.addUserFloat("caloEnergyHo", caloMuon.caloEnergyHo());
    
    // Add calorimeter energy deposits in crossed towers
    patMuon.addUserFloat("caloEnergyEmTowers", caloMuon.caloEnergyEmTowers());
    patMuon.addUserFloat("caloEnergyHadTowers", caloMuon.caloEnergyHadTowers());
    patMuon.addUserFloat("caloEnergyHoTowers", caloMuon.caloEnergyHoTowers());
    
    // Store the original calo muon as user data for full access
    patMuon.addUserData("originalCaloMuon", caloMuon);
  }
  
  return patMuon;
}

void CaloMuonToPatConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("calomuons"));
  desc.add<double>("minCaloCompatibility", 0.6);
  desc.add<bool>("addCaloInfo", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CaloMuonToPatConverter);