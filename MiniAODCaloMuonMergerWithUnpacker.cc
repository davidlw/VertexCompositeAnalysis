#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>

class MiniAODCaloMuonMergerWithUnpacker : public edm::stream::EDProducer<> {
public:
  explicit MiniAODCaloMuonMergerWithUnpacker(const edm::ParameterSet&);
  ~MiniAODCaloMuonMergerWithUnpacker() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  
  // Helper functions
  pat::Muon enhanceMuonWithCaloInfo(const pat::Muon& muon) const;
  bool isLikelyCaloMuon(const pat::Muon& muon) const;
  double calculateEnhancedCaloCompatibility(const pat::Muon& muon) const;
  
  // Input tokens
  edm::EDGetTokenT<pat::MuonCollection> unpackedMuonToken_; // From your MuonUnpacker
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandidateToken_;
  edm::EDGetTokenT<reco::TrackCollection> unpackedTrackToken_; // From your TrackAndVertexUnpacker
  
  // Configuration parameters
  double minCaloCompatibility_;
  double deltaR_;
  double minPt_;
  bool recalculateCaloCompatibility_;
  bool addMissingCaloMuons_;
};

MiniAODCaloMuonMergerWithUnpacker::MiniAODCaloMuonMergerWithUnpacker(const edm::ParameterSet& iConfig)
    : unpackedMuonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("unpackedMuons"))),
      packedCandidateToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
      unpackedTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("unpackedTracks"))),
      minCaloCompatibility_(iConfig.getParameter<double>("minCaloCompatibility")),
      deltaR_(iConfig.getParameter<double>("deltaR")),
      minPt_(iConfig.getParameter<double>("minPt")),
      recalculateCaloCompatibility_(iConfig.getParameter<bool>("recalculateCaloCompatibility")),
      addMissingCaloMuons_(iConfig.getParameter<bool>("addMissingCaloMuons")) {
  
  produces<pat::MuonCollection>();
}

void MiniAODCaloMuonMergerWithUnpacker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input collections
  edm::Handle<pat::MuonCollection> unpackedMuons;
  iEvent.getByToken(unpackedMuonToken_, unpackedMuons);
  
  edm::Handle<pat::PackedCandidateCollection> packedCandidates;
  iEvent.getByToken(packedCandidateToken_, packedCandidates);
  
  edm::Handle<reco::TrackCollection> unpackedTracks;
  iEvent.getByToken(unpackedTrackToken_, unpackedTracks);
  
  // Create output collection
  auto outputMuons = std::make_unique<pat::MuonCollection>();
  
  // Process each unpacked muon and enhance with calo info
  for (const auto& muon : *unpackedMuons) {
    pat::Muon enhancedMuon = enhanceMuonWithCaloInfo(muon);
    outputMuons->push_back(enhancedMuon);
  }
  
  if (addMissingCaloMuons_) {
    // Look for potential calo muons in packed candidates that might be missing from unpacked collection
    for (const auto& candidate : *packedCandidates) {
      if (std::abs(candidate.pdgId()) != 13) continue;
      if (candidate.pt() < minPt_) continue;
      
      // Check if this candidate is already represented in the unpacked muon collection
      bool alreadyExists = false;
      for (const auto& unpackedMuon : *outputMuons) {
        if (reco::deltaR(unpackedMuon, candidate) < deltaR_) {
          alreadyExists = true;
          break;
        }
      }
      
      if (!alreadyExists && candidate.caloFraction() > minCaloCompatibility_) {
        // Create new PAT muon for this calo candidate
        pat::Muon newMuon;
        newMuon.setP4(candidate.p4());
        newMuon.setCharge(candidate.charge());
        newMuon.setVertex(candidate.vertex());
        newMuon.setIsCaloMuon(true);
        
        // Add calo-specific information
        newMuon.addUserFloat("caloCompatibility", candidate.caloFraction());
        newMuon.addUserFloat("originalCaloFraction", candidate.caloFraction());
        newMuon.addUserInt("recoveredFromPacked", 1);
        
        outputMuons->push_back(newMuon);
      }
    }
  }
  
  // Store output collection
  iEvent.put(std::move(outputMuons));
}

pat::Muon MiniAODCaloMuonMergerWithUnpacker::enhanceMuonWithCaloInfo(const pat::Muon& muon) const {
  pat::Muon enhancedMuon = muon;
  
  // Recalculate or enhance calo compatibility if requested
  if (recalculateCaloCompatibility_) {
    double newCaloComp = calculateEnhancedCaloCompatibility(muon);
    enhancedMuon.addUserFloat("enhancedCaloCompatibility", newCaloComp);
    
    // Update calo muon flag based on enhanced calculation
    if (newCaloComp > minCaloCompatibility_) {
      enhancedMuon.setIsCaloMuon(true);
    }
  }
  
  // Add flags for analysis
  enhancedMuon.addUserInt("isLikelyCaloMuon", isLikelyCaloMuon(muon) ? 1 : 0);
  
  // Preserve existing calo compatibility if available
  if (muon.isCaloMuon()) {
    enhancedMuon.addUserFloat("originalCaloCompatibility", muon.caloCompatibility());
  }
  
  return enhancedMuon;
}

bool MiniAODCaloMuonMergerWithUnpacker::isLikelyCaloMuon(const pat::Muon& muon) const {
  // Check various criteria that might indicate a calo muon
  
  // 1. Already flagged as calo muon
  if (muon.isCaloMuon()) return true;
  
  // 2. Track quality indicators
  if (muon.innerTrack().isNonnull()) {
    const reco::Track& track = *muon.innerTrack();
    
    // Poor track quality might indicate calo-based reconstruction
    if (track.normalizedChi2() > 10.0) return true;
    if (track.hitPattern().numberOfValidTrackerHits() < 6) return true;
    if (track.hitPattern().numberOfValidPixelHits() == 0) return true;
  }
  
  // 3. Isolation characteristics (high isolation might indicate calo muon)
  if (muon.isIsolationValid()) {
    double relIso = (muon.pfIsolationR04().sumChargedHadronPt + 
                     muon.pfIsolationR04().sumNeutralHadronEt + 
                     muon.pfIsolationR04().sumPhotonEt) / muon.pt();
    if (relIso > 0.8) return true; // High isolation
  }
  
  // 4. Only tracker muon with poor matching
  if (muon.isTrackerMuon() && !muon.isGlobalMuon()) {
    if (muon.numberOfMatchedStations() <= 1) return true;
  }
  
  return false;
}

double MiniAODCaloMuonMergerWithUnpacker::calculateEnhancedCaloCompatibility(const pat::Muon& muon) const {
  // Enhanced calculation combining multiple factors
  
  double compatibility = 0.0;
  
  // Start with existing calo compatibility if available
  if (muon.isCaloMuon()) {
    compatibility = muon.caloCompatibility();
  }
  
  // Factor in track quality
  if (muon.innerTrack().isNonnull()) {
    const reco::Track& track = *muon.innerTrack();
    
    // Poor track quality increases calo likelihood
    if (track.normalizedChi2() > 5.0) compatibility += 0.2;
    if (track.hitPattern().numberOfValidTrackerHits() < 8) compatibility += 0.1;
    if (track.hitPattern().numberOfValidPixelHits() == 0) compatibility += 0.2;
  }
  
  // Factor in muon station matching
  if (muon.isTrackerMuon()) {
    if (muon.numberOfMatchedStations() <= 1) compatibility += 0.3;
    if (muon.numberOfMatches() <= 1) compatibility += 0.2;
  }
  
  // Factor in isolation
  if (muon.isIsolationValid()) {
    double relIso = (muon.pfIsolationR04().sumChargedHadronPt + 
                     muon.pfIsolationR04().sumNeutralHadronEt + 
                     muon.pfIsolationR04().sumPhotonEt) / muon.pt();
    
    // High isolation might indicate calo reconstruction
    if (relIso > 0.5) compatibility += 0.2;
  }
  
  return std::min(1.0, compatibility);
}

void MiniAODCaloMuonMergerWithUnpacker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("unpackedMuons", edm::InputTag("muonUnpacker"));  // From your MuonUnpacker
  desc.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("unpackedTracks", edm::InputTag("unpackedTracksAndVertices"));  // From your TrackAndVertexUnpacker
  desc.add<double>("minCaloCompatibility", 0.6);
  desc.add<double>("deltaR", 0.1);
  desc.add<double>("minPt", 2.0);
  desc.add<bool>("recalculateCaloCompatibility", true);
  desc.add<bool>("addMissingCaloMuons", true);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MiniAODCaloMuonMergerWithUnpacker);