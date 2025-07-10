#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/Handle.h"

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

class MiniAODCaloMuonMerger : public edm::stream::EDProducer<> {
public:
  explicit MiniAODCaloMuonMerger(const edm::ParameterSet&);
  ~MiniAODCaloMuonMerger() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  
  // Helper functions
  bool isCaloMuonCandidate(const pat::Muon& muon) const;
  bool isCaloMuonCandidate(const pat::PackedCandidate& candidate) const;
  pat::Muon createCaloMuonFromPackedCandidate(const pat::PackedCandidate& candidate) const;
  pat::Muon enhanceMuonWithCaloInfo(const pat::Muon& muon, const pat::PackedCandidate& caloCandidate) const;
  double calculateCaloCompatibility(const pat::Muon& muon) const;
  double calculateCaloCompatibility(const pat::PackedCandidate& candidate) const;
  
  // Input tokens
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandidateToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_; // optional, for additional track info
  
  // Configuration parameters
  double minCaloCompatibility_;
  double deltaR_;
  double minPt_;
  double maxEta_;
  bool addCaloMuonsFromPacked_;
  bool enhanceExistingMuons_;
  bool requireTrackerTrack_;
};

MiniAODCaloMuonMerger::MiniAODCaloMuonMerger(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      packedCandidateToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
      minCaloCompatibility_(iConfig.getParameter<double>("minCaloCompatibility")),
      deltaR_(iConfig.getParameter<double>("deltaR")),
      minPt_(iConfig.getParameter<double>("minPt")),
      maxEta_(iConfig.getParameter<double>("maxEta")),
      addCaloMuonsFromPacked_(iConfig.getParameter<bool>("addCaloMuonsFromPacked")),
      enhanceExistingMuons_(iConfig.getParameter<bool>("enhanceExistingMuons")),
      requireTrackerTrack_(iConfig.getParameter<bool>("requireTrackerTrack")) {
  
  // Optional track collection
  if (iConfig.exists("tracks")) {
    trackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  }
  
  produces<pat::MuonCollection>();
}

MiniAODCaloMuonMerger::~MiniAODCaloMuonMerger() = default;

void MiniAODCaloMuonMerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get input collections
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  
  edm::Handle<pat::PackedCandidateCollection> packedCandidates;
  iEvent.getByToken(packedCandidateToken_, packedCandidates);
  
  edm::Handle<reco::TrackCollection> tracks;
  if (!trackToken_.isUninitialized()) {
    iEvent.getByToken(trackToken_, tracks);
  }
  
  // Create output collection
  auto outputMuons = std::make_unique<pat::MuonCollection>();
  
  // Copy all input muons to output with potential enhancements
  for (const auto& muon : *muons) {
    pat::Muon enhancedMuon = muon;
    
    if (enhanceExistingMuons_) {
      // Look for matching packed candidates that could provide calo info
      for (const auto& candidate : *packedCandidates) {
        if (std::abs(candidate.pdgId()) == 13 && // is muon
            isCaloMuonCandidate(candidate) &&
            reco::deltaR(muon, candidate) < deltaR_) {
          
          enhancedMuon = enhanceMuonWithCaloInfo(muon, candidate);
          break; // Use first match
        }
      }
    }
    
    outputMuons->push_back(enhancedMuon);
  }
  
  if (addCaloMuonsFromPacked_) {
    // Look for additional calo muons in packed candidates that aren't already in muon collection
    for (const auto& candidate : *packedCandidates) {
      if (std::abs(candidate.pdgId()) != 13) continue; // not a muon
      if (candidate.pt() < minPt_) continue;
      if (std::abs(candidate.eta()) > maxEta_) continue;
      
      if (isCaloMuonCandidate(candidate)) {
        // Check if this candidate is already represented in the muon collection
        bool alreadyExists = false;
        for (const auto& existingMuon : *outputMuons) {
          if (reco::deltaR(existingMuon, candidate) < deltaR_) {
            alreadyExists = true;
            break;
          }
        }
        
        if (!alreadyExists) {
          // Create new PAT muon from this calo muon candidate
          pat::Muon newMuon = createCaloMuonFromPackedCandidate(candidate);
          outputMuons->push_back(newMuon);
        }
      }
    }
  }
  
  // Store output collection
  iEvent.put(std::move(outputMuons));
}

bool MiniAODCaloMuonMerger::isCaloMuonCandidate(const pat::Muon& muon) const {
  // Check if this muon has calo muon characteristics
  if (muon.isCaloMuon()) return true;
  
  // Additional checks for muons that might be calo-based
  double caloComp = calculateCaloCompatibility(muon);
  return caloComp > minCaloCompatibility_;
}

bool MiniAODCaloMuonMerger::isCaloMuonCandidate(const pat::PackedCandidate& candidate) const {
  // Check if this packed candidate could be a calo muon
  if (std::abs(candidate.pdgId()) != 13) return false;
  
  // Check if it has minimal tracking info (calo muons typically have some track info)
  if (requireTrackerTrack_ && !candidate.hasTrackDetails()) return false;
  
  // Calculate calo compatibility
  double caloComp = calculateCaloCompatibility(candidate);
  return caloComp > minCaloCompatibility_;
}

pat::Muon MiniAODCaloMuonMerger::createCaloMuonFromPackedCandidate(const pat::PackedCandidate& candidate) const {
  pat::Muon newMuon;
  
  // Set basic kinematics
  newMuon.setP4(candidate.p4());
  newMuon.setCharge(candidate.charge());
  newMuon.setVertex(candidate.vertex());
  
  // Set muon type flags
  newMuon.setIsCaloMuon(true);
  
  // Add information as user data
  newMuon.addUserFloat("caloCompatibility", calculateCaloCompatibility(candidate));
  newMuon.addUserFloat("originalPt", candidate.pt());
  newMuon.addUserFloat("originalEta", candidate.eta());
  newMuon.addUserFloat("originalPhi", candidate.phi());
  newMuon.addUserInt("fromPackedCandidate", 1);
  
  // Store energy information if available
  if (candidate.caloFraction() > 0) {
    newMuon.addUserFloat("caloFraction", candidate.caloFraction());
  }
  
  // Add track information if available
  if (candidate.hasTrackDetails()) {
    newMuon.addUserFloat("trackChi2", candidate.pseudoTrack().normalizedChi2());
    newMuon.addUserInt("trackHits", candidate.numberOfHits());
    newMuon.addUserInt("trackPixelHits", candidate.numberOfPixelHits());
  }
  
  return newMuon;
}

pat::Muon MiniAODCaloMuonMerger::enhanceMuonWithCaloInfo(const pat::Muon& muon, const pat::PackedCandidate& caloCandidate) const {
  pat::Muon enhancedMuon = muon;
  
  // Add calo-specific information from the packed candidate
  enhancedMuon.addUserFloat("caloCompatibilityFromPacked", calculateCaloCompatibility(caloCandidate));
  enhancedMuon.addUserFloat("matchedCaloFraction", caloCandidate.caloFraction());
  enhancedMuon.addUserFloat("matchedCaloDeltaR", reco::deltaR(muon, caloCandidate));
  
  // Set calo muon flag if not already set
  if (!muon.isCaloMuon() && calculateCaloCompatibility(caloCandidate) > minCaloCompatibility_) {
    enhancedMuon.setIsCaloMuon(true);
  }
  
  return enhancedMuon;
}

double MiniAODCaloMuonMerger::calculateCaloCompatibility(const pat::Muon& muon) const {
  // Use existing calo compatibility if available
  if (muon.isCaloMuon()) {
    return muon.caloCompatibility();
  }
  
  // Calculate based on available calorimeter information
  // This is a simplified calculation - you may want to implement more sophisticated logic
  double compatibility = 0.0;
  
  // Check isolation and energy deposits
  if (muon.isIsolationValid()) {
    double relIso = (muon.pfIsolationR04().sumChargedHadronPt + 
                     muon.pfIsolationR04().sumNeutralHadronEt + 
                     muon.pfIsolationR04().sumPhotonEt) / muon.pt();
    
    // Higher isolation might indicate calo-based reconstruction
    compatibility = std::min(1.0, relIso / 0.5);
  }
  
  return compatibility;
}

double MiniAODCaloMuonMerger::calculateCaloCompatibility(const pat::PackedCandidate& candidate) const {
  // Calculate calo compatibility for packed candidate
  double compatibility = 0.0;
  
  // Use calo fraction as primary indicator
  if (candidate.caloFraction() > 0) {
    compatibility = candidate.caloFraction();
  }
  
  // Additional criteria based on tracking quality
  if (candidate.hasTrackDetails()) {
    const reco::Track& pseudoTrack = candidate.pseudoTrack();
    
    // Lower quality tracks might indicate calo-based reconstruction
    if (pseudoTrack.normalizedChi2() > 5.0) {
      compatibility += 0.2;
    }
    
    if (candidate.numberOfHits() < 8) {
      compatibility += 0.1;
    }
  } else {
    // No track details available - likely calo-based
    compatibility += 0.5;
  }
  
  return std::min(1.0, compatibility);
}

void MiniAODCaloMuonMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  desc.addOptional<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<double>("minCaloCompatibility", 0.6);
  desc.add<double>("deltaR", 0.1);
  desc.add<double>("minPt", 2.0);
  desc.add<double>("maxEta", 2.4);
  desc.add<bool>("addCaloMuonsFromPacked", true);
  desc.add<bool>("enhanceExistingMuons", true);
  desc.add<bool>("requireTrackerTrack", false);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MiniAODCaloMuonMerger);