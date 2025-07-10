import FWCore.ParameterSet.Config as cms

# Import the original CaloMuonMerger configuration
# You may need to adjust this import based on your CMSSW setup
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons

# Configure the RecoToPatMuonAdapter
recoToPatMuonAdapter = cms.EDProducer("RecoToPatMuonAdapter",
    src = cms.InputTag("calomuons"),  # Output from CaloMuonMerger
    preserveUserData = cms.bool(True),
    addQualityFlags = cms.bool(True)
)

# Create a sequence that runs both modules
caloMuonMergerToPatSequence = cms.Sequence(
    calomuons *  # Original CaloMuonMerger (produces reco::Muon)
    recoToPatMuonAdapter  # Adapter (converts to pat::Muon)
)