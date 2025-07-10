import FWCore.ParameterSet.Config as cms

miniAODCaloMuonMerger = cms.EDProducer("MiniAODCaloMuonMerger",
    # Input collections
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    # Optional: if you have unpacked tracks from your TrackAndVertexUnpacker
    # tracks = cms.InputTag("unpackedTracksAndVertices"),
    
    # Selection criteria
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    minPt = cms.double(2.0),
    maxEta = cms.double(2.4),
    
    # Behavior configuration
    addCaloMuonsFromPacked = cms.bool(True),  # Add new calo muons from packed candidates
    enhanceExistingMuons = cms.bool(True),    # Enhance existing muons with calo info
    requireTrackerTrack = cms.bool(False)     # Whether to require track details for calo muons
)