import FWCore.ParameterSet.Config as cms

miniAODCaloMuonMergerWithUnpacker = cms.EDProducer("MiniAODCaloMuonMergerWithUnpacker",
    # Input collections (using your existing unpackers)
    unpackedMuons = cms.InputTag("muonUnpacker"),  # From your MuonUnpacker
    packedCandidates = cms.InputTag("packedPFCandidates"),
    unpackedTracks = cms.InputTag("unpackedTracksAndVertices"),  # From your TrackAndVertexUnpacker
    
    # Selection criteria
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    minPt = cms.double(2.0),
    
    # Behavior configuration
    recalculateCaloCompatibility = cms.bool(True),  # Recalculate calo compatibility with enhanced algorithm
    addMissingCaloMuons = cms.bool(True)            # Add calo muons that might be missing from unpacked collection
)