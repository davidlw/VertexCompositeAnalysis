import FWCore.ParameterSet.Config as cms

patCaloMuonMerger = cms.EDProducer("PatCaloMuonMerger",
    # Input collections
    muons = cms.InputTag("slimmedMuons"),
    caloMuons = cms.InputTag("calomuons"),
    tracks = cms.InputTag("generalTracks"),
    
    # Selection criteria
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    
    # Output configuration
    embedCaloMuon = cms.bool(True)
)