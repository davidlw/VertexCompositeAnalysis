import FWCore.ParameterSet.Config as cms

nTracks = cms.EDProducer('NTrackVertexMapper',
  primaryVertices = cms.InputTag('offlinePrimaryVertices'),
  tracks = cms.InputTag('generalTracks')
)

ntrkFilter = cms.EDFilter('MultFilter',
  primaryVertices = cms.InputTag('offlinePrimaryVertices'),
  nTracksVMap = cms.InputTag('nTracks'),
  nMultMin = cms.untracked.int32(0), # >=
  nMultMax = cms.untracked.int32(2147483647), # <
  useCent = cms.untracked.bool(False),
  centrality = cms.untracked.InputTag('hiCentrality'),
  centralityBin = cms.untracked.InputTag('centralityBin', 'HFtowers'),
  centBinMin = cms.untracked.int32(0), # >=
  centBinMax = cms.untracked.int32(2147483647) # <
)
