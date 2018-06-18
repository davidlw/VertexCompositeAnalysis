import FWCore.ParameterSet.Config as cms

d0ana = cms.EDAnalyzer('VertexCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(421),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(321),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False)
                              )

d0ana_mc = cms.EDAnalyzer('VertexCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  hasSwap = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(421),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(321),
  deltaR = cms.untracked.double(0.03),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False)
                              )
