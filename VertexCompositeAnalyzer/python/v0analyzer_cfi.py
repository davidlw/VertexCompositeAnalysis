import FWCore.ParameterSet.Config as cms

ksana = cms.EDAnalyzer('VertexCompositeNtupleProducer',

  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(310),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(211),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalV0CandidatesNew:Kshort"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False)
                              )

ksana_mc = cms.EDAnalyzer('VertexCompositeNtupleProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(310),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(211),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalV0CandidatesNew:Kshort"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False)
                              )
