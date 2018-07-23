import FWCore.ParameterSet.Config as cms

dsana = cms.EDAnalyzer('VertexCompositeNtupleProducer',
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(431),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(321),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDsToKsKCandidatesNew:DSToKsK"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),

  pTBins = cms.untracked.vdouble(0,10000.0),
  yBins = cms.untracked.vdouble(-1.0,1.0),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalDsToKsKCandidatesNew:MVAValues")
                              )

dsana_mc = cms.EDAnalyzer('VertexCompositeNtupleProducer',
  doGenMatching = cms.untracked.bool(True),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(431),
  PID_dau1 = cms.untracked.int32(310),
  PID_dau2 = cms.untracked.int32(321),
  deltaR = cms.untracked.double(0.03),

  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDsToKsKCandidatesNew:DSToKsK"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),

  pTBins = cms.untracked.vdouble(0,10000.0),
  yBins = cms.untracked.vdouble(-1.0,1.0),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalDsToKsKCandidatesNew:MVAValues")
                              )
