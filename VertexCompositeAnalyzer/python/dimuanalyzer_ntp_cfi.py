import FWCore.ParameterSet.Config as cms

dimuana = cms.EDAnalyzer('VertexCompositeNtupleProducer',
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(443),
  PID_dau1 = cms.untracked.int32(13),
  PID_dau2 = cms.untracked.int32(13),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidates:DiMu"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("muons"),
  doMuon = cms.untracked.bool(True),

  isCentrality = cms.bool(False),
  centralityBinLabel = cms.InputTag("newCentralityBin","HFtowers"),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),
  massHistPeak = cms.untracked.double(5.0),
  massHistWidth = cms.untracked.double(5.0),
  massHistBins = cms.untracked.int32(1000),

  pTBins = cms.untracked.vdouble(0.0,0.2,1.8,3.0,4.5,6.0,8.0,10.,20.),
  yBins = cms.untracked.vdouble(-2.4,-1.4,0,1.4,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightPFCandidates:MVAValues")

                              )

dimuana_mc = cms.EDAnalyzer('VertexCompositeNtupleProducer',
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  twoLayerDecay = cms.untracked.bool(False),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(443),
  PID_dau1 = cms.untracked.int32(13),
  PID_dau2 = cms.untracked.int32(13),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidates:DiMu"),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  MuonCollection = cms.untracked.InputTag("muons"),
  doMuon = cms.untracked.bool(True),

  isCentrality = cms.bool(False),
  centralityBinLabel = cms.InputTag("newCentralityBin","HFtowers"),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(True),
  massHistPeak = cms.untracked.double(5.0),
  massHistWidth = cms.untracked.double(5.0),
  massHistBins = cms.untracked.int32(1000),

  pTBins = cms.untracked.vdouble(0.0,0.2,1.8,3.0,4.5,6.0,8.0,10.,20.),
  yBins = cms.untracked.vdouble(-2.4,-1.4,0,1.4,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightPFCandidates:MVAValues")
                              )
