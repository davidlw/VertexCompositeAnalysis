import FWCore.ParameterSet.Config as cms

dimuana = cms.EDAnalyzer('PATCompositeNtupleProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID_dau = cms.untracked.vint32(13, 13),
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidatese:DiMu"),
  doMuon = cms.untracked.bool(True),
  doMuonFull = cms.untracked.bool(True),

  #Trigger info
  TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerPathNames = cms.untracked.vstring(
      #  Double muon triggers
      'HLT_HIL1DoubleMu0_v', # Dimuons
      # Single muon triggers
      'HLT_HIL3Mu12_v', # Electroweak boson
  ),
  triggerFilterNames = cms.untracked.vstring(
      'hltL1fL1sDoubleMu0L1Filtered0',
      'hltL3fL1sSingleMu7L1f0L2f0L3Filtered12',
  ),

  #Filter info
  FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_primaryVertexFilter',
      'Flag_NoScraping',
      'Flag_pileupVertexFilterCut',
      'Flag_pileupVertexFilterCutGplus',
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("",""),
  centralitySrc = cms.InputTag("pACentrality"),

  isEventPlane = cms.bool(False),
  eventplaneSrc = cms.InputTag(""),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),
  massHistPeak = cms.untracked.double(60.0),
  massHistWidth = cms.untracked.double(60.0),
  massHistBins = cms.untracked.int32(1200),

  pTBins = cms.untracked.vdouble(0.0,0.2,1.8,3.0,4.5,6.0,8.0,10.,20.),
  yBins = cms.untracked.vdouble(-2.4,-1.4,0,1.4,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalDiMuCandidates:MVAValues")
)

dimuana_mc = dimuana.clone(
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(True),
)
