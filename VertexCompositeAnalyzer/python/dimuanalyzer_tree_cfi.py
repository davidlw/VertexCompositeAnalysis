import FWCore.ParameterSet.Config as cms

dimuana = cms.EDAnalyzer('PATCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID_dau = cms.untracked.vint32(13, 13),
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVerticesRecovery"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidates:DiMu"),
  GenParticleCollection = cms.untracked.InputTag("genMuons"),
  doMuon = cms.untracked.bool(True),
  doMuonFull = cms.untracked.bool(True),

  #Trigger info
  TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerPathNames = cms.untracked.vstring(
      'HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v', # Peripheral OS dimuons
      'HLT_HIL1DoubleMuOpen_Centrality_50_100_v', # Peripheral dimuons
      'HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v', # Bottomonia
      'HLT_HIL1DoubleMu10_v', # Z bosons
      'HLT_HIUPC_DoubleMu0_NotMBHF2AND_v', # UPC dimuons
      # Single muon triggers
      'HLT_HIL1MuOpen_Centrality_80_100_v', # Peripheral muons
      'HLT_HIL3Mu12_v', # Electroweak bosons
      'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v', # UPC muons
      'HLT_HIL3Mu3_NHitQ10_v1', # Low pT muons
  ),
  triggerFilterNames = cms.untracked.vstring(
      'hltL1fL1sL1DoubleMuOpenOSCentrality40100L1Filtered0',
      'hltL1fL1sL1DoubleMuOpenCentrality50100L1Filtered0',
      'hltL3f0L3Mu2p5NHitQ10L2Mu2FilteredM7toinf',
      'hltL1fL1sL1DoubleMu10L1Filtered0',
      'hltL1sDoubleMu0NotMBHF2AND',
      'hltL1fL1sL1MuOpenCentrality80100L1Filtered0',
      'hltL3fL1sL1SingleMu*OpenL1f*L2f0L3Filtered12',
      'hltL1sSingleMuOpenNotMBHF2AND',
      'hltL3fL1sL1SingleMu*OpenL1f0L2f0L3Filtered3NHitQ10',
  ),

  #Filter info
  FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
      'Flag_clusterCompatibilityFilter',
      'Flag_hfPosFilterTh3', 
      'Flag_hfNegFilterTh3',
      'Flag_hfPosFilterTh4', 
      'Flag_hfNegFilterTh4',
      'Flag_hfPosFilterTh5', 
      'Flag_hfNegFilterTh5',
      'Flag_hfPosFilterTh6', 
      'Flag_hfNegFilterTh6',
      'Flag_hfPosFilterTh7', 
      'Flag_hfNegFilterTh7',
      'Flag_hfPosFilterTh8', 
      'Flag_hfNegFilterTh8'
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),

  isEventPlane = cms.bool(True),
  eventplaneSrc = cms.InputTag("hiEvtPlaneFlat"),

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
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(True),
)
