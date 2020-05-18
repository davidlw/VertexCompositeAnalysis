import FWCore.ParameterSet.Config as cms

generalana = cms.EDAnalyzer('PATGenericParticleTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID_dau = cms.untracked.vint32(211, 321),
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  VertexCompositeCollection = cms.untracked.InputTag("generalParticleCollection"),
  GenParticleCollection = cms.untracked.InputTag("genMuons"),
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),

  #Trigger info
  TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerPathNames = cms.untracked.vstring(
      #  Double muon triggers
      'HLT_PAL1DoubleMuOpen_v', # Dimuons
      # Single muon triggers
      'HLT_PAL3Mu12_v', # Electroweak boson
      # Other triggers
      'HLT_PAFullTracks_Multiplicity120_v', # High multiplicity
      'HLT_PAFullTracks_Multiplicity150_v', # High multiplicity
      'HLT_PAFullTracks_Multiplicity185_part', # High multiplicity
      'HLT_PAFullTracks_Multiplicity250_v', # High multiplicity
      'HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part', # Minimum bias
  ),
  triggerFilterNames = cms.untracked.vstring(
      'hltL1fL1sDoubleMuOpenBptxANDL1Filtered0',
      'hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12',
  ),

  #Filter info
  FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
      'Flag_clusterCompatibilityFilter'
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
  #MVACollection = cms.InputTag("generalDiMuCandidates:MVAValues")
)

generalana_mc = generalana.clone(
  doGenNtuple = cms.untracked.bool(True),
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(True),
)
