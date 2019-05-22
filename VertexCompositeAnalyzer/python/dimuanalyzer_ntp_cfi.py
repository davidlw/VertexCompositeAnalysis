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
  VertexCollection = cms.untracked.InputTag("hiSelectedVertex"),
  VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidatese:DiMu"),
  doMuon = cms.untracked.bool(True),
  doMuonFull = cms.untracked.bool(True),

  #Trigger info
  TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerPathNames = cms.untracked.vstring(
      #  Double muon triggers
      'HLT_HIL1DoubleMu0_v', # Dimuons
      'HLT_HIL1DoubleMu0_part', # Dimuons
      'HLT_HIL1DoubleMu0_2HF_v', # Dimuons
      'HLT_HIL1DoubleMu0_2HF0_v', # Dimuons
      'HLT_HIL1DoubleMu0_2HF_Cent30100_v', # Peripheral dimuons
      'HLT_HIL1DoubleMu0_2HF0_Cent30100_v', # Peripheral dimuons
      'HLT_HIL1DoubleMu10_v', # Z boson
      'HLT_HIUPCL1DoubleMuOpenNotZDCAND_v', # UPC dimuons
      'HLT_HIUPCL1DoubleMuOpenNotHF2_v', # UPC dimuons
      # Single muon triggers
      'HLT_HIL3Mu15_v', # Electroweak boson
      'HLT_HIUPCL1MuOpenNotZDCAND_v', # UPC muons
      'HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v', # UPC muons
  ),
  triggerFilterNames = cms.untracked.vstring(
      'hltHIDoubleMu0L1Filtered',
      'hltHIDoubleMu0L1Filtered',
      'hltHIDoubleMu0MinBiasL1Filtered',
      'hltHIDoubleMu0HFTower0Filtered',
      'hltHIDoubleMu0MinBiasCent30to100L1Filtered',
      'hltHIDoubleMu0HFTower0Cent30to100L1Filtered',
      'hltHIDoubleMu10L1Filtered',
      'hltL1MuOpenL1Filtered3',
      'hltL1MuOpenNotHF2L1Filtered2',
      'hltHISingleMu15L3Filtered',
      'hltL1MuOpenL1Filtered4',
      'hltL1MuOpenNotHF2L1Filtered0',
  ),

  #Filter info
  FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter3',
      'Flag_primaryVertexFilter',
      'Flag_clusterCompatibilityFilter',
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin", "HFtowers"),
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
  doGenMatching = cms.untracked.bool(True),
  doGenMatchingTOF = cms.untracked.bool(True),
)
