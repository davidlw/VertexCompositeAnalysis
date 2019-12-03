import FWCore.ParameterSet.Config as cms

eventinfoana = cms.EDAnalyzer('EventInfoTreeProducer',
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),

  #Trigger info
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
      # ZeroBias trigger
      'HLT_HIZeroBias_v',
      # Minimum Bias trigger
      'HLT_HIMinimumBias_v',
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
      'Flag_hfNegFilterTh8',
      'Flag_hfPosFilterTh7p3',
      'Flag_hfNegFilterTh7p6'
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),

  isEventPlane = cms.bool(False),
  eventplaneSrc = cms.InputTag("")
)

eventinfoana_mc = eventinfoana.clone()
