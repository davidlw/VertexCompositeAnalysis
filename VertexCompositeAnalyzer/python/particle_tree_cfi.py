import FWCore.ParameterSet.Config as cms

particleAna = cms.EDAnalyzer('ParticleAnalyzer',

  # reconstructed information
  beamSpot = cms.InputTag("offlineBeamSpot"),
  primaryVertices = cms.InputTag("offlinePrimaryVertices"),
  recoParticles = cms.InputTag("generalParticles"),

  # trigger information
  triggerResults = cms.untracked.InputTag("TriggerResults::HLT"),
  triggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
  triggerInfo = cms.untracked.VPSet([
      #cms.PSet(path = cms.string('HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v'), filter = cms.string('hltL1fL1sL1DoubleMuOpenOSCentrality40100L1Filtered0'), minN = cms.int32(2)),
      #cms.PSet(path = cms.string('HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v'), minN = cms.int32(2)),
      #cms.PSet(path = cms.string('HLT_HIL3Mu3_NHitQ10_v'), filter = cms.string('hltL3fL1sL1SingleMu*OpenL1f0L2f0L3Filtered3NHitQ10')),
      #cms.PSet(path = cms.string('HLT_HIL3Mu12_v'), minN = cms.int32(1)),
      cms.PSet(path = cms.string('HLT_HIMinimumBias_*'))
  ]),

  #Filter info
  eventFilterResults = cms.untracked.InputTag("TriggerResults"),
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
  selectEvents = cms.string("eventFilter_HM_step"),

  # centrality and event plane information
  centralityBin = cms.untracked.InputTag("centralityBin","HFtowers"),
  centrality    = cms.untracked.InputTag("hiCentrality"),
  eventPlane    = cms.untracked.InputTag("hiEvtPlaneFlat"),

  # luminosity information
  lumiInfo    = cms.untracked.InputTag("lumiInfo", "brilcalc"),
  lumiScalers = cms.untracked.InputTag("scalersRawToDigi"),
  lumiRecord  = cms.untracked.InputTag("onlineMetaDataDigis"),

  # options
  saveTree  = cms.untracked.bool(True),
  addTrack  = cms.untracked.bool(True),
  addSource = cms.untracked.bool(False),
  dauIDs    = cms.untracked.vint32([13]),
)

particleAna_mc = particleAna.clone(
  # generated information
  genParticles = cms.untracked.InputTag("genMuons"),
  genInfo      = cms.untracked.InputTag("generator"),
)
