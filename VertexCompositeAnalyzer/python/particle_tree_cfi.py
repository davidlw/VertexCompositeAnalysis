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
      #  Double muon triggers
      cms.PSet(path = cms.string('HLT_PAL1DoubleMuOpen_v'), minN = cms.int32(2)), # Dimuons
      # Single muon triggers
      cms.PSet(path = cms.string('HLT_PAL3Mu12_v'), minN = cms.int32(1)), # Electroweak boson
      # Other triggers
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity120_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity150_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity185_part')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAFullTracks_Multiplicity250_v')), # High multiplicity
      cms.PSet(path = cms.string('HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part')), # Minimum bias
  ]),

  #Filter info
  eventFilterResults = cms.untracked.InputTag("TriggerResults"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter',
      'Flag_primaryVertexFilterPA',
      'Flag_NoScraping',
      'Flag_pileupVertexFilterCut',
      'Flag_pileupVertexFilterCutGplus'
  ),
  selectEvents = cms.string("eventFilter_HM_step"),

  # centrality and event plane information
  centralityBin = cms.untracked.InputTag("",""),
  centrality    = cms.untracked.InputTag("pACentrality"),
  eventPlane    = cms.untracked.InputTag(""),

  # luminosity information
  lumiInfo    = cms.untracked.InputTag("lumiInfo", "brilcalc"),
  lumiScalers = cms.untracked.InputTag("scalersRawToDigi"),
  lumiRecord  = cms.untracked.InputTag("onlineMetaDataDigis"),

  # options
  saveTree  = cms.untracked.bool(True),
  addTrack  = cms.untracked.bool(True),
  addSource = cms.untracked.bool(True),
  dauIDs    = cms.untracked.vint32([13]),
)

particleAna_mc = particleAna.clone(
  # generated information
  genParticles = cms.untracked.InputTag("genMuons"),
  genInfo      = cms.untracked.InputTag("generator"),
)
