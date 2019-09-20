import FWCore.ParameterSet.Config as cms

eventinfoana = cms.EDAnalyzer('EventInfoTreeProducer',
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),

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
      #'Flag_hfCoincFilter',
      'Flag_primaryVertexFilterPA',
      'Flag_NoScraping',
      'Flag_pileupVertexFilterCut_pp5TeV',
      'Flag_pileupVertexFilterCut_pPb8TeV',
      'Flag_pileupVertexFilterCutGplus_pp5TeV',
      'Flag_pileupVertexFilterCutGplus_pPb8TeV'
  ),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("",""),
  centralitySrc = cms.InputTag("pACentrality"),

  isEventPlane = cms.bool(False),
  eventplaneSrc = cms.InputTag("")
)

eventinfoana_mc = eventinfoana.clone()
