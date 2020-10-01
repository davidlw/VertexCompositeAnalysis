import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/a/anstahll/work/Analysis/Production/DiMuonAnalysis/PbPb2018/AODSKIM/CMSSW_10_3_5/src/Analysis/Skim/test/PbPb_SKIM_v2.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_Prompt_fixEcalADCToGeV_v2')

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.centralityBin)

# Add PbPb event plane
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.hiEvtPlaneFlat.caloCentRef = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2018_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
process.evtplane_seq = cms.Sequence(process.hiEvtPlaneFlat)

# Add the VertexComposite producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles
HybridSoftIdReco2018 = "(isGlobalMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0)"
SoftIdReco = "(isTrackerMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && innerTrack.quality(\"highPurity\"))"
muonSelection = cms.string("(1.2 <= abs(eta) && abs(eta) <= 2.4 && pt >= 0.8) && ("+SoftIdReco+" || "+HybridSoftIdReco2018+")")
candidateSelection = cms.string("mass > 2.4 && mass < 4.3 && abs(rapidity) >= 1.4 && pt < 3.0 && userFloat('vertexProb')>1.E-12")

process.generalMuMuMassMin2Candidates = generalParticles.clone(
    muons = cms.InputTag('patMuonsWithTrigger'),
    finalSelection = candidateSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(13), selection = muonSelection),
        cms.PSet(pdgId = cms.int32(13), selection = muonSelection),
    ]),
)

# Add muon event selection
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("patMuonsWithTrigger"), minNumber = cms.uint32(2))
process.goodMuon = cms.EDFilter("PATMuonSelector",
            src = cms.InputTag("patMuonsWithTrigger"),
            cut = process.generalMuMuMassMin2Candidates.daughterInfo[0].selection,
            )
process.twoGoodMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodMuon"), minNumber = cms.uint32(2))
process.goodDimuon = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = process.generalMuMuMassMin2Candidates.preSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodMuon@+ goodMuon@-')
            )
process.oneGoodDimuon = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDimuon"), minNumber = cms.uint32(1))
process.dimuonEvtSel = cms.Sequence(process.twoMuons * process.goodMuon * process.twoGoodMuons * process.goodDimuon * process.oneGoodDimuon)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Double muon triggers
    'HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v*', # Peripheral OS dimuons
    'HLT_HIL1DoubleMuOpen_Centrality_50_100_v*', # Peripheral dimuons
    'HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v*', # Bottomonia
    'HLT_HIL1DoubleMu10_v*', # Z bosons
    'HLT_HIUPC_DoubleMu0_NotMBHF2AND_v*', # UPC dimuons
    # Single muon triggers
    'HLT_HIL1MuOpen_Centrality_80_100_v*', # Peripheral muons
    'HLT_HIL3Mu12_v*', # Electroweak bosons
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*', # UPC muons
    # MinimumBias
    'HLT_HIMinimumBias_*', # MinimumBias
    ]

# Luminosity producer
process.lumiInfo = cms.EDProducer('LumiProducerFromBrilcalc',
                                  lumiFile = cms.string("./lumiData.csv"),
                                  throwIfNotFound = cms.bool(False),
                                  doBunchByBunch = cms.bool(False))
process.lumi_seq = cms.Sequence(process.lumiInfo)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.dimuonEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM * process.lumi_seq )

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq * process.evtplane_seq)
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.generalMuMuMassMin2Candidates)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.dimucontana = particleAna.clone(
  nTracksVMap = cms.untracked.InputTag("nTracks"),
  triggerInfo = cms.untracked.VPSet([
    #  Double muon triggers
    cms.PSet(path = cms.string('HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v'), minN = cms.int32(2)), # Peripheral OS dimuons
    cms.PSet(path = cms.string('HLT_HIL1DoubleMuOpen_Centrality_50_100_v'), minN = cms.int32(2)), # Peripheral dimuons
    cms.PSet(path = cms.string('HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v'), minN = cms.int32(2)), # Bottomonia
    cms.PSet(path = cms.string('HLT_HIL1DoubleMu10_v'), minN = cms.int32(2)), # Z bosons
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleMu0_NotMBHF2AND_v'), minN = cms.int32(2)), # UPC dimuons
    # Single muon triggers
    cms.PSet(path = cms.string('HLT_HIL1MuOpen_Centrality_80_100_v'), minN = cms.int32(1)), # Peripheral muons
    cms.PSet(path = cms.string('HLT_HIL3Mu12_v'), minN = cms.int32(1)), # Electroweak bosons
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v'), minN = cms.int32(1)), # UPC muons
    # Other triggers
    cms.PSet(path = cms.string('HLT_HIMinimumBias_')), # MinimumBias
    cms.PSet(path = cms.string('HLT_HIL3Mu3_NHitQ10_v1'), minN = cms.int32(1)), # Low pT muons
  ]),
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
  eventFilterResults = cms.untracked.InputTag("TriggerResults::AODSKIM"),
  recoParticles = cms.InputTag("generalMuMuMassMin2Candidates"),
  selectEvents = cms.string("eventFilter_HM_step")
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana.root'))
process.p = cms.EndPath(process.dimucontana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.dimurereco_step,
    process.p
)

# Add recovery for offline primary vertex
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
