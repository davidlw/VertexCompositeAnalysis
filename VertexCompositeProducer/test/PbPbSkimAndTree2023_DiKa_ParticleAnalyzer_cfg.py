import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Modifier_highBetaStar_2018_cff import highBetaStar_2018
Run3_2023_highBetaStar = cms.ModifierChain(eras.Run3_2023, highBetaStar_2018)
process = cms.Process('ANASKIM', Run3_2023_highBetaStar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/SKIM/HBS/HIForward0/SKIM_AOD_HIGHBETASTAR_HIForward0_HIRun2023A_2023_10_13/231013_081452/0000/reco_RAW2DIGI_L1Reco_RECO_HIGHBETASTAR_73.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

# Set ZDC information
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("emap_2023_newZDC_v3.txt")))
)

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
        label = cms.untracked.string("HFtowers")
        )
    ]
)
process.cent_seq = cms.Sequence(process.centralityBin)

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

# DiKa selection
kaonSelection = cms.string("")#(pt > 0.0 && abs(eta) < 3.0) && quality(\"highPurity\")")
kaonFinalSelection = cms.string("")#abs(userFloat(\"dzSig\"))<3.0 && abs(userFloat(\"dxySig\"))<3.0")
diKaSelection = cms.string("charge==0")
process.diKa = generalParticles.clone(
    pdgId = cms.int32(333),
    preSelection = diKaSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(+1), selection = kaonSelection, finalSelection = kaonFinalSelection),
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(-1), selection = kaonSelection, finalSelection = kaonFinalSelection),
    ]),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)
process.oneDiKa = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diKa"), minNumber = cms.uint32(1))

# Add diKa event selection
process.twoTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(2))
process.hpTracks = cms.EDFilter("TrackSelector", src = cms.InputTag("generalTracks"), cut = cms.string("quality(\"highPurity\")"))
process.hpCands = cms.EDProducer("ChargedCandidateProducer", src = cms.InputTag("hpTracks"), particleType = cms.string('pi+'))
process.maxTwoHPCands = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("hpCands"), minNumber = cms.uint32(0), maxNumber = cms.uint32(2))
process.goodTracks = cms.EDFilter("TrackSelector",
            src = cms.InputTag("generalTracks"),
            cut = kaonSelection,
            )
process.twoGoodTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("goodTracks"), minNumber = cms.uint32(2))
process.goodKaons = cms.EDProducer("ChargedCandidateProducer",
            src = cms.InputTag("goodTracks"),
            particleType = cms.string('pi+')
            )
process.goodDiKaons = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = diKaSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodKaons@+ goodKaons@-')
            )
process.oneGoodDiKa = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDiKaons"), minNumber = cms.uint32(1))
process.diKaEvtSel = cms.Sequence(process.twoTracks * process.hpTracks * process.hpCands * process.maxTwoHPCands * process.goodTracks * process.twoGoodTracks * process.goodKaons * process.goodDiKaons * process.oneGoodDiKa)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC zero bias triggers
    'HLT_HIZeroBias_v*',
    'HLT_HIZeroBias_HighRate_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
    # UPC SDC triggers
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*',
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.diKaEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.diKa_rereco_step = cms.Path(process.eventFilter_HM * process.diKa * process.oneDiKa * process.cent_seq)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.diKaAna = particleAna.clone(
  recoParticles = cms.InputTag("diKa"),
  selectEvents = cms.string("diKa_rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
  ),
  triggerInfo = cms.untracked.VPSet([
    # UPC low pT triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    # UPC zero bias triggers
    cms.PSet(path = cms.string('HLT_HIZeroBias_v*')),
    cms.PSet(path = cms.string('HLT_HIZeroBias_HighRate_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC ZDC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*')),
  ]),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('diKa_ana.root'))
process.p = cms.EndPath(process.diKaAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.diKa_rereco_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter ]

process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter *
    process.hfPosFilterNTh8_seq *
    process.hfNegFilterNTh8_seq *
    process.diKaEvtSel
)

for P in eventFilterPaths:
    process.schedule.insert(0, P)
