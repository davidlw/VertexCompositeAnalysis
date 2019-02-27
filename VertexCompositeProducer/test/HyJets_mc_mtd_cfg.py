import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('MTDAnalysis',eras.Phase2C4_timing_layer_bar)

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(1)
)

process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.MessageLogger.cerr.FwkReport.reportEvery = 2

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/Hydjet_5p02TeV_TuneCP5_MTD/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/190127_085926/0000/Hydjet_RECO_99.root'
#" /store/user/anstahll/MTD/MC/NonEmbedded/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/Hydjet_5p02TeV_TuneCP5_MTD/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/190127_085926/0000/Hydjet_RECO_96.root"
#'/store/user/anstahll/MTD/MC/NonEmbedded/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/Hydjet_5p02TeV_TuneCP5_MTD/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/190127_085926/0000/Hydjet_RECO_95.root'
'/store/user/anstahll/MTD/MC/NonEmbedded/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/Hydjet_5p02TeV_TuneCP5_MTD/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190127/190127_085926/0000/Hydjet_RECO_94.root'
),
   skipEvents=cms.untracked.uint32(195)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_upgrade2023_realistic_v2'

# =============== Import Sequences =====================

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices4D"),
    cut = cms.string("!isFake && abs(z) <= 50 && position.Rho <= 5 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)


process.PAcollisionEventSelection = cms.Sequence(
                                         process.hfCoincFilter *
                                         process.PAprimaryVertexFilter
                                         )

process.eventFilter_HM = cms.Sequence(
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )


########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
#process.generalD0CandidatesNew.tkPtSumCut = cms.double(2.1)
#process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
#process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
#process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)


# centrality setup
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetTuneCP5MTD_v1040mtd4x1_mc"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = True
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = False
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices4D")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.hiCentrality.srcEBhits = cms.InputTag("HGCalRecHit","HGCHEBRecHits")
process.hiCentrality.srcEEhits = cms.InputTag("HGCalRecHit","HGCEERecHits")

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

process.d0rereco_step = cms.Path(process.cent_seq + process.eventFilter_HM * process.generalD0CandidatesNew )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.mtdanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string(
#        'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets.root'
         'hyjets.root'
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.output_HM_step
)
