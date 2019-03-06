import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('MTDAnalysis',eras.Phase2C4_timing_layer_bar)

process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/mc/PhaseIIMTDTDRAutumn18DR/LambdaC_PiKP_prompt_pt1_y4_5p5TeV_TuneCP5_Pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/80000/03FF669F-B381-8240-93B7-FFE5F0C65A00.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_upgrade2023_realistic_v2'

# =============== Import Sequences =====================

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices4D"),
    cut = cms.string("!isFake && abs(z) <= 50 && position.Rho <= 5 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
#process.NoScraping = cms.EDFilter("FilterOutScraping",
#    applyfilter = cms.untracked.bool(True),
#    debugOn = cms.untracked.bool(False),
#    numtrack = cms.untracked.uint32(10),
#    thresh = cms.untracked.double(0.25)
#)

process.PAcollisionEventSelection = cms.Sequence(
#                                         process.hfCoincFilter * 
                                         process.PAprimaryVertexFilter #*
#                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalLamC3PCandidates_cff")
process.generalLamC3PCandidatesNew = process.generalLamC3PCandidates.clone()
#process.generalLamC3PCandidatesNew.tkPtSumCut = cms.double(2.1)
#process.generalLamC3PCandidatesNew.tkEtaDiffCut = cms.double(1.0)
process.generalLamC3PCandidatesNew.tkNhitsCut = cms.int32(11)
process.generalLamC3PCandidatesNew.tkPtErrCut = cms.double(0.1)
#process.generalLamC3PCandidatesNew.tkPCut = cms.double(0.7)
#process.generalLamC3PCandidatesNew.tkEtaCut = cms.double(3.0)
#process.generalLamC3PCandidatesNew.alphaCut = cms.double(1.0)
#process.generalLamC3PCandidatesNew.alpha2DCut = cms.double(1.0)

process.generalLamC3PCandidatesNew.tkDCACut = cms.double(0.5)
#process.generalLamC3PCandidatesNew.dPt3CutMin = cms.double(0.9)
#process.generalLamC3PCandidatesNew.dPt3CutMax = cms.double(2.1)
#process.generalLamC3PCandidatesNew.dY3CutMin = cms.double(-3.2)
#process.generalLamC3PCandidatesNew.dY3CutMax = cms.double(3.2)

process.load('RecoMTD.TrackExtender.trackExtenderWithMTD_cfi')
process.load('RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi')
process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(1)
)

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
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices") 
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
process.centralityBin.Centrality = cms.InputTag("hiCentrality") 
process.centralityBin.centralityVariable = cms.string("HFtowers") 
process.centralityBin.nonDefaultGlauberModel = cms.string("") 
process.hiCentrality.srcEBhits = cms.InputTag("HGCalRecHit","HGCHEBRecHits")
process.hiCentrality.srcEEhits = cms.InputTag("HGCalRecHit","HGCEERecHits")

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)
process.cent_step = cms.Path( process.eventFilter_HM * process.cent_seq )

process.lamc3prereco_step = cms.Path( process.eventFilter_HM * process.generalLamC3PCandidatesNew )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.mtdanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('lambdac_skim.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.cent_step,
    process.lamc3prereco_step,
    process.output_HM_step
)
