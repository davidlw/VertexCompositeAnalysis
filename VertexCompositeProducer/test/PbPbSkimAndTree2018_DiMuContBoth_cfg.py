import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 20

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/hidata/HIRun2018A/HIMinimumBias1/AOD/PromptReco-v1/000/326/577/00000/532A6440-350D-2446-89FB-3CA9E8335E4F.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_dataRun2_Prompt_v3'
### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "103X_dataRun2_Prompt_v3"
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_HIMinimumBias*_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")

process.pprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

process.towersAboveThreshold.minimumE = cms.double(4.0)
process.hfPosFilter2 = process.hfPosFilter.clone(minNumber=cms.uint32(2))
process.hfNegFilter2 = process.hfNegFilter.clone(minNumber=cms.uint32(2))
process.phfCoincFilter2Th4 = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers *
    process.hfPosFilter2 *
    process.hfNegFilter2 )

process.clusterCompatibilityFilter  = cms.EDFilter('HIClusterCompatibilityFilter',
   cluscomSrc = cms.InputTag("hiClusterCompatibility"),
   minZ          = cms.double(-20.0),
   maxZ          = cms.double(20.05),
   clusterPars   = cms.vdouble(0.0,0.0045),
   nhitsTrunc    = cms.int32(150),
   clusterTrunc  = cms.double(2.0)
)

process.collisionEventSelection = cms.Sequence(
                                         process.phfCoincFilter2Th4 * 
                                         process.pprimaryVertexFilter *
                                         process.clusterCompatibilityFilter
                                         )

process.eventFilter_HM = cms.Sequence( 
#    process.hltHM *
    process.offlinePrimaryVerticesRecovery *
    process.collisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalMuMuContinuimCandidatesWrongSign = process.generalMuMuContinuimCandidates.clone(isWrongSign = cms.bool(True))

process.dimurereco_step = cms.Path( process.eventFilter_HM * process.generalMuMuContinuimCandidates )
process.dimurerecowrongsign_step = cms.Path( process.eventFilter_HM * process.generalMuMuContinuimCandidatesWrongSign )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('dimucontanatree.root')
                                   )

process.dimucontana.doMuonFull = cms.untracked.bool(False)
process.dimucontana_wrongsign.doMuonFull = cms.untracked.bool(False)
process.dimucontana.isCentrality = cms.bool(True)
process.dimucontana_wrongsign.isCentrality = cms.bool(True)
process.dimucontana_seq = cms.Sequence(process.dimucontana)
process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)

process.ptree = cms.Path(process.eventFilter_HM * process.cent_seq * process.dimucontana_seq)
process.ptree1 = cms.Path(process.eventFilter_HM * process.cent_seq * process.dimucontana_wrongsign_seq)

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb_DiMuCont.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)
process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.dimurereco_step,
#    process.dimurerecowrongsign_step,
    process.ptree,
#    process.ptree1,
    process.output_HM_step
)

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
