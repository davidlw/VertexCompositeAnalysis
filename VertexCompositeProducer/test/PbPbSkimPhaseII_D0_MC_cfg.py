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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190109/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190109/190109_203806/0000/Pythia8_TuneCP5_5TeV_D0_PiK_prompt_RECO_15.root'
'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190110/Hydjet_5p02TeV_TuneCP5_MTD/Hydjet_5p02TeV_TuneCP5_MTD_RECO_20190110/190112_184141/0000/Hydjet_RECO_99.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# enable TrigReport, TimeReport and MultiThreading
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
#    numberOfThreads = cms.untracked.uint32( 4 ),
#    numberOfStreams = cms.untracked.uint32( 4 ),
#    sizeOfStackForThreadsInKB = cms.untracked.uint32( 10*1024 )
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_HIMinimumBias*_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.pprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
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
#                                         process.phfCoincFilter2Th4 *
                                         process.pprimaryVertexFilter
#                                         process.clusterCompatibilityFilter
                                         )

process.eventFilter_HM = cms.Sequence(
    process.collisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
#process.generalD0CandidatesNew.tkPtSumCut = cms.double(1.6)
#process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(1.0)
#process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPCut = cms.double(0.7)
#process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
#process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)

process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(isWrongSign = cms.bool(True))

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb.root'),
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
