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
process.MessageLogger.cerr.FwkReport.reportEvery = 2

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_98.root',
'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_97.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_84.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_80.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_81.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_78.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_77.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_76.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_75.root',
'/store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126/190126_200650/0000/D0_RECO_73.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_upgrade2018_realistic_HI_v6'

# =============== Import Sequences =====================

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices4D"),
    cut = cms.string("!isFake && abs(z) <= 50 && position.Rho <= 5 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.PAcollisionEventSelection = cms.Sequence(
                                         process.hfCoincFilter * 
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.tkPtSumCut = cms.double(2.1)
process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)

#process.generalD0CandidatesNewWrongSign = process.generalD0Candidates.clone(isWrongSign = cms.bool(True))

#process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew * process.generalD0CandidatesNewWrongSign )

process.load('RecoMTD.TrackExtender.trackExtenderWithMTD_cfi')

process.load('RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi')
process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')

process.d0rereco_step = cms.Path(process.mtdClusters * process.mtdTrackingRecHits * process.trackExtenderWithMTD + process.eventFilter_HM * process.generalD0CandidatesNew )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.mtdanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('/afs/cern.ch/user/y/yousen/public/mtdreserach/CMSSW_10_4_0_mtd5/src/VertexCompositeAnalysis/VertexCompositeProducer/test/prompt_d0.root'),
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
