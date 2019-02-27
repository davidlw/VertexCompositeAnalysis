import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('MTDAnalysis',eras.Phase2C4_timing_layer_bar)

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

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


process.GlobalTag.globaltag = '103X_upgrade2023_realistic_v2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/mc/PhaseIIMTDTDRAutumn18DR/D0_PiK_prompt_pt0_y4_5p5TeV_TuneCP5_Pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/30000/FFDBA739-F97F-B347-BAB9-91EC1A7F2CE1.root'
                ),
                            )

# mtd information
process.load('RecoMTD.TrackExtender.trackExtenderWithMTD_cfi')
process.load('RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi')
process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')

#vertex filter
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

# D0 reconstruction
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPCut = cms.double(0.7)

#process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew )
process.d0rereco_step = cms.Path( process.generalD0CandidatesNew )

# D0 ntuple production
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('promptd0_mc_mtd.root')
                                   )

process.d0ana_mc.isUseMtd = cms.untracked.bool(True)
process.d0ana_mc.doRecoNtuple = cms.untracked.bool(True)
process.d0ana_mc.doGenNtuple = cms.untracked.bool(True)
process.d0ana_mc.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0")

process.d0ana_mc_step = cms.Path(process.d0ana_mc)
 
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.d0ana_mc_step
)
