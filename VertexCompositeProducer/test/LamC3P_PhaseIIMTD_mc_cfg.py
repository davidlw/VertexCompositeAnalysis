import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANASKIM',eras.Phase2C4_timing_layer_bar)

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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/mc/PhaseIIMTDTDRAutumn18DR/LambdaC_PiKP_prompt_pt1_y4_5p5TeV_TuneCP5_Pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/80000/03FF669F-B381-8240-93B7-FFE5F0C65A00.root')
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
                                         process.PAprimaryVertexFilter
                                         )

process.eventFilter_HM = cms.Sequence( 
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalLamC3PCandidates_cff")
process.generalLamC3PCandidatesNew = process.generalLamC3PCandidates.clone()
process.generalLamC3PCandidatesNew.tkNhitsCut = cms.int32(11)
process.generalLamC3PCandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalLamC3PCandidatesNew.tkDCACut = cms.double(0.5)

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
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices4D") 
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

# MTD RE-RECO
process.reconstruction_step = cms.Path()
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.pfPileUpIso.PFCandidates = cms.InputTag("particleFlowPtrs")
process.pfNoPileUpIso.bottomCollection = cms.InputTag("particleFlowPtrs")
process.reconstruction_step += process.mtdClusters
process.reconstruction_step += process.mtdTrackingRecHits
process.trackExtenderWithMTD.UseVertex = cms.bool(True) #run trackExtender using vertex constrain
process.trackExtenderWithMTD.DZCut = 0.3
process.reconstruction_step += process.trackExtenderWithMTD
process.tofPID.vtxsSrc = cms.InputTag('offlinePrimaryVertices4D')
process.tofPID.fixedT0Error = cms.double(0.035) #put a constant 0.035 [ns] error for each track (cannot
process.reconstruction_step += process.tofPID
                                     
process.generalLamC3PCandidatesNew.trackBeta = cms.InputTag("trackExtenderWithMTD:generalTrackBeta:ANASKIM")
process.generalLamC3PCandidatesNew.trackt0 = cms.InputTag("tofPID:t0:ANASKIM")
process.generalLamC3PCandidatesNew.trackSigmat0 = cms.InputTag("tofPID:sigmat0:ANASKIM")
process.generalLamC3PCandidatesNew.tracktmtd = cms.InputTag("trackExtenderWithMTD:generalTracktmtd:ANASKIM")
process.generalLamC3PCandidatesNew.trackSigmatmtd = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd:ANASKIM")
process.generalLamC3PCandidatesNew.trackp = cms.InputTag("trackExtenderWithMTD:generalTrackp:ANASKIM")
process.generalLamC3PCandidatesNew.trackPathLength = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength:ANASKIM")
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
    process.reconstruction_step,
    process.cent_step,
    process.lamc3prereco_step,
    process.output_HM_step
)
