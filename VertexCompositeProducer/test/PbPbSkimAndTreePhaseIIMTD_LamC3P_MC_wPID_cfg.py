import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras

process = cms.Process('ANASKIM',eras.Phase2C4_timing_layer_bar)

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.register ('pTMin',
                  0.95, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "Minimum pT of LamC3P")
options.register ('pTMax',
                  10000.0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "Maximum pT of LamC3P")
options.register ('yMin',
                  -3.05, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "Minimum pT of LamC3P")
options.register ('yMax',
                  3.05, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "Maximum pT of LamC3P")


# get and parse the command line arguments
options.parseArguments()

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
process.MessageLogger.cerr.FwkReport.reportEvery = 2

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/mc/PhaseIIMTDTDRAutumn18DR/MinBias_Hydjet_Drume5_5p5TeV_TuneCP5_Pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/30000/C14FE87B-FFB6-A043-BFE4-953DB6D2210C.root')
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
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

########## LamC3P candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalLamC3PCandidates_cff")
process.generalLamC3PCandidatesNew = process.generalLamC3PCandidates.clone()
process.generalLamC3PCandidatesNew.tkNhitsCut = cms.int32(11)
process.generalLamC3PCandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalLamC3PCandidatesNew.tkEtaCut = cms.double(3.0)
process.generalLamC3PCandidatesNew.tkPtMidCut = cms.double(0.8)
process.generalLamC3PCandidatesNew.tkPFwdCut = cms.double(0.7)
process.generalLamC3PCandidatesNew.tkPtFwdCut = cms.double(0.4)
process.generalLamC3PCandidatesNew.VtxChiProbCut = cms.double(0.05)

process.generalLamC3PCandidatesNew.tkDCACut = cms.double(0.5)
process.generalLamC3PCandidatesNew.dPt3CutMin = cms.double(options.pTMin)
process.generalLamC3PCandidatesNew.dPt3CutMax = cms.double(options.pTMax)
process.generalLamC3PCandidatesNew.dY3CutMin = cms.double(options.yMin)
process.generalLamC3PCandidatesNew.dY3CutMax = cms.double(options.yMax)

process.generalLamC3PCandidatesNew.isTOFPID = cms.bool(True)
process.generalLamC3PCandidatesNew.nSigmaTOFPID = cms.double(1.5)

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


process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('hyjets_mc_lamc3p_mtd_tree.root')
                                   )

process.lamc3pana_mc.isUseMtd = cms.untracked.bool(True)
process.lamc3pana_mc.doRecoNtuple = cms.untracked.bool(True)
process.lamc3pana_mc.doGenNtuple = cms.untracked.bool(True)
process.lamc3pana_mc.doGenMatching = cms.untracked.bool(False)
process.lamc3pana_mc.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")
process.lamc3pana_mc.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMC:LamC3P")
process.lamc3pana_mc.MVACollection = cms.InputTag("lamc3pselectorMC:MVAValuesNewLamC3P")
process.lamc3pana_mc.isCentrality = cms.bool(True)
process.lamc3pselectorMC.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")

process.lamc3pana_seq = cms.Sequence(process.lamc3pselectorMC * process.lamc3pana_mc)
process.lamc3pana_step = cms.Path( process.eventFilter_HM * process.lamc3pana_seq )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.mtdanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('hyjets_lambdac.root'),
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
    process.lamc3pana_step,
#    process.output_HM_step
)

