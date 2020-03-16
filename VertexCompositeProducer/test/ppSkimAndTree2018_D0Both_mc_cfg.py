import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/mc/RunIILowPUAutumn18DR/PrmtD0_pT-1p2_y-2p4_pp_13TeV_Pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/110000/0FCBC797-D975-F442-BEC3-394AC8D6F1F3.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')


# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
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
#                                         process.hfCoincFilter * 
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
#    process.hltHM *
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNewWrongSign = process.generalD0Candidates.clone(isWrongSign = cms.bool(True))

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew )
process.d0rereco_wrongsign_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNewWrongSign )


###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pp_HM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('d0ana_mc_tree.root')
                                   )

process.d0ana_mc.useAnyMVA = cms.bool(True)
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")
process.d0ana_mc_wrongsign = process.d0ana_mc.clone()
process.d0ana_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCWS:D0")
process.d0ana_mc_wrongsign.MVACollection = cms.InputTag("d0selectorMCWS:MVAValuesNewD0")

process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v3.root')
process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
process.d0selectorMC.GBRForestLabel = cms.string('D0Inpp')
process.d0selectorMCWS = process.d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC = process.d0selectorMC.clone()
process.npd0selectorMC.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v3.root')
process.npd0selectorMCWS = process.npd0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC1 = process.d0selectorMC.clone()
process.npd0selectorMC1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_v3.root')
process.npd0selectorMCWS1 = process.npd0selectorMC1.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0ana_mc = process.d0ana_mc.clone()
process.npd0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC:D0")
process.npd0ana_mc.MVACollection = cms.InputTag("npd0selectorMC:MVAValuesNewD0")
process.npd0ana_mc_wrongsign = process.d0ana_mc_wrongsign.clone()
process.npd0ana_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCWS:D0")
process.npd0ana_mc_wrongsign.MVACollection = cms.InputTag("npd0selectorMCWS:MVAValuesNewD0")

process.npd0ana1_mc = process.d0ana_mc.clone()
process.npd0ana1_mc.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC1:D0")
process.npd0ana1_mc.MVACollection = cms.InputTag("npd0selectorMC1:MVAValuesNewD0")
process.npd0ana1_mc_wrongsign = process.d0ana_mc_wrongsign.clone()
process.npd0ana1_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCWS1:D0")
process.npd0ana1_mc_wrongsign.MVACollection = cms.InputTag("npd0selectorMCWS1:MVAValuesNewD0")

process.d0selectorMCNew = process.d0selectorMC.clone()
process.npd0selectorMCNew = process.npd0selectorMC.clone()
process.npd0selectorMC1New = process.npd0selectorMC1.clone()
process.d0selectorMCNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selectorMCNew.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selectorMC1New.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')

process.d0ana_mc_new = process.d0ana_mc.clone()
process.npd0ana_mc_new = process.npd0ana_mc.clone()
process.npd0ana1_mc_new = process.npd0ana1_mc.clone()
process.d0ana_mc_new.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNew:D0")
process.d0ana_mc_new.MVACollection = cms.InputTag("d0selectorMCNew:MVAValuesNewD0")
process.npd0ana_mc_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCNew:D0")
process.npd0ana_mc_new.MVACollection = cms.InputTag("npd0selectorMCNew:MVAValuesNewD0")
process.npd0ana1_mc_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC1New:D0")
process.npd0ana1_mc_new.MVACollection = cms.InputTag("npd0selectorMC1New:MVAValuesNewD0")

process.d0ana_seq1 = cms.Sequence(process.d0selectorMCNew * process.d0ana_mc_new)
process.npd0ana_seq1 = cms.Sequence(process.npd0selectorMCNew * process.npd0ana_mc_new)
process.npd0ana1_seq1 = cms.Sequence(process.npd0selectorMC1New * process.npd0ana1_mc_new)
process.d0ana_seq = cms.Sequence(process.d0selectorMC * process.d0ana_mc)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorMCWS * process.d0ana_mc_wrongsign)
process.npd0ana_seq = cms.Sequence(process.npd0selectorMC * process.npd0ana_mc)
process.npd0ana_wrongsign_seq = cms.Sequence(process.npd0selectorMCWS * process.npd0ana_mc_wrongsign)
process.npd0ana1_seq = cms.Sequence(process.npd0selectorMC1 * process.npd0ana1_mc)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.npd0selectorMCWS1 * process.npd0ana1_mc_wrongsign)

process.p1 = cms.Path(process.eventFilter_HM * process.d0ana_seq)
process.p2 = cms.Path(process.eventFilter_HM * process.d0ana_wrongsign_seq)
process.p3 = cms.Path(process.eventFilter_HM * process.npd0ana_seq)
process.p4 = cms.Path(process.eventFilter_HM * process.npd0ana_wrongsign_seq)
process.p5 = cms.Path(process.eventFilter_HM * process.npd0ana1_seq)
process.p6 = cms.Path(process.eventFilter_HM * process.npd0ana1_wrongsign_seq)

process.pp1 = cms.Path(process.eventFilter_HM * process.d0ana_seq1)
process.pp2 = cms.Path(process.eventFilter_HM * process.npd0ana_seq1)
process.pp3 = cms.Path(process.eventFilter_HM * process.npd0ana1_seq1)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
#    process.d0rereco_wrongsign_step,
#    process.output_HM_step
    process.p1,
    process.p3,
    process.p5,
    process.pp1,
    process.pp2,
    process.pp3
)
