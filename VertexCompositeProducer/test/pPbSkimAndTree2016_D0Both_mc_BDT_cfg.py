import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/70000/EE4A348A-CC9D-E711-9738-28924A33AFD2.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_mcRun2_pA_v4'

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
                                         process.hfCoincFilter * 
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

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew * process.generalD0CandidatesNewWrongSign )

###############################################################################################

#process.load("RiceHIG.Skim2013.ppanalysisSkimContentFull_cff")
#process.load("RiceHIG.Skim2013.ppanalysisSkimContentSlim_cff")
process.load("RiceHIG.Skim2013.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pPb_HM.root'),
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

process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
process.d0selectorMCWS = process.d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC = process.d0selectorMC.clone()
process.npd0selectorMC.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v2.root')
process.npd0selectorMCWS = process.npd0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC1 = process.d0selectorMC.clone()
process.npd0selectorMC1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_v2.root')
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

process.d0selectorMCNoDA2D = process.d0selectorMC.clone()
process.npd0selectorMCNoDA2D = process.npd0selectorMC.clone()
process.npd0selectorMC1NoDA2D = process.npd0selectorMC1.clone()
process.d0selectorMCNoErrHitDA2D = process.d0selectorMC.clone()
process.npd0selectorMCNoErrHitDA2D = process.npd0selectorMC.clone()
process.npd0selectorMC1NoErrHitDA2D = process.npd0selectorMC1.clone()
process.d0selectorMCNoDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoDLAngle2D_v1.root')
process.npd0selectorMCNoDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoDLAngle2D_v1.root')
process.npd0selectorMC1NoDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_NoDLAngle2D_v1.root')
process.d0selectorMCNoErrHitDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v1.root')
process.npd0selectorMCNoErrHitDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v1.root')
process.npd0selectorMC1NoErrHitDA2D.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v1.root')

process.d0ana_mc_NoDA2D = process.d0ana_mc.clone()
process.npd0ana_mc_NoDA2D = process.npd0ana_mc.clone()
process.npd0ana1_mc_NoDA2D = process.npd0ana1_mc.clone()
process.d0ana_mc_NoErrHitDA2D = process.d0ana_mc.clone()
process.npd0ana_mc_NoErrHitDA2D = process.npd0ana_mc.clone()
process.npd0ana1_mc_NoErrHitDA2D = process.npd0ana1_mc.clone()
process.d0ana_mc_NoDA2D.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNoDA2D:D0")
process.d0ana_mc_NoDA2D.MVACollection = cms.InputTag("d0selectorMCNoDA2D:MVAValuesNewD0")
process.npd0ana_mc_NoDA2D.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCNoDA2D:D0")
process.npd0ana_mc_NoDA2D.MVACollection = cms.InputTag("npd0selectorMCNoDA2D:MVAValuesNewD0")
process.npd0ana1_mc_NoDA2D.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC1NoDA2D:D0")
process.npd0ana1_mc_NoDA2D.MVACollection = cms.InputTag("npd0selectorMC1NoDA2D:MVAValuesNewD0")
process.d0ana_mc_NoErrHitDA2D.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCNoErrHitDA2D:D0")
process.d0ana_mc_NoErrHitDA2D.MVACollection = cms.InputTag("d0selectorMCNoErrHitDA2D:MVAValuesNewD0")
process.npd0ana_mc_NoErrHitDA2D.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCNoErrHitDA2D:D0")
process.npd0ana_mc_NoErrHitDA2D.MVACollection = cms.InputTag("npd0selectorMCNoErrHitDA2D:MVAValuesNewD0")
process.npd0ana1_mc_NoErrHitDA2D.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC1NoErrHitDA2D:D0")
process.npd0ana1_mc_NoErrHitDA2D.MVACollection = cms.InputTag("npd0selectorMC1NoErrHitDA2D:MVAValuesNewD0")

process.d0ana_seq1 = cms.Sequence(process.d0selectorMCNoDA2D * process.d0ana_mc_NoDA2D)
process.npd0ana_seq1 = cms.Sequence(process.npd0selectorMCNoDA2D * process.npd0ana_mc_NoDA2D)
process.npd0ana1_seq1 = cms.Sequence(process.npd0selectorMC1NoDA2D * process.npd0ana1_mc_NoDA2D)
process.d0ana_seq2 = cms.Sequence(process.d0selectorMCNoErrHitDA2D * process.d0ana_mc_NoErrHitDA2D)
process.npd0ana_seq2 = cms.Sequence(process.npd0selectorMCNoErrHitDA2D * process.npd0ana_mc_NoErrHitDA2D)
process.npd0ana1_seq2 = cms.Sequence(process.npd0selectorMC1NoErrHitDA2D * process.npd0ana1_mc_NoErrHitDA2D)

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
process.pp4 = cms.Path(process.eventFilter_HM * process.d0ana_seq2)
process.pp5 = cms.Path(process.eventFilter_HM * process.npd0ana_seq2)
process.pp6 = cms.Path(process.eventFilter_HM * process.npd0ana1_seq2)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.p1,
    process.p2,
    process.p3,
    process.p4,
    process.p5,
    process.p6,
    process.pp1,
    process.pp2,
    process.pp3,
    process.pp4,
    process.pp5,
    process.pp6
#    process.output_HM_step
)
