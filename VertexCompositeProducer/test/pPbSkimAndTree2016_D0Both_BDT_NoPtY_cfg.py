import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
#'/store/hidata/PARun2016B/PAMinimumBias1/AOD/PromptReco-v1/000/285/216/00000/04214B66-30AC-E611-BFBC-FA163EBFF447.root'
# '/store/hidata/PARun2016B/PAMinimumBias1/AOD/PromptReco-v1/000/285/090/00000/DC56AAB2-5EAB-E611-A27A-FA163E251515.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(400))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_PAFullTracks_Multiplicity120_v*', # High multiplicity
    'HLT_PAFullTracks_Multiplicity150_v*', # High multiplicity
    'HLT_PAFullTracks_Multiplicity185_part*', # High multiplicity
    'HLT_PAFullTracks_Multiplicity250_v*', # High multiplicity
    'HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*', # Minimum bias
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
#process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)
#remove the default dz1p0 filter
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.colEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )
#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.trkPtSumCut = cms.double(1.6)
process.generalD0CandidatesNew.trkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)
process.generalD0CandidatesNew.dPtCut = cms.double(1.0)

process.generalD0CandidatesNewWrongSign = process.generalD0CandidatesNew.clone(isWrongSign = cms.bool(True))

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew * process.generalD0CandidatesNewWrongSign )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pPb_HM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

# produce D0 trees
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('d0ana_tree.root')
                                   )

# set up selectors
process.d0selector = process.d0selectorBDTPreCut.clone()
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt24MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.d0selector.multMin = cms.untracked.double(0)
process.d0selector.multMax = cms.untracked.double(100000)
process.d0selector.candpTMin = cms.untracked.double(2.0)
process.d0selector.candpTMax = cms.untracked.double(4.0)
process.d0selectorWS = process.d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selector = process.d0selector.clone()
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt24MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWS = process.d0selectorWS.clone()
process.npd0selectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt24MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')

process.npd0selector1 = process.d0selector.clone()
process.npd0selector1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt24MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWS1 = process.d0selectorWS.clone()
process.npd0selectorWS1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt24MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.multMin = cms.untracked.double(0)
process.d0ana.multMax = cms.untracked.double(100000)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana_wrongsign.useAnyMVA = cms.bool(True)
process.d0ana_wrongsign.multMin = cms.untracked.double(0)
process.d0ana_wrongsign.multMax = cms.untracked.double(100000)
process.d0ana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWS:D0")
process.d0ana_wrongsign.MVACollection = cms.InputTag("d0selectorWS:MVAValuesNewD0")

process.npd0ana = process.d0ana.clone()
process.npd0ana.VertexCompositeCollection = cms.untracked.InputTag("npd0selector:D0")
process.npd0ana.MVACollection = cms.InputTag("npd0selector:MVAValuesNewD0")
process.npd0ana_wrongsign = process.d0ana_wrongsign.clone()
process.npd0ana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWS:D0")
process.npd0ana_wrongsign.MVACollection = cms.InputTag("npd0selectorWS:MVAValuesNewD0")

process.npd0ana1 = process.d0ana.clone()
process.npd0ana1.VertexCompositeCollection = cms.untracked.InputTag("npd0selector1:D0")
process.npd0ana1.MVACollection = cms.InputTag("npd0selector1:MVAValuesNewD0")
process.npd0ana1_wrongsign = process.d0ana_wrongsign.clone()
process.npd0ana1_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWS1:D0")
process.npd0ana1_wrongsign.MVACollection = cms.InputTag("npd0selectorWS1:MVAValuesNewD0")

process.d0selectorNew = process.d0selector.clone()
process.npd0selectorNew = process.npd0selector.clone()
process.npd0selector1New = process.npd0selector1.clone()
process.d0selectorNewReduced = process.d0selector.clone()
process.npd0selectorNewReduced = process.npd0selector.clone()
process.npd0selector1NewReduced = process.npd0selector1.clone()
process.d0selectorNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorNew.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selector1New.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.d0selectorNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selector1NewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.d0selectorNew.candpTMin = cms.untracked.double(4.0)
process.d0selectorNew.candpTMax = cms.untracked.double(6.0)
process.npd0selectorNew.candpTMin = cms.untracked.double(4.0)
process.npd0selectorNew.candpTMax = cms.untracked.double(6.0)
process.npd0selector1New.candpTMin = cms.untracked.double(4.0)
process.npd0selector1New.candpTMax = cms.untracked.double(6.0)
process.d0selectorNewReduced.candpTMin = cms.untracked.double(6.0)
process.d0selectorNewReduced.candpTMax = cms.untracked.double(8.0)
process.npd0selectorNewReduced.candpTMin = cms.untracked.double(6.0)
process.npd0selectorNewReduced.candpTMax = cms.untracked.double(8.0)
process.npd0selector1NewReduced.candpTMin = cms.untracked.double(6.0)
process.npd0selector1NewReduced.candpTMax = cms.untracked.double(8.0)

process.d0selectorWSNew = process.d0selectorWS.clone()
process.npd0selectorWSNew = process.npd0selectorWS.clone()
process.npd0selectorWS1New = process.npd0selectorWS1.clone()
process.d0selectorWSNewReduced = process.d0selectorWS.clone()
process.npd0selectorWSNewReduced = process.npd0selectorWS.clone()
process.npd0selectorWS1NewReduced = process.npd0selectorWS1.clone()
process.d0selectorWSNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWSNew.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWS1New.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt46MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.d0selectorWSNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWSNewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWS1NewReduced.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt68MassPeak_NoPtYPtErrNHitDLAngle2D_v3.root')
process.d0selectorWSNew.candpTMin = cms.untracked.double(4.0)
process.d0selectorWSNew.candpTMax = cms.untracked.double(6.0)
process.npd0selectorWSNew.candpTMin = cms.untracked.double(4.0)
process.npd0selectorWSNew.candpTMax = cms.untracked.double(6.0)
process.npd0selectorWS1New.candpTMin = cms.untracked.double(4.0)
process.npd0selectorWS1New.candpTMax = cms.untracked.double(6.0)
process.d0selectorWSNewReduced.candpTMin = cms.untracked.double(6.0)
process.d0selectorWSNewReduced.candpTMax = cms.untracked.double(8.0)
process.npd0selectorWSNewReduced.candpTMin = cms.untracked.double(6.0)
process.npd0selectorWSNewReduced.candpTMax = cms.untracked.double(8.0)
process.npd0selectorWS1NewReduced.candpTMin = cms.untracked.double(6.0)
process.npd0selectorWS1NewReduced.candpTMax = cms.untracked.double(8.0)

process.d0ana_pt24 = process.d0ana.clone()
process.npd0ana_pt24 = process.npd0ana.clone()
process.npd0ana1_pt24 = process.npd0ana1.clone()
process.d0ana_pt46 = process.d0ana.clone()
process.npd0ana_pt46 = process.npd0ana.clone()
process.npd0ana1_pt46 = process.npd0ana1.clone()
process.d0ana_pt68 = process.d0ana.clone()
process.npd0ana_pt68 = process.npd0ana.clone()
process.npd0ana1_pt68 = process.npd0ana1.clone()
process.d0ana_pt46.VertexCompositeCollection = cms.untracked.InputTag("d0selectorNew:D0")
process.d0ana_pt46.MVACollection = cms.InputTag("d0selectorNew:MVAValuesNewD0")
process.npd0ana_pt46.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorNew:D0")
process.npd0ana_pt46.MVACollection = cms.InputTag("npd0selectorNew:MVAValuesNewD0")
process.npd0ana1_pt46.VertexCompositeCollection = cms.untracked.InputTag("npd0selector1New:D0")
process.npd0ana1_pt46.MVACollection = cms.InputTag("npd0selector1New:MVAValuesNewD0")
process.d0ana_pt68.VertexCompositeCollection = cms.untracked.InputTag("d0selectorNewReduced:D0")
process.d0ana_pt68.MVACollection = cms.InputTag("d0selectorNewReduced:MVAValuesNewD0")
process.npd0ana_pt68.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorNewReduced:D0")
process.npd0ana_pt68.MVACollection = cms.InputTag("npd0selectorNewReduced:MVAValuesNewD0")
process.npd0ana1_pt68.VertexCompositeCollection = cms.untracked.InputTag("npd0selector1NewReduced:D0")
process.npd0ana1_pt68.MVACollection = cms.InputTag("npd0selector1NewReduced:MVAValuesNewD0")

process.d0ana_wrongsign_pt24 = process.d0ana_wrongsign.clone()
process.npd0ana_wrongsign_pt24 = process.npd0ana_wrongsign.clone()
process.npd0ana1_wrongsign_pt24 = process.npd0ana1_wrongsign.clone()
process.d0ana_wrongsign_pt46 = process.d0ana_wrongsign.clone()
process.npd0ana_wrongsign_pt46 = process.npd0ana_wrongsign.clone()
process.npd0ana1_wrongsign_pt46 = process.npd0ana1_wrongsign.clone()
process.d0ana_wrongsign_pt68 = process.d0ana_wrongsign.clone()
process.npd0ana_wrongsign_pt68 = process.npd0ana_wrongsign.clone()
process.npd0ana1_wrongsign_pt68 = process.npd0ana1_wrongsign.clone()
process.d0ana_wrongsign_pt46.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWSNew:D0")
process.d0ana_wrongsign_pt46.MVACollection = cms.InputTag("d0selectorWSNew:MVAValuesNewD0")
process.npd0ana_wrongsign_pt46.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWSNew:D0")
process.npd0ana_wrongsign_pt46.MVACollection = cms.InputTag("npd0selectorWSNew:MVAValuesNewD0")
process.npd0ana1_wrongsign_pt46.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWS1New:D0")
process.npd0ana1_wrongsign_pt46.MVACollection = cms.InputTag("npd0selectorWS1New:MVAValuesNewD0")
process.d0ana_wrongsign_pt68.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWSNewReduced:D0")
process.d0ana_wrongsign_pt68.MVACollection = cms.InputTag("d0selectorWSNewReduced:MVAValuesNewD0")
process.npd0ana_wrongsign_pt68.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWSNewReduced:D0")
process.npd0ana_wrongsign_pt68.MVACollection = cms.InputTag("npd0selectorWSNewReduced:MVAValuesNewD0")
process.npd0ana1_wrongsign_pt68.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWS1NewReduced:D0")
process.npd0ana1_wrongsign_pt68.MVACollection = cms.InputTag("npd0selectorWS1NewReduced:MVAValuesNewD0")

process.d0ana_seq = cms.Sequence(process.eventFilter_HM * process.d0selector * process.d0ana_pt24)
process.npd0ana_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector * process.npd0ana_pt24)
process.npd0ana1_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector1 * process.npd0ana1_pt24)
process.d0ana_seq1 = cms.Sequence(process.eventFilter_HM * process.d0selectorNew * process.d0ana_pt46)
process.npd0ana_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorNew * process.npd0ana_pt46)
process.npd0ana1_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selector1New * process.npd0ana1_pt46)
process.d0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorNewReduced * process.d0ana_pt68)
process.npd0ana_seq2 = cms.Sequence(process.eventFilter_HM * process.npd0selectorNewReduced * process.npd0ana_pt68)
process.npd0ana1_seq2 = cms.Sequence(process.eventFilter_HM * process.npd0selector1NewReduced * process.npd0ana1_pt68)

process.d0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.d0selectorWS * process.d0ana_wrongsign_pt24)
process.npd0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS * process.npd0ana_wrongsign_pt24)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1 * process.npd0ana1_wrongsign_pt24)
process.d0ana_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.d0selectorWSNew * process.d0ana_wrongsign_pt46)
process.npd0ana_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWSNew * process.npd0ana_wrongsign_pt46)
process.npd0ana1_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1New * process.npd0ana1_wrongsign_pt46)
process.d0ana_wrongsign_seq2 = cms.Sequence(process.eventFilter_HM * process.d0selectorWSNewReduced * process.d0ana_wrongsign_pt68)
process.npd0ana_wrongsign_seq2 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWSNewReduced * process.npd0ana_wrongsign_pt68)
process.npd0ana1_wrongsign_seq2 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1NewReduced * process.npd0ana1_wrongsign_pt68)

# eventinfoana must be in EndPath, and process.eventinfoana.selectEvents must be the name of eventFilter_HM Path
process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.pevt = cms.EndPath(process.eventinfoana)

process.pa = cms.Path(process.d0ana_seq)
process.pa1 = cms.Path(process.d0ana_wrongsign_seq)
process.pa4 = cms.Path(process.npd0ana_seq)
process.pa5 = cms.Path(process.npd0ana_wrongsign_seq)
process.pb4 = cms.Path(process.npd0ana1_seq)
process.pb5 = cms.Path(process.npd0ana1_wrongsign_seq)

process.pp1 = cms.Path(process.d0ana_seq1)
process.pp2 = cms.Path(process.npd0ana_seq1)
process.pp3 = cms.Path(process.npd0ana1_seq1)
process.pp4 = cms.Path(process.d0ana_seq2)
process.pp5 = cms.Path(process.npd0ana_seq2)
process.pp6 = cms.Path(process.npd0ana1_seq2)

process.ppp1 = cms.Path(process.d0ana_wrongsign_seq1)
process.ppp2 = cms.Path(process.npd0ana_wrongsign_seq1)
process.ppp3 = cms.Path(process.npd0ana1_wrongsign_seq1)
process.ppp4 = cms.Path(process.d0ana_wrongsign_seq2)
process.ppp5 = cms.Path(process.npd0ana_wrongsign_seq2)
process.ppp6 = cms.Path(process.npd0ana1_wrongsign_seq2)

# Add the Conversion tree
process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.track_cff")
process.ptrk = cms.Path(process.eventFilter_HM * process.track_ana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.pa,
    process.pa4,
    process.pb4,
    process.pp1,
    process.pp2,
    process.pp3,
    process.pp4,
    process.pp5,
    process.pp6,
    process.ppp1,
    process.ppp2,
    process.ppp3,
    process.ppp4,
    process.ppp5,
    process.ppp6,
    process.pa1,
    process.pa5,
    process.pb5,
    process.ptrk,
    process.pevt,
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
# follow the exactly same config of process.eventinfoana.eventFilterNames
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
