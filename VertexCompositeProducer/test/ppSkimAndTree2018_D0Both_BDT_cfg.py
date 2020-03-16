import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#check
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/data/Run2018C/HighMultiplicityEOF1/AOD/PromptReco-v2/000/319/488/00000/FC9058C3-9987-E811-9CB5-FA163E4D501C.root'
#'/store/data/Run2017C/HighMultiplicityEOF1/AOD/PromptReco-v2/000/299/929/00000/345B3236-0876-E711-92A8-02163E014169.root'
),
          dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
          inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_*totem*_*_*'

   )
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(600))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v11'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_FullTrack_Multiplicity85_*', # High multiplicity
    'HLT_FullTrack_Multiplicity100_*', # High multiplicity
    'HLT_FullTrack_Multiplicity130_*', # High multiplicity
    'HLT_FullTrack_Multiplicity155_*', # High multiplicity
    'HLT_L1MinimumBiasHF_OR_*', # Minimum bias
    'HLT_ZeroBias_*', # Zero bias
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
#remove the default dz1p0 filter and hfCoincFilter
process.primaryVertexFilterPA = process.primaryVertexFilter.clone()
process.colEvtSel = cms.Sequence(process.primaryVertexFilterPA * process.NoScraping)

process.eventFilter_HM = cms.Sequence(
    process.hltFilter*
    process.colEvtSel
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

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

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew )
process.d0rereco_wrongsign_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNewWrongSign )

###############################################################################################

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pp.root'),
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
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v3.root')
process.d0selector.GBRForestLabel = cms.string('D0Inpp')
process.d0selector.multMin = cms.untracked.double(0)
process.d0selector.multMax = cms.untracked.double(100000)
process.d0selector.mvaMax = cms.untracked.double(999.9)
process.d0selector.mvaMin = cms.untracked.double(0.4)
process.d0selectorWS = process.d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selector = process.d0selector.clone()
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v3.root')
process.npd0selectorWS = process.d0selectorWS.clone()
process.npd0selectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v3.root')
process.npd0selector1 = process.d0selector.clone()
process.npd0selector1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_v3.root')
process.npd0selectorWS1 = process.d0selectorWS.clone()
process.npd0selectorWS1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_v3.root')

process.d0selectorNew = process.d0selector.clone()
process.npd0selectorNew = process.npd0selector.clone()
process.npd0selector1New = process.npd0selector1.clone()
process.d0selectorNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selectorNew.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selector1New.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.d0selectorWSNew = process.d0selectorWS.clone()
process.npd0selectorWSNew = process.npd0selectorWS.clone()
process.npd0selectorWS1New = process.npd0selectorWS1.clone()
process.d0selectorWSNew.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWSNew.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')
process.npd0selectorWS1New.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_NoPtErrNHitDLAngle2D_v3.root')

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.multMin = cms.untracked.double(0)
process.d0ana.multMax = cms.untracked.double(100000)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana_wrongsign = process.d0ana.clone()
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

process.d0ana_new = process.d0ana.clone()
process.npd0ana_new = process.npd0ana.clone()
process.npd0ana1_new = process.npd0ana1.clone()
process.d0ana_new.VertexCompositeCollection = cms.untracked.InputTag("d0selectorNew:D0")
process.d0ana_new.MVACollection = cms.InputTag("d0selectorNew:MVAValuesNewD0")
process.npd0ana_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorNew:D0")
process.npd0ana_new.MVACollection = cms.InputTag("npd0selectorNew:MVAValuesNewD0")
process.npd0ana1_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selector1New:D0")
process.npd0ana1_new.MVACollection = cms.InputTag("npd0selector1New:MVAValuesNewD0")
process.d0ana_wrongsign_new = process.d0ana_wrongsign.clone()
process.npd0ana_wrongsign_new = process.npd0ana_wrongsign.clone()
process.npd0ana1_wrongsign_new = process.npd0ana1_wrongsign.clone()
process.d0ana_wrongsign_new.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWSNew:D0")
process.d0ana_wrongsign_new.MVACollection = cms.InputTag("d0selectorWSNew:MVAValuesNewD0")
process.npd0ana_wrongsign_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWSNew:D0")
process.npd0ana_wrongsign_new.MVACollection = cms.InputTag("npd0selectorWSNew:MVAValuesNewD0")
process.npd0ana1_wrongsign_new.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorWS1New:D0")
process.npd0ana1_wrongsign_new.MVACollection = cms.InputTag("npd0selectorWS1New:MVAValuesNewD0")

process.d0ana_seq = cms.Sequence(process.eventFilter_HM * process.d0selector * process.d0ana)
process.npd0ana_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector * process.npd0ana)
process.npd0ana1_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector1 * process.npd0ana1)
process.d0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.d0selectorWS * process.d0ana_wrongsign)
process.npd0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS * process.npd0ana_wrongsign)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1 * process.npd0ana1_wrongsign)
process.d0ana_seq1 = cms.Sequence(process.eventFilter_HM * process.d0selectorNew * process.d0ana_new)
process.npd0ana_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorNew * process.npd0ana_new)
process.npd0ana1_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selector1New * process.npd0ana1_new)
process.d0ana_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.d0selectorWSNew * process.d0ana_wrongsign_new)
process.npd0ana_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWSNew * process.npd0ana_wrongsign_new)
process.npd0ana1_wrongsign_seq1 = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1New * process.npd0ana1_wrongsign_new)

# eventinfoana must be in EndPath, and process.eventinfoana.selectEvents must be the name of eventFilter_HM Path
process.eventinfoana.selectEvents = cms.untracked.string('eventFilter_HM_step')
process.pevt = cms.EndPath(process.eventinfoana)

process.pp = cms.Path(process.d0ana_seq)
process.pp1 = cms.Path(process.npd0ana_seq)
process.pp2 = cms.Path(process.npd0ana1_seq)
process.pp_ws = cms.Path(process.d0ana_wrongsign_seq)
process.pp1_ws = cms.Path(process.npd0ana_wrongsign_seq)
process.pp2_ws = cms.Path(process.npd0ana1_wrongsign_seq)

process.ppp = cms.Path(process.d0ana_seq1)
process.ppp1 = cms.Path(process.npd0ana_seq1)
process.ppp2 = cms.Path(process.npd0ana1_seq1)
process.ppp_ws = cms.Path(process.d0ana_wrongsign_seq1)
process.ppp1_ws = cms.Path(process.npd0ana_wrongsign_seq1)
process.ppp2_ws = cms.Path(process.npd0ana1_wrongsign_seq1)

# Add the Conversion tree
process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.track_cff")
process.ptrk = cms.Path(process.eventFilter_HM * process.track_ana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
#    process.d0rereco_wrongsign_step,
    process.pp,
    process.pp1,
    process.pp2,
#    process.pp_ws,
#    process.pp1_ws,
#    process.pp2_ws,
    process.ppp,
    process.ppp1,
    process.ppp2,
#    process.ppp_ws,
#    process.ppp1_ws,
#    process.ppp2_ws,
    process.ptrk,
    process.pevt,
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
# no hfCoincFilter available here
#process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
# here, if one want to use other config of olvFilter, Gplus filter, one need to change process.eventinfoana.eventFilterNames
process.Flag_pileupVertexFilterCut_pp5TeV = cms.Path(process.eventFilter_HM * process.olvFilter_pp5TeV_dz1p0)
process.Flag_pileupVertexFilterCut_pPb8TeV = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus_pp5TeV = cms.Path(process.eventFilter_HM * process.pileUpFilter_pp5TeV_Gplus)
process.Flag_pileupVertexFilterCutGplus_pPb8TeV = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
# follow the exactly same config of process.eventinfoana.eventFilterNames
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut_pp5TeV,  process.Flag_pileupVertexFilterCut_pPb8TeV, process.Flag_pileupVertexFilterCutGplus_pp5TeV, process.Flag_pileupVertexFilterCutGplus_pPb8TeV ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
