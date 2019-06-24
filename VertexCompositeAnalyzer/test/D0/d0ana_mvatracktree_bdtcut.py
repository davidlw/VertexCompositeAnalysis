import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/user/davidlw/HighMultiplicityEOF0/pp_Skim_D0Both_default_v1/180920_103737/0000/pPb_HM_99.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'/store/data/Run2018C/HighMultiplicityEOF0/AOD/PromptReco-v2/000/319/468/00000/AE4AFC07-2D87-E811-A403-FA163E405ADB.root',
'/store/data/Run2018C/HighMultiplicityEOF0/AOD/PromptReco-v2/000/319/468/00000/5670AD4D-2787-E811-B599-FA163E296FCB.root',
'/store/data/Run2018C/HighMultiplicityEOF0/AOD/PromptReco-v2/000/319/468/00000/A49CFA74-2487-E811-8EE2-FA163EAE92EA.root'
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

# set up selectors
process.d0selector = process.d0selectorBDTPreCut.clone()
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v1.root')
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
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v1.root')
process.npd0selectorWS = process.d0selectorWS.clone()
process.npd0selectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_WS_Pt1p5MassPeak_v1.root')

process.npd0selector1 = process.d0selector.clone()
process.npd0selector1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_v1.root')
process.npd0selectorWS1 = process.d0selectorWS.clone()
process.npd0selectorWS1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0Inpp_default_pplowpuHLT100_RS_Pt1p5MassPeak_v1.root')

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

process.d0ana_seq = cms.Sequence(process.hltHM * process.d0selector * process.d0ana)
process.npd0ana_seq = cms.Sequence(process.hltHM * process.npd0selector * process.npd0ana)
process.npd0ana1_seq = cms.Sequence(process.hltHM * process.npd0selector1 * process.npd0ana1)
process.d0ana_wrongsign_seq = cms.Sequence(process.hltHM * process.d0selectorWS * process.d0ana_wrongsign)
process.npd0ana_wrongsign_seq = cms.Sequence(process.hltHM * process.npd0selectorWS * process.npd0ana_wrongsign)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.hltHM * process.npd0selectorWS1 * process.npd0ana1_wrongsign)

process.pa = cms.Path(process.d0ana_seq)
process.pa1 = cms.Path(process.d0ana_wrongsign_seq)
process.pa4 = cms.Path(process.npd0ana_seq)
process.pa5 = cms.Path(process.npd0ana_wrongsign_seq)
process.pb4 = cms.Path(process.npd0ana1_seq)
process.pb5 = cms.Path(process.npd0ana1_wrongsign_seq)

# Add the Conversion tree
process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.track_cff")
process.ptrk = cms.Path(process.hltHM * process.track_ana)
