import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('d0ana',eras.Run2_2017)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')

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

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/HighMultiplicityEOF5/pp2017E_Skim_D0Both_default_v1/190704_014637/0000/pp_HM_6.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/F6BEC50B-8CEA-E711-812E-C45444922BB0.root',
'root://cms-xrd-global.cern.ch//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/EE03FFF7-8DEA-E711-AEE2-0CC47A706DC0.root',
'root://cms-xrd-global.cern.ch//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/C43B18C2-92EA-E711-8DAD-0025901D46E4.root',
'root://cms-xrd-global.cern.ch//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/A666580D-8CEA-E711-BAB6-0025901D46E4.root',
'root://cms-xrd-global.cern.ch//store/data/Run2017E/HighMultiplicityEOF5/AOD/17Nov2017-v1/50000/805AC504-87EA-E711-B5AC-002590DE71DA.root',
)
                            )

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Other triggers
    'HLT_FullTrack_Multiplicity85_*', # High multiplicity
    'HLT_FullTrack_Multiplicity105_*', # High multiplicity
    'HLT_FullTrack_Multiplicity155_*', # High multiplicity
    'HLT_L1MinimumBiasHF_OR_*', # Minimum bias
    'HLT_ZeroBias_*', # Zero bias
    # Dimuon triggers
    'HLT_L1DoubleMu0_v*',
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.NoScraping)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")

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

process.d0ana_seq = cms.Sequence(process.eventFilter_HM * process.d0selector * process.d0ana)
process.npd0ana_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector * process.npd0ana)
process.npd0ana1_seq = cms.Sequence(process.eventFilter_HM * process.npd0selector1 * process.npd0ana1)
process.d0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.d0selectorWS * process.d0ana_wrongsign)
process.npd0ana_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS * process.npd0ana_wrongsign)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.eventFilter_HM * process.npd0selectorWS1 * process.npd0ana1_wrongsign)

process.pevt = cms.Path(process.eventFilter_HM * process.eventinfoana)
process.pa = cms.Path(process.d0ana_seq)
process.pa1 = cms.Path(process.d0ana_wrongsign_seq)
process.pa4 = cms.Path(process.npd0ana_seq)
process.pa5 = cms.Path(process.npd0ana_wrongsign_seq)
process.pb4 = cms.Path(process.npd0ana1_seq)
process.pb5 = cms.Path(process.npd0ana1_wrongsign_seq)

# Add the Conversion tree
process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.track_cff")
process.ptrk = cms.Path(process.eventFilter_HM * process.track_ana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pevt,
    process.pa,
#    process.pa1,
#    process.pa4,
#    process.pa5,
    process.pb4,
#    process.pb5,
    process.ptrk
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pp5TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pp5TeV_Gplus)
eventFilterPaths = [ process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
