import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('d0ana',eras.Run2_2016_pA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_dataRun2_v19')

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/user/davidlw/PAHighMultiplicity1/Pbp_Skim_D0Both_default_v1/180827_205121/0000/pPb_HM_685.root',
                ),
secondaryFileNames = cms.untracked.vstring(
'/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/301/00000/6A955CD0-7BBA-E611-BB15-02163E011C00.root',
)
                            )

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
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)

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
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v1.root')
process.d0selector.multMin = cms.untracked.double(0)
process.d0selector.multMax = cms.untracked.double(100000)
process.d0selectorWS = process.d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selector = process.d0selector.clone()
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v1.root')
process.npd0selectorWS = process.d0selectorWS.clone()
process.npd0selectorWS.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v1.root')

process.npd0selector1 = process.d0selector.clone()
process.npd0selector1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_v1.root')
process.npd0selectorWS1 = process.d0selectorWS.clone()
process.npd0selectorWS1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_v1.root')

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
process.Flag_hfCoincFilter = cms.Path(process.eventFilter_HM * process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.eventFilter_HM * process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.eventFilter_HM * process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.eventFilter_HM * process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.eventFilter_HM * process.pileUpFilter_pPb8TeV_Gplus)
#eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
eventFilterPaths = [ process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
