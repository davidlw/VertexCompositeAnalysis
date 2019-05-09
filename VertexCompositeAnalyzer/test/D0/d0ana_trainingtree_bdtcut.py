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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_D0Both_v7_test/180704_132822/0000/pPb_HM_401.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/4C0D189A-1BAF-E611-B5CA-FA163EAF1F45.root',
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/84DA1928-90AF-E611-A838-02163E01240B.root'
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana_training.root')
                                   )

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana_wrongsign.useAnyMVA = cms.bool(True)
process.d0ana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWS:D0")
process.d0ana_wrongsign.MVACollection = cms.InputTag("d0selectorWS:MVAValuesNewD0")

#process.d0selector = process.d0selectorBDTPrompt.clone()
process.d0selector = process.d0selectorBDTPreCut.clone()
#process.d0selector.useAnyMVA = cms.bool(False)
#process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_zVtxStep5_v3.root')
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5_v1.root')
process.d0selector.multMin = cms.untracked.double(185)
process.d0selector.multMax = cms.untracked.double(250)
#process.d0selector.mvaMax = cms.untracked.double(999.9)
#process.d0selector.mvaMin = cms.untracked.double(0.3)

process.d0selectorWS = process.d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selector = process.d0selector.clone()
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5_v1.root')
process.npd0selectorWS = process.npd0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selector1 = process.d0selector.clone()
process.npd0selector1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5_v1.root')
process.npd0selectorWS1 = process.npd0selector1.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

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

process.p = cms.Path(process.d0ana_seq)
process.p1 = cms.Path(process.d0ana_wrongsign_seq)
process.p2 = cms.Path(process.npd0ana_seq)
process.p3 = cms.Path(process.npd0ana_wrongsign_seq)
process.p4 = cms.Path(process.npd0ana1_seq)
process.p5 = cms.Path(process.npd0ana1_wrongsign_seq)
