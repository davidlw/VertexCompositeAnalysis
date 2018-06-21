import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_D0Both_v1_test/180523_221148/0000/pPb_HM_81.root',
'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_8_0_28_jpsi/src/VertexCompositeAnalysis/VertexCompositeProducer/test/pPb_HM.root'
                ),
secondaryFileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/0C73EB2B-22AF-E611-85F1-02163E013657.root'
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
)
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana_wrongsign.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWS:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana_wrongsign.MVACollection = cms.InputTag("d0selectorWS:MVAValuesNewD0")
process.d0ana.isSkimMVA = cms.untracked.bool(True)
process.d0ana_wrongsign.isSkimMVA = cms.untracked.bool(True)

process.d0selector.useAnyMVA = cms.bool(True)
process.d0selectorWS.useAnyMVA = cms.bool(True)
process.d0selector.mvaMin = cms.untracked.double(0.2)
process.d0selectorWS.mvaMin = cms.untracked.double(0.2)

process.d0ana_seq = cms.Sequence(process.d0selector * process.d0ana)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorWS * process.d0ana_wrongsign)

process.p = cms.Path(process.d0ana_seq)
process.p1 = cms.Path(process.d0ana_wrongsign_seq)
