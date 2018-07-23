import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

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
'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_D0_v1/180621_163728/0000/pPb_HM_720.root'
#'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_8_0_28_jpsi/src/VertexCompositeAnalysis/VertexCompositeProducer/test/pPb_HM.root'
                ),
secondaryFileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/0C97B6B7-50B0-E611-B38F-FA163E41A46B.root'
)
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selectorCut:D0")
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveAllHistogram = cms.untracked.bool(False)
process.d0ana.saveTree = cms.untracked.bool(False)

process.d0ananew = process.d0ana.clone()
process.d0ananew.VertexCompositeCollection = cms.untracked.InputTag("d0selectorCutNew:D0")
process.d0ananew3 = process.d0ana.clone()
process.d0ananew3.VertexCompositeCollection = cms.untracked.InputTag("d0selectorCutNew3:D0")
process.d0anapid = process.d0ana.clone()
process.d0anapid.VertexCompositeCollection = cms.untracked.InputTag("d0selectorPID:D0")
process.d0anapid2 = process.d0ana.clone()
process.d0anapid2.VertexCompositeCollection = cms.untracked.InputTag("d0selectorPID2:D0")

process.d0ana_seq = cms.Sequence(process.d0selectorCut * process.d0ana)
process.d0ananew3_seq = cms.Sequence(process.d0selectorCutNew3 * process.d0ananew3)
process.d0anapid_seq = cms.Sequence(process.d0selectorPID * process.d0anapid)
process.d0anapid2_seq = cms.Sequence(process.d0selectorPID2 * process.d0anapid2)

process.p = cms.Path(process.d0ana_seq)
process.p3 = cms.Path(process.d0ananew3_seq)
process.p1 = cms.Path(process.d0anapid_seq)
process.p2 = cms.Path(process.d0anapid2_seq)
