import FWCore.ParameterSet.Config as cms

process = cms.Process("dimucontana")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_DiMuContBoth_default_v1/180911_150852/0000/pPb_HM_DiMuCont_99.root'
'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_8_0_28_jpsi/src/VertexCompositeAnalysis/VertexCompositeProducer/test/pPb_HM_DiMuCont.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/E4752A95-9DB0-E611-BF16-02163E01441F.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/D66C33BF-9AB0-E611-9FE8-FA163E9EAB09.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/9A38ACFB-A0B0-E611-904C-FA163EC60378.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/8C21991D-99B0-E611-B908-FA163EA20FBE.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/700FB2C2-9AB0-E611-9B6E-02163E01364A.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/38D04DAF-9DB0-E611-841E-02163E014418.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/530/00000/1C2D2DEC-7CB0-E611-9C54-02163E0120BA.root'
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('dimucontana_training.root')
                                   )

process.dimucontana_seq = cms.Sequence(process.hltHM * process.dimucontana)
process.dimucontana1_seq = cms.Sequence(process.hltHM * process.dimucontana1)
process.dimucontana2_seq = cms.Sequence(process.hltHM * process.dimucontana2)

process.dimucontana_wrongsign_seq = cms.Sequence(process.hltHM * process.dimucontana_wrongsign)
process.dimucontana1_wrongsign_seq = cms.Sequence(process.hltHM * process.dimucontana1_wrongsign)
process.dimucontana2_wrongsign_seq = cms.Sequence(process.hltHM * process.dimucontana2_wrongsign)

#process.dimucontana_seq = cms.Sequence(process.dimucontana)
#process.dimucontana1_seq = cms.Sequence(process.dimucontana1)
#process.dimucontana2_seq = cms.Sequence(process.dimucontana2)

#process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)
#process.dimucontana1_wrongsign_seq = cms.Sequence(process.dimucontana1_wrongsign)
#process.dimucontana2_wrongsign_seq = cms.Sequence(process.dimucontana2_wrongsign)

process.p = cms.Path(process.dimucontana_seq)
process.p1 = cms.Path(process.dimucontana1_seq)
process.p2 = cms.Path(process.dimucontana2_seq)
process.p3 = cms.Path(process.dimucontana_wrongsign_seq)
process.p4 = cms.Path(process.dimucontana1_wrongsign_seq)
process.p5 = cms.Path(process.dimucontana2_wrongsign_seq)
