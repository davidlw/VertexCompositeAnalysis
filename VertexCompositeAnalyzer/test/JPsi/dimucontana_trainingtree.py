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
'root://cms-xrd-global.cern.ch//store/user/davidlw/PASingleMuon/pPb_Skim_DiMuContBoth_default_v1/181103_103656/0000/pPb_HM_DiMuCont_99.root',
                ),
secondaryFileNames = cms.untracked.vstring(
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/8EAA7A20-89B0-E611-B231-02163E0146A5.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/A8011762-8AB0-E611-A142-02163E014455.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/AAEDB7DD-89B0-E611-AF57-FA163E30331F.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/B0E5D1A6-8AB0-E611-816C-02163E013749.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/C4AB91D2-88B0-E611-A8E9-02163E0134CE.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/DA6C982B-89B0-E611-A6E6-02163E01439F.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/E6113F5F-8AB0-E611-847E-FA163E314C09.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/F6B5DBEE-8BB0-E611-9ACD-02163E011CAC.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/4C573EFF-8AB0-E611-B795-02163E014343.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/4E263033-8AB0-E611-ACC2-02163E014210.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/705E7FD0-89B0-E611-8A96-02163E01441F.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/7EE00112-8BB0-E611-BF96-02163E01455B.root',
'/store/hidata/PARun2016C/PASingleMuon/AOD/PromptReco-v1/000/285/530/00000/3258AC9D-88B0-E611-9964-FA163E6337E8.root'
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

process.dimucontana_seq = cms.Sequence(process.dimucontana)
process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)

process.p = cms.Path(process.dimucontana_seq)
process.p1 = cms.Path(process.dimucontana_wrongsign_seq)
