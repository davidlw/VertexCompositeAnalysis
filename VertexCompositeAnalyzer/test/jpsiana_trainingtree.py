import FWCore.ParameterSet.Config as cms

process = cms.Process("jpsiana")

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
'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/pPb_Skim_JPsiBoth_default_v1/180827_220926/0000/pPb_HM_JPsi_27.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/4C0D189A-1BAF-E611-B5CA-FA163EAF1F45.root',
'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/6A19A7DC-92AF-E611-AB81-02163E0142DF.root',
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

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('jpsiana_training.root')
                                   )

process.jpsiana_seq = cms.Sequence(process.hltHM * process.jpsiana)
process.jpsiana1_seq = cms.Sequence(process.hltHM * process.jpsiana1)
process.jpsiana2_seq = cms.Sequence(process.hltHM * process.jpsiana2)

process.jpsiana_wrongsign_seq = cms.Sequence(process.hltHM * process.jpsiana_wrongsign)
process.jpsiana1_wrongsign_seq = cms.Sequence(process.hltHM * process.jpsiana1_wrongsign)
process.jpsiana2_wrongsign_seq = cms.Sequence(process.hltHM * process.jpsiana2_wrongsign)

#process.jpsiana_seq = cms.Sequence(process.jpsiana)
#process.jpsiana1_seq = cms.Sequence(process.jpsiana1)
#process.jpsiana2_seq = cms.Sequence(process.jpsiana2)

#process.jpsiana_wrongsign_seq = cms.Sequence(process.jpsiana_wrongsign)
#process.jpsiana1_wrongsign_seq = cms.Sequence(process.jpsiana1_wrongsign)
#process.jpsiana2_wrongsign_seq = cms.Sequence(process.jpsiana2_wrongsign)

process.p = cms.Path(process.jpsiana_seq)
process.p1 = cms.Path(process.jpsiana1_seq)
process.p2 = cms.Path(process.jpsiana2_seq)
process.p3 = cms.Path(process.jpsiana_wrongsign_seq)
process.p4 = cms.Path(process.jpsiana1_wrongsign_seq)
process.p5 = cms.Path(process.jpsiana2_wrongsign_seq)
