import FWCore.ParameterSet.Config as cms

process = cms.Process("dimucontana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) 
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = '75X_mcRun2_HeavyIon_v7'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v12', '')

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_7_5_5_patch3/src/VertexCompositeAnalysis/VertexCompositeProducer/test/PbPb_DiMuCont.root',
                ),
secondaryFileNames = cms.untracked.vstring(
'/store/hidata/HIRun2015/HIOniaL1DoubleMu0/AOD/PromptReco-v1/000/263/757/00000/F2BD4D30-ACBB-E511-9FC2-02163E012647.root'
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

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")

process.newCentralityBin = process.centralityBin.clone()

process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.hiCentrality.produceHFhits = cms.bool(False)
process.hiCentrality.produceEcalhits = cms.bool(False)
process.hiCentrality.produceZDChits = cms.bool(False)
process.hiCentrality.produceETmidRapidity = cms.bool(False)
process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceTracks = cms.bool(False)
process.hiCentrality.producePixelTracks = cms.bool(False)

process.cent_seq = cms.Sequence(process.hiCentrality * process.newCentralityBin)

process.dimucontana.isCentrality = cms.bool(True)
process.dimucontana.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
process.dimucontana.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")
process.dimucontana.VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightGlobalCandidates:DiMu")
#process.dimucontana1.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
#process.dimucontana1.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")
#process.dimucontana2.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
#process.dimucontana2.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")
process.dimucontana_wrongsign.isCentrality = cms.bool(True)
process.dimucontana_wrongsign.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
process.dimucontana_wrongsign.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")
process.dimucontana_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightGlobalCandidatesWrongSign:DiMu")
#process.dimucontana1_wrongsign.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
#process.dimucontana1_wrongsign.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")
#process.dimucontana2_wrongsign.VertexCollection = cms.untracked.InputTag("hiSelectedVertex")
#process.dimucontana2_wrongsign.TrackCollection = cms.untracked.InputTag("hiGeneralTracks")

process.dimucontana_seq = cms.Sequence(process.dimucontana)
#process.dimucontana1_seq = cms.Sequence(process.dimucontana1)
#process.dimucontana2_seq = cms.Sequence(process.dimucontana2)
process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)
#process.dimucontana1_wrongsign_seq = cms.Sequence(process.dimucontana1_wrongsign)
#process.dimucontana2_wrongsign_seq = cms.Sequence(process.dimucontana2_wrongsign)

process.p = cms.Path(process.cent_seq * process.dimucontana_seq)
#process.p1 = cms.Path(process.dimucontana1_seq)
#process.p2 = cms.Path(process.dimucontana2_seq)
process.p3 = cms.Path(process.cent_seq * process.dimucontana_wrongsign_seq)
#process.p4 = cms.Path(process.dimucontana1_wrongsign_seq)
#process.p5 = cms.Path(process.dimucontana2_wrongsign_seq)

