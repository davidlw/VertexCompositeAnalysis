import FWCore.ParameterSet.Config as cms

process = cms.Process("dimucontana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) 
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "103X_dataRun2_Prompt_v3"
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/user/davidlw/HIDoubleMuon/Skim_DiMuContBoth_DCSOnly327327_v1/181127_104136/0000/PbPb_DiMuCont_13.root '
                ),
secondaryFileNames = cms.untracked.vstring(
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/816A5339-7E58-9F4E-85B2-64981FF42FB6.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/9DFA327B-02F6-9849-9169-BDAF1057E9D2.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/551C06DD-8A73-D944-8CA1-6B812CC3CCD9.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/5E760C32-A180-FB4B-976F-5D42972F20CF.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/163B4299-C86C-E74A-99FF-36365029CCDF.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/1CE70355-A1BC-E649-82E8-2CDB8B2154D0.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/06D10739-8DE8-504C-B872-12449FEAE4E5.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/10E3662B-BB42-E94B-A012-D8446CE44885.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/02468590-6F76-F24E-9022-E8F6BDF4076F.root',
'/store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/326/941/00000/064283D3-7F97-D248-B002-DC100F9456C0.root'
)
                            )

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(False),
  useJobReport = cms.untracked.bool(True)
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
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.hiCentrality.produceHFhits = cms.bool(False)
process.hiCentrality.produceEcalhits = cms.bool(False)
process.hiCentrality.produceZDChits = cms.bool(False)
process.hiCentrality.produceETmidRapidity = cms.bool(False)
process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceTracks = cms.bool(False)
process.hiCentrality.producePixelTracks = cms.bool(False)

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

process.dimucontana.isCentrality = cms.bool(True)
process.dimucontana_wrongsign.isCentrality = cms.bool(True)

process.dimucontana_seq = cms.Sequence(process.dimucontana)
process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)

process.p = cms.Path(process.cent_seq * process.dimucontana_seq)
process.p1 = cms.Path(process.cent_seq * process.dimucontana_wrongsign_seq)

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
