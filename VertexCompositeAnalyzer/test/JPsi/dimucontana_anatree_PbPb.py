import FWCore.ParameterSet.Config as cms

process = cms.Process("dimucontana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(200)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) 
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "103X_dataRun2_Prompt_v3"

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/user/davidlw/HIDoubleMuon/2018Skim_DiMuCont_MuonPhysics_v2_NoEvtSel/190128_232428/0000/PbPb_DiMuCont_99.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/CCD474D5-F2A6-1548-991F-F8F6CB4BF6D3.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/C121619F-F2A7-7545-8389-5F81D95377AA.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/B92F8304-76FA-FF4D-89C3-92FB8A3BCF00.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/B8B4E278-E587-144D-A953-D8C5220843EF.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/AB31BA6C-81F5-1544-A6B6-04D96CFB6B1C.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/F6CFA5C5-BFB7-394E-8F39-590A1ABE9638.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/F963B5B5-8C6C-7346-B0C6-E362E0D35B85.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/DFE46C89-8497-9849-B8FA-E95957A0734C.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/DB4FD845-293D-CE42-AC7D-7D3596F9FFFC.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/1C632DD7-B80A-0745-A83B-76DEAF9A5222.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/2643A583-7C1E-4648-BCB5-9D99021B0056.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/47FC0F6A-DDDD-924B-8439-9D43462A63D0.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/4EBB00E9-CCE5-9D47-BC54-CC3E4003CAE7.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/57977A2C-ED36-A043-87FB-C182E1CF2542.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/60330E2C-4569-7247-993D-670E0D92A5D0.root',
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/004/00000/66954445-5658-5142-9209-CB3CF65652CE.root',
)
                            )

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHM.HLTPaths = ['HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v*']
process.hltHM.HLTPaths = ['HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('dimucontana_training.root')
                                   )

process.dimucontana.isCentrality = cms.bool(True)

process.dimucontana_seq = cms.Sequence(process.dimucontana)

process.p = cms.Path(process.hltHM * process.dimucontana_seq)
#process.p1 = cms.Path(process.dimucontana_wrongsign_seq)


from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
