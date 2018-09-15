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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_8160GeV_Ds_KsKaon/pPb_Skim_DsToKsK_v1/180709_052511/0000/pPb_HM_DS_1.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_102.root',
'root://cms-xrd-global.cern.ch//store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_101.root',
'root://cms-xrd-global.cern.ch//store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_100.root', 
'root://cms-xrd-global.cern.ch//store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_10.root', 
'root://cms-xrd-global.cern.ch//store/user/zhchen/Pythia8_8160GeV_Ds_KsKaon/Pythia8_8160GeV_Ds_KsKaon_AODSIM_v1-batch1/170709_180925/0000/Pythia8_8160GeV_Ds_KsKaon_step2_1.root'
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dsselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dsanalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('dsana_mc.root')
                                   )

process.dsana_mc_genmatch = process.dsana_mc.clone()
process.dsana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("dsselectorMCGenMatch:DSToKsK")
process.dsana_mc_genunmatch = process.dsana_mc.clone()
process.dsana_mc_genunmatch.VertexCompositeCollection = cms.untracked.InputTag("dsselectorMCGenUnMatch:DSToKsK")

process.dsana_genmatch_seq = cms.Sequence(process.dsselectorMCGenMatch * process.dsana_mc_genmatch)
process.dsana_genunmatch_seq = cms.Sequence(process.dsselectorMCGenUnMatch * process.dsana_mc_genunmatch)

process.p1 = cms.Path(process.dsana_genmatch_seq)
process.p2 = cms.Path(process.dsana_genunmatch_seq)
