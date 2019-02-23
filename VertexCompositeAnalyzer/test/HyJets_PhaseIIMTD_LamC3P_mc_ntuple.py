import FWCore.ParameterSet.Config as cms

process = cms.Process("lamc3pana")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'file:/eos/cms/store/group/phys_heavyions/davidlw/hyjets_lambdac.root'
                ),
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('hyjets_mc_mtd.root')
                                   )

process.lamc3pana_mc.isUseMtd = cms.untracked.bool(True);
process.lamc3pana_mc.doRecoNtuple = cms.untracked.bool(True);
process.lamc3pana_mc.doGenNtuple = cms.untracked.bool(True);
process.lamc3pana_mc.doGenMatching = cms.untracked.bool(False);
process.lamc3pana_mc.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")
process.lamc3pana_mc.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMC:LamC3P")
process.lamc3pana_mc.MVACollection = cms.InputTag("lamc3pselectorMC:MVAValuesNewLamC3P")
process.lamc3pana_mc.isCentrality = cms.bool(True)

process.lamc3pselectorMC.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")

process.lamc3pana_seq = cms.Sequence(process.lamc3pselectorMC * process.lamc3pana_mc)

process.p = cms.Path(process.lamc3pana_seq)
