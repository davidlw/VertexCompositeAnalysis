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
'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_10_4_0_mtd3/src/VertexCompositeAnalysis/VertexCompositeProducer/test/PbPb.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/user/anstahll/MTD/MC/NonEmbedded/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190109/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190109/190109_203806/0000/Pythia8_TuneCP5_5TeV_D0_PiK_prompt_RECO_15.root'
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana_mc.root')
                                   )

process.d0ana_mc_genmatch = process.d0ana_mc.clone()
process.d0ana_mc_genunmatch = process.d0ana_mc.clone()
process.d0ana_mc_genmatchswap = process.d0ana_mc.clone()
process.d0ana_mc_genmatchunswap = process.d0ana_mc.clone()

process.d0ana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatch:D0")
process.d0ana_mc_genunmatch.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenUnMatch:D0")
process.d0ana_mc_genmatchswap.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatchSwap:D0")
process.d0ana_mc_genmatchunswap.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatchUnSwap:D0")
process.d0ana_wrongsign_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWSMC:D0")

process.d0ana_genmatch_seq = cms.Sequence(process.d0selectorMCGenMatch * process.d0ana_mc_genmatch)
process.d0ana_genunmatch_seq = cms.Sequence(process.d0selectorMCGenUnMatch * process.d0ana_mc_genunmatch)
process.d0ana_genmatchswap_seq = cms.Sequence(process.d0selectorMCGenMatchSwap * process.d0ana_mc_genmatchswap)
process.d0ana_genmatchunswap_seq = cms.Sequence(process.d0selectorMCGenMatchUnSwap * process.d0ana_mc_genmatchunswap)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorWSMC * process.d0ana_wrongsign)

#process.p1 = cms.Path(process.d0ana_genmatch_seq)
process.p2 = cms.Path(process.d0ana_genmatchswap_seq)
process.p3 = cms.Path(process.d0ana_genmatchunswap_seq)
process.p4 = cms.Path(process.d0ana_genunmatch_seq)
