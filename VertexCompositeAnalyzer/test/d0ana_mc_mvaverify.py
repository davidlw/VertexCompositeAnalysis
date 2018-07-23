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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/pPb_Skim_D0Both_v1/180612_110439/0000/pPb_HM_99.root'
                ),

secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/64FAD506-089C-E711-B8B3-FA163E8FAA46.root'
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

process.d0ana_mc_genmatch.useAnyMVA = cms.bool(True)
process.d0ana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatch:D0")
process.d0ana_mc_genmatch.MVACollection = cms.InputTag("d0selectorMCGenMatch:MVAValuesNewD0")
process.d0ana_mc_genmatch.isSkimMVA = cms.untracked.bool(True)
process.d0ana_mc_genmatch.saveHistogram = cms.untracked.bool(True)
process.d0ana_mc_genmatch.saveTree = cms.untracked.bool(False)
process.d0ana_mc_genmatch.pTBins = cms.untracked.vdouble(1.5,2.4)
process.d0ana_mc_genmatch.yBins = cms.untracked.vdouble(-1.0,1.0)

process.d0selectorMCGenMatch.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_pT15_24_scenario1.root')
process.d0selectorMCGenMatch.useAnyMVA = cms.bool(True)

process.d0ana_genmatch_seq = cms.Sequence(process.d0selectorMCGenMatch * process.d0ana_mc_genmatch)
process.d0ana_genunmatch_seq = cms.Sequence(process.d0selectorMCGenUnMatch * process.d0ana_mc_genunmatch)
process.d0ana_genmatchswap_seq = cms.Sequence(process.d0selectorMCGenMatchSwap * process.d0ana_mc_genmatchswap)
process.d0ana_genmatchunswap_seq = cms.Sequence(process.d0selectorMCGenMatchUnSwap * process.d0ana_mc_genmatchunswap)
#process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorWSMC * process.d0ana_wrongsign)

process.p3 = cms.Path(process.d0ana_genmatch_seq)
