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
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/Skim_LamC3P_v3/181212_205641/0000/PbPb_6.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_260.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_26.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_259.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_258.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_257.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_256.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_255.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_254.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_253.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_252.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_251.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_250.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_25.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_249.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_248.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_247.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_246.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_245.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_244.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_243.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_242.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_241.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_240.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_24.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_239.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_238.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_237.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_236.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_235.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/Pythia8_TuneCUETP8M1_5020GeV_PhaseII/LambdaC_PiKP_prompt_RECO_102X_test_v2/181022_104205/0000/step3_234.root',
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('lamc3pana_mc.root')
                                   )

process.lamc3pana_mc_genmatch = process.lamc3pana_mc.clone()
process.lamc3pana_mc_genunmatch = process.lamc3pana_mc.clone()
process.lamc3pana_mc_genmatchswap = process.lamc3pana_mc.clone()
process.lamc3pana_mc_genmatchunswap = process.lamc3pana_mc.clone()

process.lamc3pana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMCGenMatch:LamC3P")
process.lamc3pana_mc_genunmatch.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMCGenUnMatch:LamC3P")
process.lamc3pana_mc_genmatchswap.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMCGenMatchSwap:LamC3P")
process.lamc3pana_mc_genmatchunswap.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorMCGenMatchUnSwap:LamC3P")
process.lamc3pana_wrongsign_mc.VertexCompositeCollection = cms.untracked.InputTag("lamc3pselectorWSMC:LamC3P")

process.lamc3pana_genmatch_seq = cms.Sequence(process.lamc3pselectorMCGenMatch * process.lamc3pana_mc_genmatch)
process.lamc3pana_genunmatch_seq = cms.Sequence(process.lamc3pselectorMCGenUnMatch * process.lamc3pana_mc_genunmatch)
process.lamc3pana_genmatchswap_seq = cms.Sequence(process.lamc3pselectorMCGenMatchSwap * process.lamc3pana_mc_genmatchswap)
process.lamc3pana_genmatchunswap_seq = cms.Sequence(process.lamc3pselectorMCGenMatchUnSwap * process.lamc3pana_mc_genmatchunswap)
#process.lamc3pana_wrongsign_seq = cms.Sequence(process.lamc3pselectorWSMC * process.lamc3pana_wrongsign)

process.p1 = cms.Path(process.lamc3pana_genmatch_seq)
process.p2 = cms.Path(process.lamc3pana_genmatchswap_seq)
process.p3 = cms.Path(process.lamc3pana_genmatchunswap_seq)
