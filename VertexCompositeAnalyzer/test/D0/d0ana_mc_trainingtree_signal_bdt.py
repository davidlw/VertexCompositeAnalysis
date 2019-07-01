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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/pPb_Skim_D0Both_v1/180612_110439/0000/pPb_HM_96.root',
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/64778055-CE9B-E711-A567-0025904C66F4.root'
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana_mc.root')
                                   )

process.d0ana_mc.useAnyMVA = cms.bool(True)
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")
process.d0ana_mc_wrongsign = process.d0ana_mc.clone()
process.d0ana_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCWS:D0")
process.d0ana_mc_wrongsign.MVACollection = cms.InputTag("d0selectorMCWS:MVAValuesNewD0")

#process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5_v1.root')
#process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v1.root')
process.d0selectorMCBDTPreCut.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InpPb_default_HLT185_SB_Pt1p5_v1.root')
process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
process.d0selectorMCWS = process.d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC = process.d0selectorMC.clone()
#process.npd0selectorMC.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5_v1.root')
#process.npd0selectorMC.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_WS_Pt1p5MassPeak_v1.root')
process.npd0selectorMC.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_SB_Pt1p5_v1.root')
process.npd0selectorMCWS = process.npd0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0selectorMC1 = process.d0selectorMC.clone()
#process.npd0selectorMC1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5_v1.root')
process.npd0selectorMC1.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InpPb_default_HLT185_RS_Pt1p5MassPeak_v1.root')
process.npd0selectorMCWS1 = process.npd0selectorMC1.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

process.npd0ana_mc = process.d0ana_mc.clone()
process.npd0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC:D0")
process.npd0ana_mc.MVACollection = cms.InputTag("npd0selectorMC:MVAValuesNewD0")
process.npd0ana_mc_wrongsign = process.d0ana_mc_wrongsign.clone()
process.npd0ana_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCWS:D0")
process.npd0ana_mc_wrongsign.MVACollection = cms.InputTag("npd0selectorMCWS:MVAValuesNewD0")

process.npd0ana1_mc = process.d0ana_mc.clone()
process.npd0ana1_mc.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMC1:D0")
process.npd0ana1_mc.MVACollection = cms.InputTag("npd0selectorMC1:MVAValuesNewD0")
process.npd0ana1_mc_wrongsign = process.d0ana_mc_wrongsign.clone()
process.npd0ana1_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("npd0selectorMCWS1:D0")
process.npd0ana1_mc_wrongsign.MVACollection = cms.InputTag("npd0selectorMCWS1:MVAValuesNewD0")

process.d0ana_seq = cms.Sequence(process.d0selectorMC * process.d0ana_mc)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorMCWS * process.d0ana_mc_wrongsign)
process.npd0ana_seq = cms.Sequence(process.npd0selectorMC * process.npd0ana_mc)
process.npd0ana_wrongsign_seq = cms.Sequence(process.npd0selectorMCWS * process.npd0ana_mc_wrongsign)
process.npd0ana1_seq = cms.Sequence(process.npd0selectorMC1 * process.npd0ana1_mc)
process.npd0ana1_wrongsign_seq = cms.Sequence(process.npd0selectorMCWS1 * process.npd0ana1_mc_wrongsign)

process.p1 = cms.Path(process.d0ana_seq)
process.p2 = cms.Path(process.d0ana_wrongsign_seq)
process.p3 = cms.Path(process.npd0ana_seq)
process.p4 = cms.Path(process.npd0ana_wrongsign_seq)
#process.p5 = cms.Path(process.npd0ana1_seq)
#process.p6 = cms.Path(process.npd0ana1_wrongsign_seq)
