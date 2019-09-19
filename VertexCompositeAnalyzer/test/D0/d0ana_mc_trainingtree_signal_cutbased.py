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

process.d0ana_mc.useAnyMVA = cms.bool(False)
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc_wrongsign = process.d0ana_mc.clone()
process.d0ana_mc_wrongsign.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCWS:D0")

#process.d0selectorMC = process.d0selectorMCBDTPreCut.clone()
#process.d0selectorMC.cand3DDecayLengthSigMin = cms.untracked.double(3.5)
#process.d0selectorMC.cand3DPointingAngleMax = cms.untracked.double(0.3)
process.d0selectorMC.useAnyMVA = cms.bool(False)
process.d0selectorMCWS = process.d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
)

process.d0ana_seq = cms.Sequence(process.d0selectorMC * process.d0ana_mc)
process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorMCWS * process.d0ana_mc_wrongsign)

process.p1 = cms.Path(process.d0ana_seq)
process.p2 = cms.Path(process.d0ana_wrongsign_seq)
