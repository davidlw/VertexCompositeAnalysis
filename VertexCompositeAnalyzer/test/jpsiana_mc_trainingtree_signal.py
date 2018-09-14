import FWCore.ParameterSet.Config as cms

process = cms.Process("jpsiana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/user/davidlw/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/RecoSkim2016_pPb_JPsi_v1/171218_174354/0000/pPb_HM_JPsi_1.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/00000/0408EB99-1C19-E711-9838-FA163E8D37EC.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/00000/0AD80CDA-1B19-E711-AA6D-FA163EAABFAA.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/00000/0CB0F3B9-5C19-E711-97C1-FA163EC49F79.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/00000/142D3DDC-1B19-E711-841E-FA163EBFBC7A.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/00000/1E25422B-1B19-E711-85D4-FA163ED76972.root'
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('jpsiana_mc.root')
                                   )
process.jpsiana_mc_genmatch = process.jpsiana_mc.clone()
process.jpsiana_mc_genunmatch = process.jpsiana_mc.clone()
process.jpsiana_wrongsign_mc = process.jpsiana_mc.clone()

process.jpsiselectorMCGenMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidates:JPsiMuMu")
process.jpsiselectorMCGenUnMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidates:JPsiMuMu")
process.jpsiselectorWSMC.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidates:JPsiMuMu")
process.jpsiselector1MCGenMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:JPsiMuMu")
process.jpsiselector1MCGenUnMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:JPsiMuMu")
process.jpsiselector1WSMC.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:JPsiMuMu")
process.jpsiselector2MCGenMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:JPsiMuMu")
process.jpsiselector2MCGenUnMatch.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:JPsiMuMu")
process.jpsiselector2WSMC.VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:JPsiMuMu")

process.jpsiana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("jpsiselectorMCGenMatch:JPsiMuMu")
process.jpsiana_mc_genunmatch.VertexCompositeCollection = cms.untracked.InputTag("jpsiselectorMCGenUnMatch:JPsiMuMu")
process.jpsiana_wrongsign_mc.VertexCompositeCollection = cms.untracked.InputTag("jpsiselectorWSMC:JPsiMuMu")

process.jpsiana_genmatch_seq = cms.Sequence(process.jpsiselectorMCGenMatch * process.jpsiana_mc_genmatch)
process.jpsiana_genunmatch_seq = cms.Sequence(process.jpsiselectorMCGenUnMatch * process.jpsiana_mc_genunmatch)
process.jpsiana_wrongsign_seq = cms.Sequence(process.jpsiselectorWSMC * process.jpsiana_wrongsign_mc)

process.jpsiana1_mc_genmatch = process.jpsiana_mc.clone()
process.jpsiana1_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("jpsiselector1MCGenMatch:JPsiMuMu")
process.jpsiana1_genmatch_seq = cms.Sequence(process.jpsiselector1MCGenMatch * process.jpsiana1_mc_genmatch)

process.jpsiana2_mc_genmatch = process.jpsiana_mc.clone()
process.jpsiana2_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("jpsiselector2MCGenMatch:JPsiMuMu")
process.jpsiana2_genmatch_seq = cms.Sequence(process.jpsiselector2MCGenMatch * process.jpsiana2_mc_genmatch)

process.p1 = cms.Path(process.jpsiana_genmatch_seq)
process.p2 = cms.Path(process.jpsiana1_genmatch_seq)
process.p3 = cms.Path(process.jpsiana2_genmatch_seq)
