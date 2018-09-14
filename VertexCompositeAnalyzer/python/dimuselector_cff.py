import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cfi import *
jpsiselector = dimuselector.clone()
jpsiselectorMC = dimuselectorMC.clone()

jpsiselectorWS = jpsiselector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:MVAValues")
)

jpsiselectorWSMC = jpsiselectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:MVAValues")
)

jpsiselectorMCGenMatch = jpsiselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

jpsiselectorMCGenUnMatch = jpsiselectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

jpsiselectorCut = jpsiselector.clone(
  useAnyMVA = cms.bool(False),

  trkPMin = cms.untracked.double(3.),
  trkNHitMin = cms.untracked.int32(6),
  candVtxProbMin = cms.untracked.double(0.01)
)

########## PF Candidates only
jpsiselector1 = jpsiselector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuPFCandidates:MVAValues")
)

jpsiselector1MC = jpsiselectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuPFCandidates:MVAValues")
)

jpsiselector1WS = jpsiselector1.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuPFCandidatesWrongSign:MVAValues")
)

jpsiselector1WSMC = jpsiselector1MC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuPFCandidatesWrongSign:MVAValues")
)

jpsiselector1MCGenMatch = jpsiselector1MC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

jpsiselector1MCGenUnMatch = jpsiselector1MC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

########## one station tight only
jpsiselector2 = jpsiselector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightCandidates:MVAValues")
)

jpsiselector2MC = jpsiselectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightCandidates:MVAValues")
)

jpsiselector2WS = jpsiselector2.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:MVAValues")
)

jpsiselector2WSMC = jpsiselector2MC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:DiMu"),
  MVACollection = cms.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:MVAValues")
)

jpsiselector2MCGenMatch = jpsiselector2MC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

jpsiselector2MCGenUnMatch = jpsiselector2MC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

