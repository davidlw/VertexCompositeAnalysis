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

jpsiselectorBDTPrompt = jpsiselector.clone(
  useAnyMVA = cms.bool(True),
  trkPMin = cms.untracked.double(3.),
  trkNHitMin = cms.untracked.int32(6),
  trkPSumMin = cms.untracked.double(8.0),
  mvaCuts = cms.vdouble(-1.993e-01,1.222e-02,-4.700e-03,7.764e-02,-1.003e-01,3.025e-02,-2.079e-03)
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

