import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cfi import *

d0selectorBDTPreCut = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkNHitMin = cms.untracked.int32(11),
)

d0selectorBDTNonPrompt = d0selectorBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  mvaCuts = cms.vdouble(0.45,-0.0047,-0.023,0.23,-0.087,0.011,-0.0005)
)

d0selectorBDTPrompt = d0selectorBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  mvaCuts = cms.vdouble(0.4,0.19,-0.088,0.14,-0.0054,-0.0016,0.0001)
)

d0selectorMCBDTPreCut = d0selectorMC.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkNHitMin = cms.untracked.int32(11),
)

d0selectorMCBDTNonPrompt = d0selectorMCBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  mvaCuts = cms.untracked.vdouble(0.45,-0.0046,-0.023,0.28,-0.083,0.01,-0.0005)
)

d0selectorMCBDTPrompt = d0selectorMCBDTPreCut.clone(
  useAnyMVA = cms.bool(True),
  mvaCuts = cms.untracked.vdouble(0.42,0.21,-0.09,0.1,-0.006,0,0)
)     

d0selectorPID = d0selector.clone(
  useAnyMVA = cms.bool(False),
  userPID = cms.bool(True),

  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorPID2 = d0selector.clone(
  useAnyMVA = cms.bool(False),
  userPID = cms.bool(True),

  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.2),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCut = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutNew = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutNew2 = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)


d0selectorCutNew3 = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.7),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutNew4 = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.8),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutNew5 = d0selector.clone(
  useAnyMVA = cms.bool(False),
  
  trkPtSumMin = cms.untracked.double(1.6),
  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutMC = d0selector.clone(
  useAnyMVA = cms.bool(False),
  
  trkPtMin = cms.untracked.double(0.7),
  trkEtaMax = cms.untracked.double(1.5),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorCutNewMC = d0selector.clone(
  useAnyMVA = cms.bool(False),

  trkPtMin = cms.untracked.double(0.7),
  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.),
  trkPtErrMax = cms.untracked.double(0.1),
  trkNHitMin = cms.untracked.int32(11),
  cand3DDecayLengthSigMin = cms.untracked.double(3.5),
  cand3DPointingAngleMax = cms.untracked.double(0.15),
  candVtxProbMin = cms.untracked.double(0.15)
)

d0selectorWS = d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

d0selectorWSMC = d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

d0selectorMCGenMatch = d0selectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

d0selectorMCGenUnMatch = d0selectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

d0selectorMCGenMatchSwap = d0selectorMC.clone(
  selectGenMatchSwap = cms.untracked.bool(True)
)

d0selectorMCGenMatchUnSwap = d0selectorMC.clone(
  selectGenMatchUnSwap = cms.untracked.bool(True)
)
