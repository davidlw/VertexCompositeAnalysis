import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cfi import *

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

  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.5),
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

  trkPtSumMin = cms.untracked.double(1.6),
  trkEtaDiffMax = cms.untracked.double(1.5),
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
