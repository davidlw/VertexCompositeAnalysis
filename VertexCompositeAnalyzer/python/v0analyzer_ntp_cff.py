import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.v0analyzer_ntp_cfi import *

lamana = ksana.clone(
  PID = cms.untracked.int32(310),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(211),
  VertexCompositeCollection = cms.untracked.InputTag("generalV0CandidatesNew:Lambda")
)

lamana_mc = ksana_mc.clone(
  PID = cms.untracked.int32(310),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(211),
  VertexCompositeCollection = cms.untracked.InputTag("generalV0CandidatesNew:Lambda")
)

xiana = ksana.clone(
  twoLayerDecay = cms.untracked.bool(True),
  PID = cms.untracked.int32(3312),
  PID_dau1 = cms.untracked.int32(3122),
  PID_dau2 = cms.untracked.int32(211),
  VertexCompositeCollection = cms.untracked.InputTag("generalCascadeCandidatesNew:Xi")
)

xiana_mc = ksana.clone(
  twoLayerDecay = cms.untracked.bool(True),
  PID = cms.untracked.int32(3312),
  PID_dau1 = cms.untracked.int32(3122),
  PID_dau2 = cms.untracked.int32(211),
  VertexCompositeCollection = cms.untracked.InputTag("generalCascadeCandidatesNew:Xi")
)

omana = ksana.clone(
  twoLayerDecay = cms.untracked.bool(True),
  PID = cms.untracked.int32(3334),
  PID_dau1 = cms.untracked.int32(3122),
  PID_dau2 = cms.untracked.int32(321),
  VertexCompositeCollection = cms.untracked.InputTag("generalCascadeCandidatesNew:Omega")
)

omana_mc = ksana.clone(
  twoLayerDecay = cms.untracked.bool(True),
  PID = cms.untracked.int32(3334),
  PID_dau1 = cms.untracked.int32(3122),
  PID_dau2 = cms.untracked.int32(321),
  VertexCompositeCollection = cms.untracked.InputTag("generalCascadeCandidatesNew:Omega")
)
