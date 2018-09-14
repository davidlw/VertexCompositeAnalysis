import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.lambdacselector_cfi import *

lambdacselectorMCGenMatch = lambdacselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)

lambdacselectorMCGenUnMatch = lambdacselectorMC.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

lambdacselectorKsP = lambdacselector.clone()
lambdacselectorLambdaPi = lambdacselector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalLambdaCCandidatesNew:LambdaCToLamPi")
)

lambdacselectorKsPMCGenMatch = lambdacselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True)
)
lambdacselectorKsPMCGenUnMatch = lambdacselectorKsPMCGenMatch.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)

lambdacselectorLamPiMCGenMatch = lambdacselectorMC.clone(
  selectGenMatch = cms.untracked.bool(True),

  VertexCompositeCollection = cms.untracked.InputTag("generalLambdaCCandidatesNew:LambdaCToLamPi"),

  PID_dau1 = cms.untracked.int32(3122),
  PID_dau2 = cms.untracked.int32(211)
)
lambdacselectorLamPiMCGenUnMatch = lambdacselectorLamPiMCGenMatch.clone(
  selectGenUnMatch = cms.untracked.bool(True)
)
