import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cfi import *

d0selectorWS = d0selector.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

d0selectorWSMC = d0selectorMC.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)
