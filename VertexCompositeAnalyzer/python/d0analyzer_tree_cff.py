import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cfi import *

d0ana_wrongsign = d0ana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0")
)

d0ana_wrongsign_mc = d0ana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0")
)