import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cfi import *

d0ana_wrongsign = d0ana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

d0ana_wrongsign_mc = d0ana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNewWrongSign:D0"),
  MVACollection = cms.InputTag("generalD0CandidatesNewWrongSign:MVAValues")
)

d0ana_tof = d0ana.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

d0ana_tof_wrongsign = d0ana_wrongsign.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

d0ana_tof_mc = d0ana_mc.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

d0ana_tof_wrongsign_mc = d0ana_wrongsign_mc.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

