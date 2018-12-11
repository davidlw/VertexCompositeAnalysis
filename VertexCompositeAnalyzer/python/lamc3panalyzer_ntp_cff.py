import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cfi import *

lamc3pana_wrongsign = lamc3pana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalLamC3PCandidatesNewWrongSign:LamC3P"),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNewWrongSign:MVAValues")
)

lamc3pana_wrongsign_mc = lamc3pana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalLamC3PCandidatesNewWrongSign:LamC3P"),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNewWrongSign:MVAValues")
)

lamc3pana_tof = lamc3pana.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

lamc3pana_tof_wrongsign = lamc3pana_wrongsign.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

lamc3pana_tof_mc = lamc3pana_mc.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)

lamc3pana_tof_wrongsign_mc = lamc3pana_wrongsign_mc.clone(
  doGenMatchingTOF = cms.untracked.bool(True)
)
