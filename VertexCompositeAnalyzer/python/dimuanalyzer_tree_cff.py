import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cfi import *

dimucontana = dimuana.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin7Candidates:DiMu")
)
dimucontana_mc = dimuana_mc.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidates:DiMu")
)

dimucontana_wrongsign = dimuana.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin7CandidatesWrongSign:DiMu")
)
dimucontana_wrongsign_mc = dimuana.clone(
     VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidatesWrongSign:DiMu")
)
