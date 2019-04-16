import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cfi import *

dimucontana = dimuana.clone(
    massHistPeak = cms.untracked.double(62.),
    massHistWidth = cms.untracked.double(60.),
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2Candidates:DiMu")
)
dimucontana_mc = dimuana_mc.clone(
    massHistPeak = cms.untracked.double(62.),
    massHistWidth = cms.untracked.double(60.),
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2Candidates:DiMu")
)

dimucontana_wrongsign = dimuana.clone(
    massHistPeak = cms.untracked.double(62.),
    massHistWidth = cms.untracked.double(60.),
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2CandidatesWrongSign:DiMu")
)
dimucontana_wrongsign_mc = dimuana_mc.clone(
    massHistPeak = cms.untracked.double(62.),
    massHistWidth = cms.untracked.double(60.),
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2CandidatesWrongSign:DiMu")
)
