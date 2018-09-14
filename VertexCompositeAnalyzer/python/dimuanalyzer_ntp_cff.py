import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_ntp_cfi import *

dimucontana = dimuana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightPFCandidates:DiMu")
)
dimucontana_mc = dimuana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightPFCandidates:DiMu")
)

dimucontana1 = dimucontana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimPFCandidates:DiMu")
)
dimucontana1_mc = dimucontana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimPFCandidates:DiMu")
)

dimucontana2 = dimucontana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightCandidates:DiMu")
)
dimucontana2_mc = dimucontana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightCandidates:DiMu")
)

dimucontana_wrongsign = dimuana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightPFCandidatesWrongSign:DiMu")
)
dimucontana_wrongsign_mc = dimuana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightPFCandidatesWrongSign:DiMu")
)

dimucontana1_wrongsign = dimucontana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimPFCandidatesWrongSign:DiMu")
)
dimucontana1_wrongsign_mc = dimucontana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimPFCandidatesWrongSign:DiMu")
)

dimucontana2_wrongsign = dimucontana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightCandidatesWrongSign:DiMu")
)
dimucontana2_wrongsign_mc = dimucontana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalMuMuContinuimOneStTightCandidatesWrongSign:DiMu")
)

jpsiana = dimuana.clone(
  massHistPeak = cms.untracked.double(3.1),
  massHistWidth = cms.untracked.double(0.5),
  massHistBins = cms.untracked.int32(100)
)
jpsiana_mc = dimuana_mc.clone(
  massHistPeak = cms.untracked.double(3.1),
  massHistWidth = cms.untracked.double(0.5),
  massHistBins = cms.untracked.int32(100)
)

jpsiana1 = jpsiana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:DiMu")
)
jpsiana1_mc = jpsiana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidates:DiMu")
)

jpsiana2 = jpsiana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:DiMu")
)
jpsiana2_mc = jpsiana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidates:DiMu")
)

jpsiana_wrongsign = jpsiana.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:DiMu")
)
jpsiana_wrongsign_mc = jpsiana_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightPFCandidatesWrongSign:DiMu")
)

jpsiana1_wrongsign = jpsiana1.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidatesWrongSign:DiMu")
)
jpsiana1_wrongsign_mc = jpsiana1_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuPFCandidatesWrongSign:DiMu")
)

jpsiana2_wrongsign = jpsiana2.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:DiMu")
)
jpsiana2_wrongsign_mc = jpsiana2_mc.clone(
  VertexCompositeCollection = cms.untracked.InputTag("generalJPsiMuMuOneStTightCandidatesWrongSign:DiMu")
)
