import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

generalJPsiMuMuCandidates = generalDiMuCandidates.clone()

generalJPsiMuMuOneStTightCandidates = generalJPsiMuMuCandidates.clone(
    isMuonId = cms.bool(True),
    muonId = cms.string('TMOneStationTight')
)

generalJPsiMuMuGlobalCandidates = generalJPsiMuMuCandidates.clone(
    isMuonId = cms.bool(True),
    isGlobalMuon = cms.bool(True),
    trackQualities = cms.vstring(''),
    muonId = cms.string('TMOneStationTight')
)

generalUpsilonMuMuGlobalCandidates = generalJPsiMuMuCandidates.clone(
    isMuonId = cms.bool(True),
    isGlobalMuon = cms.bool(True),
    trackQualities = cms.vstring('loose'),
    muonId = cms.string('TMOneStationTight'),

    mllCutMin = cms.double(7.0),
    mllCutMax = cms.double(15.0),

    DiMuMassCut = cms.double(999999.),
)

generalJPsiMuMuPFCandidates = generalJPsiMuMuCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalJPsiMuMuOneStTightPFCandidates = generalJPsiMuMuOneStTightCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalMuMuContinuimCandidates = generalDiMuCandidates.clone(

    mllCutMin = cms.double(0.0),
    mllCutMax = cms.double(10000.0),

    DiMuMassCut = cms.double(999999.),
)

generalMuMuContinuimOneStTightCandidates = generalMuMuContinuimCandidates.clone(
    isMuonId = cms.bool(True),
    muonId = cms.string('TMOneStationTight')
)

generalMuMuContinuimPFCandidates = generalMuMuContinuimCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalMuMuContinuimOneStTightPFCandidates = generalMuMuContinuimOneStTightCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalMuMuContinuimGlobalCandidates = generalMuMuContinuimCandidates.clone(
    isGlobalMuon = cms.bool(True)
)

generalMuMuContinuimPFGlobalCandidates = generalMuMuContinuimPFCandidates.clone(
    isGlobalMuon = cms.bool(True)
)
