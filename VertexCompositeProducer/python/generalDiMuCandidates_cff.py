import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

generalJPsiMuMuCandidates = generalDiMuCandidates.clone()

generalJPsiMuMuOneStTightCandidates = generalJPsiMuMuCandidates.clone(
    isMuonId = cms.bool(True),
    muonId = cms.string('TMOneStationTight')
)

generalJPsiMuMuPFCandidates = generalJPsiMuMuCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalJPsiMuMuOneStTightPFCandidates = generalJPsiMuMuOneStTightCandidates.clone(
    isPFMuon = cms.bool(True)
)

generalMuMuContinuimCandidates = generalDiMuCandidates.clone(

    mllCutMin = cms.double(0.0),
    mllCutMax = cms.double(100.0),

    DiMuMassCut = cms.double(9999.),
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
