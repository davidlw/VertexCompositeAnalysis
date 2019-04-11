import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

HybridSoftIdReco2018 = "(isGlobalMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0)"
TightIdReco = "(isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0)"

generalMuMuMassMin7Candidates = generalDiMuCandidates.clone(
    mllCutMin = cms.double(7.0),
    muonSelection = cms.string("(p > 3.5 && abs(eta) < 2.4) && ("+HybridSoftIdReco2018+" || "+TightIdReco+")"),
    tkdXYCut = cms.double(0.3),
    tkdZCut = cms.double(20.0),
    candidateSelection = cms.string("mass > 7.0")
)
