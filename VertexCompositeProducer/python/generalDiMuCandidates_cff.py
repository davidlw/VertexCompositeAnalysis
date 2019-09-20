import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

SoftIdReco = "(innerTrack.isNonnull && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && innerTrack.quality(\"highPurity\"))"
TightIdReco = "(isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0)"

generalMuMuMassMin2Candidates = generalDiMuCandidates.clone(
    mllCutMin = cms.double(2.0),
    muonSelection = cms.string("(p > 2.9 && abs(eta) < 2.4) && ("+SoftIdReco+" || "+TightIdReco+")"),
    candidateSelection = cms.string("mass > 2.0")
)
