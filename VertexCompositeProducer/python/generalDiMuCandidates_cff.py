import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

HybridSoftIdReco2015 = "(isGlobalMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0)"
TightIdReco2015 = "(isGlobalMuon && globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0)"

generalMuMuMassMin7Candidates = generalDiMuCandidates.clone(
    mllCutMin = cms.double(7.0),
    muonSelection = cms.string("(p > 3.5 && abs(eta) < 2.4) && ("+HybridSoftIdReco2015+" || "+TightIdReco2015+")"),
    candidateSelection = cms.string("mass > 7.0")
)

generalMuMuMassMin0Candidates = generalDiMuCandidates.clone(
    mllCutMin = cms.double(0.0),
    muonSelection = cms.string("(p > 3.5 && abs(eta) < 2.4) && (isTrackerMuon || isGlobalMuon)"),
    candidateSelection = cms.string("mass > 0.0")
)
