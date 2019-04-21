import FWCore.ParameterSet.Config as cms

generalDiMuCandidates = cms.EDProducer("DiMuProducer",

    # InputTag that tells which TrackCollection to use for vertexing
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),
    muonRecoAlgorithm = cms.InputTag('patMuonsWithTrigger'),

    # Muon selection
    muonSelection = cms.string("isGlobalMuon || muonID('TMOneStationTight')"),
    tkdXYCut = cms.double(99999.), #|dXY| <
    tkdZCut = cms.double(99999.), #|dZ| <

    # Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.0),
    dauLongImpactSigCut = cms.double(0.0),

    mllCutMin = cms.double(0.0),
    mllCutMax = cms.double(10000.0),

    # PCA distance between tracks <
    tkDCACut = cms.double(99999.),
    vtxChi2Cut = cms.double(99999.), #vtxChi2 <
    VtxChiProbCut = cms.double(0.0), #vtx prob >
    collinearityCut = cms.double(-2.0), #cos(pointAngle) >
    alphaCut = cms.double(99999.), #pointAngle <
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    candidateSelection = cms.string(""),

    isWrongSign = cms.bool(False)
)
