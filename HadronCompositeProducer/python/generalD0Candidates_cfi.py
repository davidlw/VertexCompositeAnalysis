import FWCore.ParameterSet.Config as cms

generalD0Candidates = cms.EDProducer("D0Producer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),

    trackQualities = cms.vstring('highPurity'),
                                     
    tkChi2Cut = cms.double(9999.0), #trk Chi2 <
    tkNhitsCut = cms.int32(0), #trk Nhits >=
    tkPtCut = cms.double(0.0), #trk pT >
    tkEtaCut = cms.double(999.0), #trk abs(eta) <

    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.),
    dauLongImpactSigCut = cms.double(0.),

    #   PCA distance between tracks <
    tkDCACut = cms.double(9999.),
    vtxChi2Cut = cms.double(9999.0), #vtxChi2 <
    VtxChiProbCut = cms.double(0.0001), #vtx prob >
    collinearityCut = cms.double(-2.0), #cos(pointAngle) >
    alphaCut = cms.double(999.0), #pointAngle <
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    d0MassCut = cms.double(0.2),
    dPtCut = cms.double(0.3),

    isWrongSign = cms.bool(False)
)
