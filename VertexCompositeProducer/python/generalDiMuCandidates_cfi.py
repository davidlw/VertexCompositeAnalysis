import FWCore.ParameterSet.Config as cms

generalDiMuCandidates = cms.EDProducer("DiMuProducer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),
    muonRecoAlgorithm = cms.InputTag('muons'),
    pfCandAlgorithm = cms.InputTag('particleFlow'),

    trackQualities = cms.vstring('highPurity'),
                                     
    tkChi2Cut = cms.double(9999.0), #trk Chi2 <
    tkNhitsCut = cms.int32(0), #trk Nhits >=
    tkPtCut = cms.double(0.0), #trk pT >
    tkEtaCut = cms.double(999.0), #trk abs(eta) <

    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.),
    dauLongImpactSigCut = cms.double(0.),

    mllCutMin = cms.double(2.0),
    mllCutMax = cms.double(5.0),

    #   PCA distance between tracks <
    tkDCACut = cms.double(9999.),
    vtxChi2Cut = cms.double(9999.0), #vtxChi2 <
    VtxChiProbCut = cms.double(0.0000001), #vtx prob >
    collinearityCut = cms.double(-2.0), #cos(pointAngle) >
    alphaCut = cms.double(999.0), #pointAngle <
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    DiMuMassCut = cms.double(1.4),
    dPtCut = cms.double(0.0),

    muonId = cms.string('TMOneStationTight'),

    isMuonId = cms.bool(False), 
    isPFMuon = cms.bool(False),
    isWrongSign = cms.bool(False)
)
