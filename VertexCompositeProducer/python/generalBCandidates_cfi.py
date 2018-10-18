import FWCore.ParameterSet.Config as cms

generalBCandidates = cms.EDProducer("BProducer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),
    d0RecoAlgorithm = cms.InputTag('generalD0CandidatesNew','D0'),

    trackQualities = cms.vstring('highPurity'),
                                     
    batTkChi2Cut = cms.double(7), #trk Chi2 <
    batTkNhitsCut = cms.int32(5), #trk Nhits >=
    batTkPtErrCut = cms.double(9999.0), #trk pT err <
    batTkPtCut = cms.double(0.3), #trk pT >
    batTkEtaCut = cms.double(999.0), #trk abs(eta) <
    batTkPtSumCut = cms.double(0.0), 
    batTkEtaDiffCut = cms.double(999.0), 

    mPiDCutMin = cms.double(5.279-0.35),
    mPiDCutMax = cms.double(5.279+0.35),

    #   Track impact parameter significance >
    batDauTransImpactSigCut = cms.double(0.),
    batDauLongImpactSigCut = cms.double(0.),

    #   PCA distance between tracks <
    batTkDCACut = cms.double(9999.),
    bVtxChi2Cut = cms.double(9999.0), #vtxChi2 <
    bVtxChiProbCut = cms.double(0.0001), #vtx prob >
    bCollinCut2D = cms.double(-2.0), #cos(pointAngle) >
    bCollinCut3D = cms.double(-2.0), #cos(pointAngle) >
    bAlphaCut = cms.double(999.0), #pointAngle <
    bAlpha2DCut = cms.double(999.0), #pointAngle2D <
    bVtxSignificance2DCut = cms.double(0.0),
    bVtxSignificance3DCut = cms.double(0.0),
    bVtx2DCut = cms.double(0.0),
    bVtx3DCut = cms.double(0.0),
    bMassCut = cms.double(0.3),
    bPtCut = cms.double(0.0),

    isWrongSignB = cms.bool(False)

# MVA 
#    useAnyMVA = cms.bool(False),
#    mvaType = cms.string('BDT'), 
#    GBRForestLabel = cms.string('D0InpPb'),
#    GBRForestFileName = cms.string('GBRForestfile.root'),
)
