import FWCore.ParameterSet.Config as cms

generalD0Candidates = cms.EDProducer("D0Producer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices4D'),

    trackQualities = cms.vstring('highPurity'),
                                     
    tkChi2Cut = cms.double(7), #trk Chi2 <
    tkNhitsCut = cms.int32(5), #trk Nhits >=
    tkPtErrCut = cms.double(9999.0), #trk pT err <
    tkPtMidCut = cms.double(0.3), #trk pT >
    tkPMidCut = cms.double(0.3), #trk p >
    tkPtFwdCut = cms.double(0.3), #trk pT >
    tkPFwdCut = cms.double(0.3), #trk p >
    tkEtaBound = cms.double(1.5), 
#    tkPtTrapCutPar0 = cms.double(0.75), #trk pT > a-b*(abs(eta)-c)
#    tkPtTrapCutPar1 = cms.double(0.28125), #trk pT > a-b*(abs(eta)-c), 0.3 at eta=3
#    tkPtTrapCutPar2 = cms.double(1.4), #trk pT > a-b*abs(eta)-c)
    tkEtaCut = cms.double(999.0), #trk abs(eta) <
    tkPtSumCut = cms.double(0.0), 
    tkEtaDiffCut = cms.double(999.0), 

    mPiKCutMin = cms.double(1.72),
    mPiKCutMax = cms.double(2.01),

    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.),
    dauLongImpactSigCut = cms.double(0.),

    #   PCA distance between tracks <
    tkDCACut = cms.double(9999.),
    vtxChi2Cut = cms.double(9999.0), #vtxChi2 <
    VtxChiProbCut = cms.double(0.0001), #vtx prob >
    collinearityCut2D = cms.double(-2.0), #cos(pointAngle) >
    collinearityCut3D = cms.double(-2.0), #cos(pointAngle) >
    alphaCut = cms.double(999.0), #pointAngle <
    alpha2DCut = cms.double(999.0), #pointAngle2D <
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    d0MassCut = cms.double(0.15),
    dPtCut = cms.double(0.0),

    isWrongSign = cms.bool(False),

# MVA 

    useAnyMVA = cms.bool(False),
    mvaType = cms.string('BDT'), 
    GBRForestLabel = cms.string('D0InpPb'),
    GBRForestFileName = cms.string('GBRForestfile.root'),
)
