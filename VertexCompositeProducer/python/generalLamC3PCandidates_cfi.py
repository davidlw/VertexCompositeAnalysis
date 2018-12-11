import FWCore.ParameterSet.Config as cms

generalLamC3PCandidates = cms.EDProducer("LamC3PProducer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),

    trackQualities = cms.vstring('highPurity'),
                                     
    tkChi2Cut = cms.double(7), #trk Chi2 <
    tkNhitsCut = cms.int32(5), #trk Nhits >=
    tkPtErrCut = cms.double(9999.0), #trk pT err <
    tkPtCut = cms.double(0.3), #trk pT >
    tkEtaCut = cms.double(999.0), #trk abs(eta) <
#    tkPtSumCut = cms.double(0.0), 
#    tkEtaDiffCut = cms.double(999.0), 

    mPiKPCutMin = cms.double(2.13),
    mPiKPCutMax = cms.double(2.45),
    mKPCutMin = cms.double(0.938+0.494),
    mKPCutMax = cms.double(2.45),

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
    lamCMassCut = cms.double(0.15),
    dPt3Cut = cms.double(1.0),

    isWrongSign = cms.bool(False),

# MVA 

    useAnyMVA = cms.bool(False),
    mvaType = cms.string('BDT'), 
    GBRForestLabel = cms.string('D0InpPb'),
    GBRForestFileName = cms.string('GBRForestfile.root'),
)
