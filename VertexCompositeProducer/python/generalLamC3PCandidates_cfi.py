import FWCore.ParameterSet.Config as cms

generalLamC3PCandidates = cms.EDProducer("LamC3PProducer",
                                     
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
#    tkPtTrapCutPar0 = cms.double(0.3), #trk pT > a-b*(abs(eta)-c)
#    tkPtTrapCutPar1 = cms.double(0.28125), #trk pT > a-b*(abs(eta)-c)
#    tkPtTrapCutPar2 = cms.double(1.4), #trk pT > a-b*abs(eta)-c)
    tkEtaCut = cms.double(999.0), #trk abs(eta) <
#    tkPtSumCut = cms.double(0.0), 
#    tkEtaDiffCut = cms.double(999.0), 

    mPiKPCutMin = cms.double(2.13),
    mPiKPCutMax = cms.double(2.45),
    mKPCutMin = cms.double(0.938+0.140),
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

    dPt3CutMin = cms.double(0.0),
    dPt3CutMax = cms.double(10000.0),
    dY3CutMin = cms.double(-10000.0),
    dY3CutMax = cms.double(10000.0),

    isWrongSign = cms.bool(False),

# MVA 

    useAnyMVA = cms.bool(False),
    mvaType = cms.string('BDT'), 
    GBRForestLabel = cms.string('D0InpPb'),
    GBRForestFileName = cms.string('GBRForestfile.root'),
)
