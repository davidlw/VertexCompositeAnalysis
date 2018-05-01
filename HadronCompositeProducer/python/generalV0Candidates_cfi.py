import FWCore.ParameterSet.Config as cms

generalV0Candidates = cms.EDProducer("V0Producer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),

    # These bools decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),
    selectD0s = cms.bool(False),
    selectDSs = cms.bool(False),
    selectDPMs = cms.bool(False),
    selectLambdaCToLamPis = cms.bool(False),
    selectLambdaCToKsPs = cms.bool(False),
    selectXis = cms.bool(True),
    selectOmegas = cms.bool(True),

    # Recommend leaving this one as is.
    vertexFitter = cms.InputTag('KalmanVertexFitter'),

    # set to true, uses tracks refit by the KVF for V0Candidate kinematics
    #  NOTE: useSmoothing is automatically set to FALSE
    #  if using the AdaptiveVertexFitter (which is NOT recommended)
    useSmoothing = cms.bool(True),
                                     
    # Select tracks using TrackBase::TrackQuality.
    # Select ALL tracks by leaving this vstring empty, which
    #   is equivalent to using 'loose'
    #trackQualities = cms.vstring('highPurity', 'goodIterative'),
    trackQualities = cms.vstring('loose'),
                                     
    # The next parameters are cut values.
    # All distances are in cm, all energies in GeV, as usual.

    # --Track quality/compatibility cuts--
    #   Normalized track Chi2 <
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track >=
    tkNhitsCut = cms.int32(6),
    #   Number of valid hits on track >=
    tkPtCut = cms.double(0.0),
    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(2.),
    dauLongImpactSigCut = cms.double(2.),
    batDauTransImpactSigCut = cms.double(2.),
    batDauLongImpactSigCut = cms.double(2.),
    # We calculate the PCA of the tracks quickly in RPhi, extrapolating
    # the z position as well, before vertexing.  Used in the following 2 cuts:
    #   m_pipi calculated at PCA of tracks <
    mPiPiCutMin = cms.double(0.0),
    mPiPiCutMax = cms.double(0.6),
    #   PCA distance between tracks <
    tkDCACut = cms.double(1.),

    # --V0 Vertex cuts--
    #   Vertex chi2 < 
    vtxChi2Cut = cms.double(7.0),
    #   Lambda collinearity cut
    #   (UNUSED)
    collinearityCut = cms.double(-2.0),
    #   Vertex radius cut >
    #   (UNUSED)
    rVtxCut = cms.double(0.0),
    #   V0 decay length from primary cut >
    #   (UNUSED)
    lVtxCut = cms.double(0.0),
    #   Radial vertex significance >
    vtxSignificance2DCut = cms.double(15.0),
    #   3D vertex significance using primary vertex
    #   (UNUSED)
    vtxSignificance3DCut = cms.double(0.0),
    
    # cuts for Xi 
    xiVtxChi2Cut = cms.double(7.0),
    xiCollinearityCut = cms.double(-2.0),
    xiRVtxCut = cms.double(0.0),
    xiLVtxCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(15.0),
    xiVtxSignificance3DCut = cms.double(0.0),

    #   V0 mass window, Candidate mass must be within these values of
    #     the PDG mass to be stored in the collection
    kShortMassCut = cms.double(0.07),
    lambdaMassCut = cms.double(0.05),
    d0MassCut = cms.double(0.2),
    dsMassCut = cms.double(0.2),
    dpmMassCut = cms.double(0.2),
    lambdaCMassCut = cms.double(0.3),
    xiMassCut = cms.double(0.10),
    omegaMassCut = cms.double(0.10),
    #   Mass window cut using normalized mass (mass / massError)
    #   (UNUSED)
    kShortNormalizedMassCut = cms.double(0.0),
    lambdaNormalizedMassCut = cms.double(0.0),
    # We check if either track has a hit inside (radially) the vertex position
    #  minus this number times the sigma of the vertex fit
    #  NOTE: Set this to -1 to disable this cut, which MUST be done
    #  if you want to run V0Producer on the AOD track collection!
    innerHitPosCut = cms.double(4.)
)
