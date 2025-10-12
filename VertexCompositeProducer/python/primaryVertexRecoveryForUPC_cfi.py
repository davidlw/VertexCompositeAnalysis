import FWCore.ParameterSet.Config as cms

primaryVertexRecoveryForUPC = cms.EDProducer('PrimaryVertexRecoveryForUPC',
  primaryVertexLabel = cms.InputTag('offlinePrimaryVertices'),
  beamSpotLabel = cms.InputTag('offlineBeamSpot'),
  TrackLabel = cms.InputTag('generalTracks'),
  TkFilterParameters = cms.VPSet(
    cms.PSet(
      algorithm = cms.string('filter'),
      maxD0Error = cms.double(10),
      maxD0Significance = cms.double(2),
      maxDzError = cms.double(10),
      maxEta = cms.double(2.4),
      maxNormalizedChi2 = cms.double(10),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(500),
      minPixelLayersWithHits = cms.int32(3),
      minPt = cms.double(1),
      minSiliconLayersWithHits = cms.int32(5),
      minValidStripHits = cms.int32(0),      
      trackQuality = cms.string('highPurity')
    ),
    cms.PSet(
      algorithm = cms.string('filter'),
      maxD0Error = cms.double(10),
      maxD0Significance = cms.double(2),
      maxDzError = cms.double(10),
      maxEta = cms.double(2.4),
      maxNormalizedChi2 = cms.double(10),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(10),
      minPixelLayersWithHits = cms.int32(3),
      minPt = cms.double(0.7),
      minSiliconLayersWithHits = cms.int32(5),
      minValidStripHits = cms.int32(0),      
      trackQuality = cms.string('highPurity')
    ),
    cms.PSet(
      algorithm = cms.string('filter'),
      maxD0Error = cms.double(1),
      maxD0Significance = cms.double(4),
      maxDzError = cms.double(1),
      maxEta = cms.double(2.4),
      maxNormalizedChi2 = cms.double(10),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(3),
      minPixelLayersWithHits = cms.int32(2),
      minPt = cms.double(0),
      minSiliconLayersWithHits = cms.int32(5),
      minValidStripHits = cms.int32(0),      
      trackQuality = cms.string('any')
    ),
    cms.PSet(
      algorithm = cms.string('filter'),
      maxD0Error = cms.double(10),
      maxD0Significance = cms.double(7),
      maxDzError = cms.double(10),
      maxEta = cms.double(2.5),
      maxNormalizedChi2 = cms.double(80),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(3),
      minPixelLayersWithHits = cms.int32(1),
      minPt = cms.double(0),
      minSiliconLayersWithHits = cms.int32(3),
      minValidStripHits = cms.int32(0),      
      trackQuality = cms.string('any')
    ),
    cms.PSet(
      algorithm = cms.string('none'),
      maxNtracks = cms.int32(10),
      minNtracks = cms.int32(2)
    )
  ),
  TkClusParameters = cms.VPSet(
    cms.PSet(
      algorithm = cms.string('gap'),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(3),
      zSeparation = cms.double(1)
    ),
    cms.PSet(
      Tmin = cms.double(2),
      Tpurge = cms.double(2),
      Tstop = cms.double(0.5),
      algorithm = cms.string('DA_vect'),
      block_size = cms.uint32(10000),
      convergence_mode = cms.int32(0),
      coolingFactor = cms.double(0.6),
      d0CutOff = cms.double(3),
      delta_highT = cms.double(0.01),
      delta_lowT = cms.double(0.001),
      dzCutOff = cms.double(3),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(3),
      overlap_frac = cms.double(0),
      runInBlocks = cms.bool(False),
      uniquetrkminp = cms.double(0),
      uniquetrkweight = cms.double(0.8),
      vertexSize = cms.double(0.006),
      zmerge = cms.double(0.01),
      zrange = cms.double(4)
    ),
    cms.PSet(
      Tmin = cms.double(4),
      Tpurge = cms.double(1),
      Tstop = cms.double(1),
      algorithm = cms.string('DA_vect'),
      block_size = cms.uint32(10000),
      convergence_mode = cms.int32(0),
      coolingFactor = cms.double(0.6),
      d0CutOff = cms.double(4),
      delta_highT = cms.double(0.01),
      delta_lowT = cms.double(0.001),
      dzCutOff = cms.double(5),
      maxNtracks = cms.int32(1000000000),
      minNtracks = cms.int32(3),
      overlap_frac = cms.double(0),
      runInBlocks = cms.bool(False),
      uniquetrkminp = cms.double(0),
      uniquetrkweight = cms.double(0.9),
      vertexSize = cms.double(0.01),
      zmerge = cms.double(0.02),
      zrange = cms.double(4)
    ),
    cms.PSet(
      algorithm = cms.string('none'),
      maxNtracks = cms.int32(4),
      minNtracks = cms.int32(2)
    )
  ),
  VtxFitParameters = cms.VPSet(
    cms.PSet(
      algorithm = cms.string('AdaptiveVertexFitter'),
      chi2cutoff = cms.double(2.5),
      maxDistanceToBeam = cms.double(1),
      minNclusters = cms.int32(2),
      minNdof = cms.double(0),
      useBeamConstraint = cms.bool(False)
    ),
    cms.PSet(
      algorithm = cms.string('AdaptiveVertexFitter'),
      chi2cutoff = cms.double(4),
      maxDistanceToBeam = cms.double(1),
      minNclusters = cms.int32(2),
      minNdof = cms.double(-1.1),
      useBeamConstraint = cms.bool(False)
    ),
    cms.PSet(
      algorithm = cms.string('AdaptiveVertexFitter'),
      chi2cutoff = cms.double(4),
      maxDistanceToBeam = cms.double(2),
      minNclusters = cms.int32(1),
      minNdof = cms.double(-1.1),
      useBeamConstraint = cms.bool(False)
    ),
    cms.PSet(
      algorithm = cms.string('KalmanVertexFitter'),
      maxDistanceToBeam = cms.double(2),
      minNclusters = cms.int32(1),
      minNdof = cms.double(0),
      useBeamConstraint = cms.bool(False)
    ),
    cms.PSet(
      algorithm = cms.string('AdaptiveVertexFitter'),
      chi2cutoff = cms.double(4),
      maxDistanceToBeam = cms.double(1),
      minNclusters = cms.int32(1),
      minNdof = cms.double(1),
      useBeamConstraint = cms.bool(True)
    )
  ),
  mightGet = cms.optional.untracked.vstring
)
