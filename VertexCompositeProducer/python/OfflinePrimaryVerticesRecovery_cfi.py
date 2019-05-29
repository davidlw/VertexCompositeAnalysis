import FWCore.ParameterSet.Config as cms

offlinePrimaryVerticesRecovery = cms.EDProducer(
    "PrimaryVertexRecoveryProducer",

    redoAllVertices = cms.bool(False),#Only turn this on if you want to spend lots of CPU time, nt recommended
    oldVertexLabel = cms.InputTag("offlinePrimaryVertices"), 

    verbose = cms.untracked.bool(False),
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    
    TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(999.0),
        minPixelLayersWithHits=cms.int32(0),
        minSiliconLayersWithHits = cms.int32(0),
        maxD0Significance = cms.double(999.0), 
        minPt = cms.double(0.0),
        maxEta = cms.double(999.0),
        trackQuality = cms.string("any")
    ),

    TkClusParameters = cms.PSet(
        algorithm = cms.string("gap"),
        TkGapClusParameters = cms.PSet(
            zSeparation = cms.double(1.0)        
        )
    ),

    vertexCollections = cms.VPSet(
     [cms.PSet(label=cms.string(""),
               algorithm=cms.string("AdaptiveVertexFitter"),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(0.0),
               useBeamConstraint = cms.bool(False),
               maxDistanceToBeam = cms.double(1.0)
               ),
      #do not run vertexing with beamspot constraint
      #cms.PSet(label=cms.string("WithBS"),
      #         algorithm = cms.string('AdaptiveVertexFitter'),
      #         chi2cutoff = cms.double(2.5),
      #         minNdof=cms.double(2.0),
      #         useBeamConstraint = cms.bool(True),
      #         maxDistanceToBeam = cms.double(1.0)
      #         )
      ]
    )
)
