import FWCore.ParameterSet.Config as cms

# Coincidence of HF towers above threshold
from VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff import *

# Cluster-shape filter
from VertexCompositeAnalysis.VertexCompositeProducer.HIClusterCompatibilityFilter_cfi import *

# Selection of at least a two-track fitted vertex
primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("hiSelectedVertex"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

collisionEventSelection = cms.Sequence(hfCoincFilter3 *
                                       primaryVertexFilter *
                                       clusterCompatibilityFilter
                                      )
