import FWCore.ParameterSet.Config as cms

#Reject events with pile-up, optionally
from VertexCompositeAnalysis.VertexCompositeProducer.pileUpFilter_cff import *

# Selection of at least a two-track fitted vertex
primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

#Reject beam scraping events standard pp and pA configuration
NoScraping = cms.EDFilter("FilterOutScraping",
 applyfilter = cms.untracked.bool(True),
 debugOn = cms.untracked.bool(False),
 numtrack = cms.untracked.uint32(10),
 thresh = cms.untracked.double(0.25)
)

collisionEventSelection = cms.Sequence(primaryVertexFilter *
                                       NoScraping
                                      )

collisionEventSelection_rejectPU = cms.Sequence(collisionEventSelection *
                                                olvFilter_pPb8TeV_dz1p0)

