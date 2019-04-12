import FWCore.ParameterSet.Config as cms

# Coincidence of HF towers above threshold
from VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff import *

#Reject events with pile-up, optionally
from VertexCompositeAnalysis.VertexCompositeProducer.pileUpFilter_cff import *

# Selection of at least a two-track fitted vertex
primaryVertexFilterPA = cms.EDFilter("VertexSelector",
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

collisionEventSelectionPA = cms.Sequence(hfCoincFilter *
                                         primaryVertexFilterPA *
                                         NoScraping
                                         )

collisionEventSelectionPA_rejectPU = cms.Sequence(collisionEventSelectionPA *
                                                  olvFilter_pPb8TeV_dz1p0)

