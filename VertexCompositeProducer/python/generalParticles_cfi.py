import FWCore.ParameterSet.Config as cms

generalParticles = cms.EDProducer("ParticleProducer",

    pdgId = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(999999999.),
    fitAlgo = cms.vuint32([0]),

    # particle selection
    preSelection = cms.string(""),
    preMassSelection = cms.string(""),
    pocaSelection = cms.string(""),
    postSelection = cms.string(""),
    finalSelection = cms.string(""),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(0), charge = cms.int32(0), selection = cms.string("")),
    ]),

    # input collections
    primaryVertices = cms.InputTag('offlinePrimaryVertices'),
    tracks = cms.InputTag('generalTracks'),
    muons = cms.InputTag('patMuonsWithTrigger'),
    electrons = cms.InputTag(''),
    taus = cms.InputTag(''),
    photons = cms.InputTag(''),
    pfParticles = cms.InputTag(''),
    jets = cms.InputTag(''),
    conversions = cms.InputTag(''),
    mva = cms.InputTag(''),
    dedxHarmonic2 = cms.InputTag(''),
)
