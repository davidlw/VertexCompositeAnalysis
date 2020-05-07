import FWCore.ParameterSet.Config as cms

generalParticles = cms.EDProducer("ParticleProducer",

    pdgId = cms.int32(455),

    # particle selection
    preSelection = cms.string(""),
    pocaSelection = cms.string(""),
    postSelection = cms.string(""),
    finalSelection = cms.string(""),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(13), charge = cms.int32(-1), selection = cms.string("pt>1.0 && abs(eta)<2.4")),
        cms.PSet(pdgId = cms.int32(13), charge = cms.int32(+1), selection = cms.string("pt>1.0 && abs(eta)<2.4")),
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
    mva = cms.InputTag(''),
    dedxHarmonic2 = cms.InputTag(''),
)
