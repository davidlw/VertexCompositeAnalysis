import FWCore.ParameterSet.Config as cms

generalParticles = cms.EDProducer("ParticleProducer",

    pdgId = cms.int32(0),
    doSwap = cms.bool(False),
    vtxSortByTrkSize = cms.bool(True),

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

    # general settingss
    fitAlgo = cms.vuint32([0]),
    matchVertex = cms.bool(False),
    puMap = cms.vdouble(999., 999., 999., 999., 999.0, 4.0, 1.5, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0),

    # input collections
    primaryVertices = cms.InputTag('offlinePrimaryVertices'),
    tracks = cms.InputTag('generalTracks'),
    muons = cms.InputTag('patMuons'),
    electrons = cms.InputTag(''),
    taus = cms.InputTag(''),
    photons = cms.InputTag(''),
    pfParticles = cms.InputTag(''),
    jets = cms.InputTag(''),
    conversions = cms.InputTag(''),
    mva = cms.InputTag(''),
    dEdxInputs = cms.vstring(),
)
