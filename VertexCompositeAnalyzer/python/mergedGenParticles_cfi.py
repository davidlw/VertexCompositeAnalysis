import FWCore.ParameterSet.Config as cms

mergedGenParticles = cms.EDProducer("MergedGenParticleProducerRice",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
