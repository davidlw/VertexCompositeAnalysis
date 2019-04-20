import FWCore.ParameterSet.Config as cms

analysisSkimContent = cms.PSet(
        outputCommands = cms.untracked.vstring(
        'drop *',
        # event
        'keep *_offlineBeamSpot_*_*',
        'keep *_TriggerResults_*_*',

        # centrality
        'keep *_hiCentrality_*_*',
        'keep *_hiEvtPlane_*_*',
        'keep *_hiEvtPlaneFlat_*_*',
        'keep *_centralityBin_*_*',

        # V0
        'keep patCompositeCandidates_*_*_*',

        # L1 Trigger Prescales
        'keep *_gtDigis__*',

        # vertex
        'keep *_hiSelectedVertex_*_*',

        # mc (if present)
        'keep recoGenParticles_genParticles_*_*',
        'keep *_genMuons_*_*',
        'keep *_generator_*_*',
        )
)
