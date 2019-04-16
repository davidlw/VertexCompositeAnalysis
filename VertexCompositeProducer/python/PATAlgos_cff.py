import FWCore.ParameterSet.Config as cms

def doPATMuons(process, MC=False):
    # Check if sequence already ran
    if hasattr(process, 'patMuonsWithTrigger'): return

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")

    # Add trigger matching
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import switchOffAmbiguityResolution, addHLTL1Passthrough, useL1Stage2Candidates
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
    # Use L1 Stage2 (Run2 2016-onwards)
    addHLTL1Passthrough(process)
    useL1Stage2Candidates(process)
    if "hltL1extraParticles" in process.patTrigger.collections:
        process.patTrigger.collections.remove("hltL1extraParticles")
    if "hltGtStage2Digis:Muon" not in process.patTrigger.collections:
        process.patTrigger.collections.append("hltGtStage2Digis:Muon")
    process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGtStage2Digis:Muon")')
    process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
    process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
    process.muonMatchHLTL1.preselection = cms.string("")
    # Selection used to match online-offline muons
    process.muonL1Info.maxDeltaR = cms.double(0.3)
    process.muonL1Info.maxDeltaEta   = cms.double(0.2)
    process.muonL1Info.fallbackToME1 = cms.bool(True)
    process.muonMatchHLTL1.maxDeltaR = cms.double(0.3)
    process.muonMatchHLTL1.maxDeltaEta   = cms.double(0.2)
    process.muonMatchHLTL1.fallbackToME1 = cms.bool(True)
    process.muonMatchHLTL2.maxDeltaR = cms.double(0.3)
    process.muonMatchHLTL2.maxDPtRel = cms.double(10.0)
    process.muonMatchHLTL3.maxDeltaR = cms.double(0.1)
    process.muonMatchHLTL3.maxDPtRel = cms.double(10.0)

    # Add MC gen matching
    if MC:
        # Prune generated particles to muons and their parents
        process.genMuons = cms.EDProducer("GenParticlePruner",
            src = cms.InputTag("genParticles"),
            select = cms.vstring(
                "drop  *  ",                      # this is the default
                "++keep abs(pdgId) = 13"          # keep muons and their parents
            )
        )
        from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = cms.double(0.05)
        process.muonMatch.resolveByMatchQuality = cms.bool(True)
        process.muonMatch.matched = cms.string("genMuons")

    # Make a sequence
    if MC:
        process.patMuonSequence = cms.Sequence( process.genMuons + process.patMuonsWithTriggerSequence )
    else:
        process.patMuonSequence = cms.Sequence( process.patMuonsWithTriggerSequence )
