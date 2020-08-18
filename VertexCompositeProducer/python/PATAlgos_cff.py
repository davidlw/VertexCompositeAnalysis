import FWCore.ParameterSet.Config as cms

def doPATMuons(process):
    # Check if sequence already ran
    if hasattr(process, 'patMuonsWithTrigger'): return

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
    process.patMuonsWithoutTrigger.embedTrack = cms.bool(True)
    process.patMuonsWithoutTrigger.embedMuonBestTrack = cms.bool(True)
    process.patMuonsWithoutTrigger.forceBestTrackEmbedding = cms.bool(False)
    process.patMuonsWithoutTrigger.embedTunePMuonBestTrack = cms.bool(False)
    process.patMuonsWithoutTrigger.embedCombinedMuon = cms.bool(True)
    process.patMuonsWithoutTrigger.embedStandAloneMuon = cms.bool(True)
    process.patMuonsWithoutTrigger.embedPFCandidate = cms.bool(False)
    process.patMuonsWithoutTrigger.embedCaloMETMuonCorrs = cms.bool(False)
    process.patMuonsWithoutTrigger.embedTcMETMuonCorrs = cms.bool(False)
    process.patMuonsWithoutTrigger.embedPfEcalEnergy = cms.bool(False)
    process.patMuonsWithoutTrigger.embedPickyMuon = cms.bool(False)
    process.patMuonsWithoutTrigger.embedTpfmsMuon = cms.bool(False)
    process.patMuonsWithoutTrigger.embedDytMuon = cms.bool(False)
    process.patMuonsWithoutTrigger.embedHighLevelSelection = cms.bool(False)

    # Add trigger matching
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import switchOffAmbiguityResolution, addHLTL1Passthrough, useL1Stage2Candidates
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
    # Use L1 Stage2 (Run2 2016-onwards)
    addHLTL1Passthrough(process)
    useL1Stage2Candidates(process)
    if "hltL1extraParticles" in process.patTrigger.collections:
        process.patTrigger.collections.remove("hltL1extraParticles")
    if "hltGmtStage2Digis:Muon" not in process.patTrigger.collections:
        process.patTrigger.collections.append("hltGmtStage2Digis:Muon")
    process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGmtStage2Digis:Muon")')
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

    # Make a sequence
    process.patMuonSequence = cms.Sequence( process.patMuonsWithTriggerSequence )

    ### Temporal fix for the PAT Trigger prescale warnings.
    process.patTriggerFull.l1GtReadoutRecordInputTag = cms.InputTag("gtDigis","","RECO")
    process.patTriggerFull.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
    process.patTriggerFull.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
