import FWCore.ParameterSet.Config as cms

def doPATMuons(process):
    # Check if sequence already ran
    if hasattr(process, 'patMuons'): return

    # Make PAT Muons
    from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons
    process.patMuons = patMuons.clone(
        muonSource = cms.InputTag('muons'),
        # track information
        embedTrack = cms.bool(True),
        embedCombinedMuon = cms.bool(True),
        embedStandAloneMuon = cms.bool(True),
        embedMuonBestTrack = cms.bool(True),
        forceBestTrackEmbedding = cms.bool(False),
        embedTunePMuonBestTrack = cms.bool(False),
        embedCaloMETMuonCorrs = cms.bool(False),
        embedTcMETMuonCorrs = cms.bool(False),
        embedPickyMuon = cms.bool(False),
        embedTpfmsMuon = cms.bool(False),
        embedDytMuon = cms.bool(False),
        # high level information
        embedHighLevelSelection = cms.bool(False),
        # gen information
        addGenMatch = cms.bool(False),
        embedGenMatch = cms.bool(False),
        # PF information
        useParticleFlow = cms.bool(False),
        embedPFCandidate = cms.bool(False),
        embedPfEcalEnergy = cms.bool(False),
        # extra information
        addEfficiencies = cms.bool(False),
        addResolutions = cms.bool(False),
        userIsolation = cms.PSet(),
        isoDeposits = cms.PSet(),
        isolationValues = cms.PSet(),
    )

    patMuons.userData.userInts.src    = []
    patMuons.userData.userFloats.src  = []
    patMuons.userData.userCands.src   = []
    patMuons.userData.userClasses.src = []

    # Make a sequence
    process.patMuonSequence = cms.Sequence( process.patMuons )

def changeToMiniAOD(process):

    process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
        patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
        triggerResults              = cms.InputTag('TriggerResults::HLT'),
        unpackFilterLabels          = cms.bool(True)
    )
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.eventFilter.insert(0, process.unpackedTracksAndVertices)

    if hasattr(process, "patMuons"):
        process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
        process.patMuons = process.unpackedMuons.clone()

    from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
    process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")



def changeToMiniAODMC(process):

    process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
        patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
        triggerResults              = cms.InputTag('TriggerResults::HLT'),
        unpackFilterLabels          = cms.bool(True)
    )
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    #process.load('GeneratorInterface.RivetInterface.mergedGenParticles_cfi')
    process.load('VertexCompositeAnalysis.VertexCompositeAnalyzer.mergedGenParticles_cfi')
    eventFilters = [ process.unpackedTracksAndVertices , process.mergedGenParticles ]
    for P in eventFilters:
        process.eventFilter.insert(0, P)

    if hasattr(process, "patMuons"):
        process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
        process.patMuons = process.unpackedMuons.clone()

    from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
    process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"genParticles","mergedGenParticles")



def changeToMiniAODJet(process):

    process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
        patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
        triggerResults              = cms.InputTag('TriggerResults::HLT'),
        unpackFilterLabels          = cms.bool(True)
    )

    process.packedJetCSPFCandidates = cms.EDProducer("PATJetConstituentSelector",
        src = cms.InputTag('slimmedJets'),
        cut = cms.string('') #pt>500 && abs(eta)<2.0')
    )

    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.unpackedTracksAndVertices.packedCandidates = cms.InputTag("packedJetCSPFCandidates")

    eventFilters = [ process.unpackedTracksAndVertices , process.packedJetCSPFCandidates ]
#    eventFilters = [ process.packedJetCSPFCandidates ]
    for P in eventFilters:
        process.eventFilter.insert(0, P)

    from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
    process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")

def changeToMiniAODJetMC(process):

    process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
        patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
        triggerResults              = cms.InputTag('TriggerResults::HLT'),
        unpackFilterLabels          = cms.bool(True)
    )

    process.unpackedGenJetPFCandidate = cms.EDProducer("PATGenJetConstituentSelector",
        src = cms.InputTag('slimmedGenJets'),
        cut = cms.string('pt>500 && abs(eta)<1.6')
    )

    process.packedJetCSPFCandidates = cms.EDProducer("PATJetConstituentSelector",
        src = cms.InputTag('slimmedJets'),
        cut = cms.string('pt>500 && abs(eta)<1.6')
    )

    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.unpackedTracksAndVertices.packedCandidates = cms.InputTag("packedJetCSPFCandidates")

    eventFilters = [ process.unpackedTracksAndVertices, process.packedJetCSPFCandidates, process.unpackedGenJetPFCandidate]
    for P in eventFilters:
        process.eventFilter.insert(0, P)

    from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
    process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"genParticles","unpackedGenJetPFCandidate")
