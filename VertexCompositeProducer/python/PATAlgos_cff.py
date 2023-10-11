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
        embedHighLevelSelection = cms.bool(True),
        # gen information
        addGenMatch = cms.bool(False),
        embedGenMatch = cms.bool(False),
        # PF information
        useParticleFlow = cms.bool(False),
        embedPFCandidate = cms.bool(False),
        embedPfEcalEnergy = cms.bool(False),
        # extra information
        addInverseBeta = cms.bool(False),
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


def doPATElectrons(process):
    # Check if sequence already ran
    if hasattr(process, 'patElectrons'): return

    # Make PAT Electrons
    from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import patElectrons
    process.patElectrons = patElectrons.clone(
        electronSource = cms.InputTag('lowPtGsfElectrons'),
        # track information
        embedTrack                  = cms.bool(True),
        embedGsfElectronCore        = cms.bool(True),
        embedGsfTrack               = cms.bool(True),
        embedSuperCluster           = cms.bool(True),
        embedSeedCluster            = cms.bool(True),
        embedBasicClusters          = cms.bool(True),
        embedPreshowerClusters      = cms.bool(False),
        embedRecHits                = cms.bool(False),
        # high level information
        embedHighLevelSelection = cms.bool(False),
        # gen information
        addGenMatch = cms.bool(False),
        embedGenMatch = cms.bool(False),
        # PF information
        embedPflowSuperCluster      = cms.bool(False),
        embedPflowBasicClusters     = cms.bool(False),
        embedPflowPreshowerClusters = cms.bool(False),
        embedPFCandidate            = cms.bool(False),
        # extra information
        addMVAVariables         = cms.bool(False),
        isoDeposits             = cms.PSet(),
        isolationValues         = cms.PSet(),
        isolationValuesNoPFId   = cms.PSet(),
    )

    patElectrons.userData.userInts.src    = []
    patElectrons.userData.userFloats.src  = []
    patElectrons.userData.userCands.src   = []
    patElectrons.userData.userClasses.src = []

    # Make a sequence
    process.patElectronSequence = cms.Sequence( process.patElectrons )


def changeToMiniAOD(process):

    process.patTriggerFull = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
        patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
        triggerResults              = cms.InputTag('TriggerResults::HLT'),
        unpackFilterLabels          = cms.bool(True)
    )
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.eventFilter_HM.insert(0, process.unpackedTracksAndVertices)

    if hasattr(process, "patMuons"):
        process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
        process.patMuons = process.unpackedMuons.clone()

    from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
    process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
    process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
