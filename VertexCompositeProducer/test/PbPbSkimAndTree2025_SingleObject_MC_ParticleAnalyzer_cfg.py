import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2025_OXY)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it//store/user/srdas/Hijing_NoPU_1kEvents_OO_5360GeV_GenSim_030925/Hijing_NoPU_1kEvents_OO_5360GeV_AOD_030925/250310_031931/0000/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root")
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_mcRun3_2024_realistic_HI_v13')

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

process.muons = generalParticles.clone(
    pdgId = cms.uint32(13),
    muons = cms.InputTag('patMuons'),
)

process.electrons = generalParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patElectrons')
)

process.lowPtElectrons = generalParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patLowPtElectrons')
)

process.photons = generalParticles.clone(
    pdgId = cms.uint32(22),
    photons = cms.InputTag('patPhotons')
)

process.convertedPhotons = generalParticles.clone(
    pdgId = cms.uint32(22),
    conversions = cms.InputTag('allConversions')
)

process.tracks = generalParticles.clone(
    tracks = cms.InputTag('generalTracks'),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)

process.pixelTracks = generalParticles.clone(
    tracks = cms.InputTag('hiConformalPixelTracks')
)

process.pfCandidates = generalParticles.clone(
    pfParticles = cms.InputTag('particleFlow'),
    tracks = cms.InputTag('')
)

# Add PAT objects
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons, doPATElectrons, doPATPhotons
doPATMuons(process)
doPATElectrons(process)
doPATPhotons(process)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility * process.primaryVertexFilter)

# Define the analysis steps
process.muon_step = cms.Path(process.patMuonSequence * process.muons)
process.electron_step = cms.Path(process.patElectronSequence * process.electrons * process.lowPtElectrons)
process.photon_step = cms.Path(process.patPhotonSequence * process.photons * process.convertedPhotons)
process.track_step = cms.Path(process.tracks  * process.pixelTracks * process.pfCandidates)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc

process.muonAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("muons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.elecAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("electrons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.lowPtElecAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("lowPtElectrons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.phoAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("photons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.convAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("convertedPhotons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.trackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("tracks"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.pixelTrackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("pixelTracks"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.pfAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("pfCandidates"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('obj_ana_mc.root'))
process.p = cms.EndPath(process.muonAna * process.elecAna * process.lowPtElecAna * process.phoAna * process.convAna * process.trackAna * process.pixelTrackAna * process.pfAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.muon_step,
    process.electron_step,
    process.photon_step,
    process.track_step,
    process.p
)
