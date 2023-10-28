import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023)

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
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/ChiC_SUPERCHIC_5p36TeV_2023Run3_HIFwd_RECO_AOD_2023_10_23/231022_221941/0000/SUPERCHIC_chiC_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/a/anstahll/work/Run3_2023/ERA/CMSSW_13_2_6_patch2/src/STARLIGHT/STARLIGHT_QED_elec_RECO_DEBUG.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2023_10_23/STARLIGHT_5p36TeV_2023Run3/double_diff_STARLIGHT_5p36TeV_2023Run3_pp_RECO_AOD_2023_10_23/231022_214327/0000/STARLIGHT_double_diff_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2023_10_26/STARLIGHT_5p36TeV_2023Run3/double_diff_STARLIGHT_5p36TeV_2023Run3_pp_RECO_AOD_2023_10_26/231026_194253/0000/STARLIGHT_double_diff_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/LbyL_SUPERCHIC_5p36TeV_2023Run3_HIFwd_RECO_AOD_2023_10_23/231022_221920/0000/SUPERCHIC_LbyL_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/ChiC_SUPERCHIC_5p36TeV_2023Run3_HIFwd_RECO_AOD_2023_10_23/231022_221941/0000/SUPERCHIC_chiC_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/ChiC_SUPERCHIC_5p36TeV_2023Run3_HIFwd_RECO_MOD_AOD_2023_10_23/231022_222620/0000/SUPERCHIC_chiC_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/ChiC_SUPERCHIC_5p36TeV_2023Run3_pp_RECO_AOD_2023_10_23/231022_215005/0000/SUPERCHIC_chiC_RECO_1.root"),
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/a/anstahll/work/Run3_2023/ERA/CMSSW_13_2_6_patch2/src/SUPERCHIC/SUPERCHIC_chiC_RECO.root"),
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2023_10_23/STARLIGHT_5p36TeV_2023Run3/QED_muon_STARLIGHT_5p36TeV_2023Run3_pp_RECO_AOD_2023_10_23/231022_214310/0000/STARLIGHT_QED_muon_RECO_1.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_mcRun3_2023_realistic_HI_v5')

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

process.muons = generalParticles.clone(
    pdgId = cms.int32(13),
    muons = cms.InputTag('patMuons')
)

process.electrons = generalParticles.clone(
    pdgId = cms.int32(11),
    electrons = cms.InputTag('patElectrons')
)

process.lowPtElectrons = generalParticles.clone(
    pdgId = cms.int32(11),
    electrons = cms.InputTag('patLowPtElectrons')
)

process.photons = generalParticles.clone(
    pdgId = cms.int32(22),
    photons = cms.InputTag('patPhotons')
)

process.convertedPhotons = generalParticles.clone(
    pdgId = cms.int32(-22),
    conversions = cms.InputTag('allConversions')
)

process.tracks = generalParticles.clone(
    tracks = cms.InputTag('generalTracks')
)

process.pixelTracks = generalParticles.clone(
    tracks = cms.InputTag('hiConformalPixelTracks')
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
process.track_step = cms.Path(process.tracks * process.pixelTracks)

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
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.lowPtElecAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("lowPtElectrons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(10.0),
)

process.phoAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("photons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(10.0),
)

process.convAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("convertedPhotons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.trackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("tracks"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.pixelTrackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("pixelTracks"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('obj_ana_mc.root'))
process.p = cms.EndPath(process.muonAna * process.elecAna * process.lowPtElecAna * process.phoAna * process.convAna * process.trackAna * process.pixelTrackAna)

# Add fix for SuperChic
process.genParticles = cms.EDProducer('GenParticleSuperChicFixer', genPars = cms.InputTag('genParticles::SIM'))
process.gen_step = cms.Path(process.genParticles)

# Define the process schedule
process.schedule = cms.Schedule(
    process.gen_step,
    process.muon_step,
    process.electron_step,
    process.photon_step,
    process.track_step,
    process.p
)
