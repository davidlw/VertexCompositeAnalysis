import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANA',eras.Run2_2016_pA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:skim_mc.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_pA_v4')

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.dimucontana_mc = particleAna_mc.clone(
  recoParticles = cms.InputTag("generalMuMuCandidates"),
  genParticles = cms.untracked.InputTag("genParticles"),
  selectEvents = cms.string("")
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana_mc.root'))
process.p = cms.EndPath(process.dimucontana_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.p
)
