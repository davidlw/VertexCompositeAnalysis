import FWCore.ParameterSet.Config as cms

process = cms.Process("PATCALOMUONMERGER")

# Load the standard sequences
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Add your input files here
        '/store/mc/RunIISummer16MiniAODv2/...',
    )
)

# Maximum number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Load the PAT CaloMuon merger module
process.load("YourPackage.YourSubsystem.patCaloMuonMerger_cfi")

# Customize the module if needed
process.patCaloMuonMerger.minCaloCompatibility = cms.double(0.7)
process.patCaloMuonMerger.deltaR = cms.double(0.05)

# Output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patCaloMuonMerged_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_patCaloMuonMerger_*_*',
        'keep *_slimmedMuons_*_*',  # Keep original PAT muons for comparison
    )
)

# Path and EndPath definitions
process.p = cms.Path(process.patCaloMuonMerger)
process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.p, process.outpath)