import FWCore.ParameterSet.Config as cms

process = cms.Process("CALOTOPATADAPTER")

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

# Load the original CaloMuonMerger
process.load("RecoMuon.MuonIdentification.calomuons_cfi")

# Configure the RecoToPatMuonAdapter
process.recoToPatMuonAdapter = cms.EDProducer("RecoToPatMuonAdapter",
    src = cms.InputTag("calomuons"),  # Output from CaloMuonMerger
    preserveUserData = cms.bool(True),
    addQualityFlags = cms.bool(True)
)

# Output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('caloMuonMergerToPat_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_recoToPatMuonAdapter_*_*',  # Keep PAT muons from adapter
        'keep *_calomuons_*_*',  # Keep original reco muons for comparison
    )
)

# Path and EndPath definitions
process.p = cms.Path(
    process.calomuons *  # Original CaloMuonMerger
    process.recoToPatMuonAdapter  # Adapter to convert to PAT
)
process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.p, process.outpath)

# Print summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)