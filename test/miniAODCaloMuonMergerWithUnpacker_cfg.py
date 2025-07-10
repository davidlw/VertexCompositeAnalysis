import FWCore.ParameterSet.Config as cms

process = cms.Process("MINIAODCALOMUONWITHUNPACKER")

# Load the standard sequences
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Global tag - adjust for your data/MC and era
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Input source - replace with your miniAOD files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Example miniAOD file - replace with your actual files
        'file:/path/to/your/miniAOD.root',
    )
)

# Maximum number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Message logger configuration
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Load your existing unpacker modules
# Note: Adjust these paths according to your actual package structure
# You may need to modify these imports based on where you have the unpackers

# Track and Vertex Unpacker
process.unpackedTracksAndVertices = cms.EDProducer("TrackAndVertexUnpacker",
    # Configure according to your TrackAndVertexUnpacker.cc
    packedCandidates = cms.InputTag("packedPFCandidates"),
    # Add other parameters from your implementation
)

# Muon Unpacker
process.muonUnpacker = cms.EDProducer("MuonUnpacker",
    # Configure according to your MuonUnpacker.cc
    slimmedMuons = cms.InputTag("slimmedMuons"),
    # Add other parameters from your implementation
)

# CaloMuon Merger with Unpacker integration
process.miniAODCaloMuonMergerWithUnpacker = cms.EDProducer("MiniAODCaloMuonMergerWithUnpacker",
    unpackedMuons = cms.InputTag("muonUnpacker"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    unpackedTracks = cms.InputTag("unpackedTracksAndVertices"),
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    minPt = cms.double(2.0),
    recalculateCaloCompatibility = cms.bool(True),
    addMissingCaloMuons = cms.bool(True)
)

# Output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('miniAODCaloMuonMergerWithUnpacker_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_miniAODCaloMuonMergerWithUnpacker_*_*',
        'keep *_muonUnpacker_*_*',  # Keep unpacked muons for comparison
        'keep *_slimmedMuons_*_*',  # Keep original slimmed muons
    )
)

# Path and EndPath definitions
process.p = cms.Path(
    process.unpackedTracksAndVertices *
    process.muonUnpacker *
    process.miniAODCaloMuonMergerWithUnpacker
)

process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.p, process.outpath)

# Print summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(False)
)

# Optional: Add some debugging output
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)
process.MessageLogger.cerr.MiniAODCaloMuonMergerWithUnpacker = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)