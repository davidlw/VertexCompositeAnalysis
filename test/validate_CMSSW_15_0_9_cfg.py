import FWCore.ParameterSet.Config as cms

process = cms.Process("VALIDATECMSSW1509")

# Load standard sequences for CMSSW_15_0_9
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Run 3 Global tag for CMSSW_15_0_9
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# Input source - using Run 3 miniAOD format
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Example Run 3 MC file - replace with available file
        'file:/path/to/Run3/miniAOD.root',
        # For testing, you can try:
        # '/store/mc/Run3Summer23MiniAODv4/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2023_realistic_v14-v2/2810000/...'
    )
)

# Limit events for validation
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Enhanced message logging for validation
process.MessageLogger.cerr.FwkReport.reportEvery = 5
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')

# Load both CaloMuon merger modules
process.miniAODCaloMuonMerger = cms.EDProducer("MiniAODCaloMuonMerger",
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    minPt = cms.double(2.0),
    maxEta = cms.double(2.4),
    addCaloMuonsFromPacked = cms.bool(True),
    enhanceExistingMuons = cms.bool(True),
    requireTrackerTrack = cms.bool(False)
)

# Validation analyzer to check output
process.validationAnalyzer = cms.EDAnalyzer("CandViewCountAnalyzer",
    src = cms.InputTag("miniAODCaloMuonMerger"),
    verbose = cms.untracked.bool(True)
)

# Additional validation analyzer for detailed checks
process.muonValidator = cms.EDAnalyzer("CandViewHistoAnalyzer",
    src = cms.InputTag("miniAODCaloMuonMerger"),
    histograms = cms.VPSet(
        cms.PSet(
            min = cms.untracked.double(0.0),
            max = cms.untracked.double(200.0),
            nbins = cms.untracked.int32(100),
            name = cms.untracked.string("pt"),
            description = cms.untracked.string("Muon pT [GeV]"),
            plotquantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            min = cms.untracked.double(-3.0),
            max = cms.untracked.double(3.0),
            nbins = cms.untracked.int32(60),
            name = cms.untracked.string("eta"),
            description = cms.untracked.string("Muon eta"),
            plotquantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            min = cms.untracked.double(-3.2),
            max = cms.untracked.double(3.2),
            nbins = cms.untracked.int32(64),
            name = cms.untracked.string("phi"),
            description = cms.untracked.string("Muon phi"),
            plotquantity = cms.untracked.string("phi")
        )
    )
)

# Compare with original slimmed muons
process.originalMuonValidator = cms.EDAnalyzer("CandViewCountAnalyzer",
    src = cms.InputTag("slimmedMuons"),
    verbose = cms.untracked.bool(True)
)

# TFileService for validation histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("validation_CMSSW_15_0_9.root")
)

# Output module for detailed inspection
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('validation_CMSSW_15_0_9_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_miniAODCaloMuonMerger_*_*',
        'keep *_slimmedMuons_*_*',
        'keep *_packedPFCandidates_*_*',
        # Keep some additional collections for debugging
        'keep recoVertexs_offlineSlimmedPrimaryVertices_*_*',
        'keep *_slimmedMETs_*_*'
    )
)

# Validation path
process.validation_path = cms.Path(
    process.originalMuonValidator *
    process.miniAODCaloMuonMerger *
    process.validationAnalyzer *
    process.muonValidator
)

process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.validation_path, process.outpath)

# CMSSW_15_0_9 specific options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)

# Print configuration summary
print("=== CMSSW_15_0_9 Validation Configuration ===")
print("Global Tag:", process.GlobalTag.globaltag.value())
print("Input events:", process.maxEvents.input.value())
print("Threading: {} threads, {} streams".format(
    process.options.numberOfThreads.value(),
    process.options.numberOfStreams.value()
))
print("===========================================")

# Additional debugging for validation
process.MessageLogger.categories.append('MiniAODCaloMuonMerger')
process.MessageLogger.cerr.MiniAODCaloMuonMerger = cms.untracked.PSet(
    limit = cms.untracked.int32(10)
)