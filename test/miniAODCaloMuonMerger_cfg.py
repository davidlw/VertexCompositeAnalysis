import FWCore.ParameterSet.Config as cms

process = cms.Process("MINIAODCALOMUON")

# Load the standard sequences
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Global tag - adjust for your data/MC and era
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# Input source - replace with your miniAOD files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Example miniAOD file - replace with your actual files
        'file:/path/to/your/miniAOD.root',
        # For testing, you can also use:
        # '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrackingTest_80X_mcRun2_asymptotic_2016_TrackingTest-v6/50000/...'
    )
)

# Maximum number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Message logger configuration
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Option 1: Use the direct miniAOD approach
process.load("YourPackage.YourSubsystem.miniAODCaloMuonMerger_cfi")

# Customize if needed
process.miniAODCaloMuonMerger.minCaloCompatibility = cms.double(0.7)
process.miniAODCaloMuonMerger.addCaloMuonsFromPacked = cms.bool(True)
process.miniAODCaloMuonMerger.enhanceExistingMuons = cms.bool(True)

# Option 2: Use with your existing unpackers (comment out Option 1 and uncomment this)
# process.load("YourPackage.YourSubsystem.miniAODCaloMuonWorkflow_cff")

# EDAnalyzer to check the output (optional)
process.muonAnalyzer = cms.EDAnalyzer("CandViewHistoAnalyzer",
    src = cms.InputTag("miniAODCaloMuonMerger"),
    # or cms.InputTag("miniAODCaloMuonMergerWithUnpacker") if using the unpacker version
    histograms = cms.VPSet(
        cms.PSet(
            min = cms.untracked.double(0.0),
            max = cms.untracked.double(100.0),
            nbins = cms.untracked.int32(100),
            name = cms.untracked.string("pt"),
            description = cms.untracked.string("pt [GeV/c]"),
            plotquantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            min = cms.untracked.double(-2.5),
            max = cms.untracked.double(2.5),
            nbins = cms.untracked.int32(50),
            name = cms.untracked.string("eta"),
            description = cms.untracked.string("eta"),
            plotquantity = cms.untracked.string("eta")
        )
    )
)

# TFileService for histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("miniAODCaloMuon_histos.root")
)

# Output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('miniAODCaloMuonMerger_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_miniAODCaloMuonMerger_*_*',
        'keep *_slimmedMuons_*_*',  # Keep original for comparison
        'keep *_packedPFCandidates_*_*',  # Keep for debugging
    )
)

# Path and EndPath definitions
process.p = cms.Path(
    process.miniAODCaloMuonMerger *
    process.muonAnalyzer
)

# Uncomment if using the unpacker version:
# process.p = cms.Path(
#     process.miniAODCaloMuonSequence *
#     process.muonAnalyzer
# )

process.outpath = cms.EndPath(process.out)

# Schedule
process.schedule = cms.Schedule(process.p, process.outpath)

# Print summary and enable threaded framework
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)