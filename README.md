# CaloMuonMerger MiniAOD Adaptation Solutions

This package provides solutions for adapting the CMSSW CaloMuonMerger functionality to work with miniAOD files where muons are stored as `slimmedMuons` (packed candidates).

## Problem Description

The original `CaloMuonMerger.cc` module in CMSSW works with:
- `reco::Muon` collections as input
- `reco::CaloMuon` collections 
- `reco::Track` collections

However, in miniAOD files:
- Muons are stored as `slimmedMuons` (packed candidates)
- There is no separate `reco::CaloMuon` collection
- Tracks are also packed in `packedPFCandidates`
- All objects are in PAT format for analysis

This creates a mismatch between the CaloMuonMerger expectations and the miniAOD format.

## Solutions Provided

### Solution 1: Direct MiniAOD CaloMuon Merger (`MiniAODCaloMuonMerger.cc`)

A standalone module that works directly with miniAOD format, extracting calo muon information from packed candidates.

**Features:**
- Takes `slimmedMuons` as primary input
- Analyzes `packedPFCandidates` for calo muon signatures
- Produces enhanced `pat::MuonCollection` output
- Identifies potential calo muons based on calo fraction and track quality
- No dependency on unpacking modules

**Usage:**
```python
process.miniAODCaloMuonMerger = cms.EDProducer("MiniAODCaloMuonMerger",
    muons = cms.InputTag("slimmedMuons"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    addCaloMuonsFromPacked = cms.bool(True),
    enhanceExistingMuons = cms.bool(True)
)
```

### Solution 2: Integration with Existing Unpackers (`MiniAODCaloMuonMergerWithUnpacker.cc`)

A module that integrates with your existing `MuonUnpacker` and `TrackAndVertexUnpacker` modules.

**Features:**
- Works with output from your existing `MuonUnpacker.cc`
- Uses unpacked tracks from your `TrackAndVertexUnpacker.cc`
- Enhanced calo compatibility calculation
- Recovers missing calo muons from packed candidates
- Seamless integration with your current workflow

**Usage:**
```python
# Your existing sequence
process.load("VertexCompositeProducer.VertexCompositeProducer.trackAndVertexUnpacker_cfi")
process.load("VertexCompositeProducer.VertexCompositeProducer.muonUnpacker_cfi")

# Add calo muon merger
process.miniAODCaloMuonMergerWithUnpacker = cms.EDProducer("MiniAODCaloMuonMergerWithUnpacker",
    unpackedMuons = cms.InputTag("muonUnpacker"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    unpackedTracks = cms.InputTag("unpackedTracksAndVertices"),
    recalculateCaloCompatibility = cms.bool(True),
    addMissingCaloMuons = cms.bool(True)
)

# Complete sequence
process.p = cms.Path(
    process.unpackedTracksAndVertices *
    process.muonUnpacker *
    process.miniAODCaloMuonMergerWithUnpacker
)
```

## Installation Instructions

1. **Setup CMSSW environment:**
   ```bash
   cmsrel CMSSW_X_Y_Z
   cd CMSSW_X_Y_Z/src
   cmsenv
   ```

2. **Create your package structure:**
   ```bash
   mkdir -p YourPackage/YourSubsystem
   cd YourPackage/YourSubsystem
   ```

3. **Copy the source files:**
   - Copy the `.cc` files to `plugins/`
   - Copy the `BuildFile.xml` to `plugins/`
   - Copy the Python configs to `python/`
   - Copy the test config to `test/`

4. **Compile:**
   ```bash
   cd $CMSSW_BASE/src
   scram b -j8
   ```

## Key Differences Between Solutions

| Aspect | MiniAODCaloMuonMerger | MiniAODCaloMuonMergerWithUnpacker |
|--------|----------------------|----------------------------------|
| **Complexity** | Moderate, self-contained | Simple, builds on existing |
| **Dependencies** | Only miniAOD collections | Requires your unpacker modules |
| **Performance** | Direct miniAOD processing | Uses pre-unpacked objects |
| **Integration** | Standalone solution | Seamless with existing workflow |
| **Calo Detection** | Based on packed candidate properties | Enhanced algorithm using unpacked info |
| **Recommended for** | New analyses, simple setup | Existing workflows using unpackers |

## Advanced Configuration

### MiniAODCaloMuonMerger Parameters

- `muons`: Input slimmed muon collection (default: "slimmedMuons")
- `packedCandidates`: Input packed candidate collection (default: "packedPFCandidates")
- `tracks`: Optional unpacked track collection
- `minCaloCompatibility`: Minimum calo compatibility threshold (default: 0.6)
- `deltaR`: Matching cone size for duplicate removal (default: 0.1)
- `minPt`: Minimum pT threshold for new calo muons (default: 2.0)
- `maxEta`: Maximum |η| for new calo muons (default: 2.4)
- `addCaloMuonsFromPacked`: Add new calo muons from packed candidates (default: true)
- `enhanceExistingMuons`: Enhance existing muons with calo info (default: true)
- `requireTrackerTrack`: Require track details for calo muon candidates (default: false)

### MiniAODCaloMuonMergerWithUnpacker Parameters

- `unpackedMuons`: Input from MuonUnpacker (default: "muonUnpacker")
- `packedCandidates`: Input packed candidate collection (default: "packedPFCandidates")
- `unpackedTracks`: Input from TrackAndVertexUnpacker (default: "unpackedTracksAndVertices")
- `minCaloCompatibility`: Minimum calo compatibility threshold (default: 0.6)
- `deltaR`: Matching cone size (default: 0.1)
- `minPt`: Minimum pT threshold (default: 2.0)
- `recalculateCaloCompatibility`: Use enhanced calo compatibility calculation (default: true)
- `addMissingCaloMuons`: Add calo muons missing from unpacked collection (default: true)

## Accessing Added Information

When using either solution, additional information is stored as user data in the PAT muons:

```cpp
// C++ example
const pat::Muon& muon = ...;

// Check if enhanced calo muon information is available
if (muon.hasUserFloat("caloCompatibility")) {
    float caloComp = muon.userFloat("caloCompatibility");
}

// For MiniAODCaloMuonMerger
if (muon.hasUserFloat("caloFraction")) {
    float caloFrac = muon.userFloat("caloFraction");
}

if (muon.hasUserInt("fromPackedCandidate")) {
    bool fromPacked = muon.userInt("fromPackedCandidate") > 0;
}

// For MiniAODCaloMuonMergerWithUnpacker
if (muon.hasUserFloat("enhancedCaloCompatibility")) {
    float enhancedComp = muon.userFloat("enhancedCaloCompatibility");
}

if (muon.hasUserInt("isLikelyCaloMuon")) {
    bool likelyCalo = muon.userInt("isLikelyCaloMuon") > 0;
}
```

```python
# Python/PyROOT example
# Basic calo information
caloCompatibility = muon.userFloat("caloCompatibility")
caloFraction = muon.userFloat("caloFraction")

# Enhanced information (with unpacker)
enhancedCaloComp = muon.userFloat("enhancedCaloCompatibility")
isLikelyCaloMuon = muon.userInt("isLikelyCaloMuon") > 0
```

## Testing

Use the provided test configuration files:

```bash
# Test direct miniAOD approach
cmsRun test/miniAODCaloMuonMerger_cfg.py

# Test with your existing unpackers (adjust paths as needed)
cmsRun test/miniAODCaloMuonMergerWithUnpacker_cfg.py
```

## Troubleshooting

### Common Issues

1. **Compilation errors**: Ensure all dependencies are properly declared in `BuildFile.xml`
2. **Missing input collections**: Check that input tags match your data format (AOD vs MiniAOD)
3. **Empty output collections**: Verify calo compatibility thresholds and deltaR cuts

### Performance Considerations

- For large-scale processing, prefer `PatCaloMuonMerger` for better performance
- For quick prototyping or existing workflows, `RecoToPatMuonAdapter` is sufficient
- Monitor memory usage when embedding large amounts of user data

## Support and Contributing

For questions or contributions:
1. Check CMSSW documentation for PAT and muon systems
2. Consult the CMS Muon POG for analysis recommendations
3. Test thoroughly with your specific use case before production use

## References

- [CMSSW Muon Analysis Workbook](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis)
- [PAT Documentation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePAT)
- [CMS Muon Identification](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId) 
