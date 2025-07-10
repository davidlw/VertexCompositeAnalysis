# CaloMuonMerger PAT Adaptation Solutions

This package provides solutions for adapting the CMSSW CaloMuonMerger module to work with `pat::Muon` objects instead of `reco::Muon` objects.

## Problem Description

The original `CaloMuonMerger.cc` module in CMSSW works with `reco::Muon` objects as input, but many analyses require `pat::Muon` objects which provide additional analysis-level information including:

- Enhanced isolation information
- MC truth matching
- Trigger matching
- Embedded tracks and objects
- User data storage capabilities

## Solutions Provided

### Solution 1: Native PAT CaloMuon Merger (`PatCaloMuonMerger.cc`)

A complete rewrite of the CaloMuonMerger specifically designed to work with PAT muons.

**Features:**
- Takes `pat::MuonCollection` as primary input
- Merges with `reco::CaloMuonCollection` 
- Produces enhanced `pat::MuonCollection` output
- Preserves all PAT muon functionality
- Adds calo muon information as user data

**Usage:**
```python
process.patCaloMuonMerger = cms.EDProducer("PatCaloMuonMerger",
    muons = cms.InputTag("slimmedMuons"),
    caloMuons = cms.InputTag("calomuons"),
    tracks = cms.InputTag("generalTracks"),
    minCaloCompatibility = cms.double(0.6),
    deltaR = cms.double(0.1),
    embedCaloMuon = cms.bool(True)
)
```

### Solution 2: Adapter Approach (`RecoToPatMuonAdapter.cc`)

An adapter module that converts the output of the existing CaloMuonMerger to PAT format.

**Features:**
- Works with existing CaloMuonMerger without modifications
- Converts `reco::MuonCollection` to `pat::MuonCollection`
- Preserves all muon information
- Adds quality flags as user data
- Minimal code changes required

**Usage:**
```python
# Use existing CaloMuonMerger
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons

# Add adapter to convert to PAT
process.recoToPatMuonAdapter = cms.EDProducer("RecoToPatMuonAdapter",
    src = cms.InputTag("calomuons"),
    preserveUserData = cms.bool(True),
    addQualityFlags = cms.bool(True)
)

# Run in sequence
process.p = cms.Path(calomuons * process.recoToPatMuonAdapter)
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

| Aspect | PatCaloMuonMerger | RecoToPatMuonAdapter |
|--------|-------------------|---------------------|
| **Complexity** | More complex, native PAT | Simple adapter |
| **Performance** | Single-pass processing | Two-pass (reco→PAT) |
| **Maintenance** | Independent of original | Depends on original |
| **Features** | Full PAT integration | Basic conversion |
| **Recommended for** | New analyses | Existing workflows |

## Advanced Configuration

### PatCaloMuonMerger Parameters

- `muons`: Input PAT muon collection (default: "slimmedMuons")
- `caloMuons`: Input calo muon collection (default: "calomuons")
- `tracks`: Input track collection (default: "generalTracks")
- `minCaloCompatibility`: Minimum calo compatibility threshold (default: 0.6)
- `deltaR`: Matching cone size for duplicate removal (default: 0.1)
- `embedCaloMuon`: Whether to embed calo muon info as user data (default: true)

### RecoToPatMuonAdapter Parameters

- `src`: Input reco muon collection (default: "muons")
- `preserveUserData`: Preserve existing user data (default: true)
- `addQualityFlags`: Add muon quality flags as user data (default: true)

## Accessing Added Information

When using either solution, additional information is stored as user data in the PAT muons:

```cpp
// C++ example
const pat::Muon& muon = ...;

// Check if calo muon information is available
if (muon.hasUserFloat("caloCompatibility")) {
    float caloComp = muon.userFloat("caloCompatibility");
}

// Access quality flags (from adapter)
if (muon.hasUserInt("isGlobalMuon")) {
    bool isGlobal = muon.userInt("isGlobalMuon") > 0;
}
```

```python
# Python/PyROOT example
caloCompatibility = muon.userFloat("caloCompatibility")
isGlobalMuon = muon.userInt("isGlobalMuon") > 0
```

## Testing

Use the provided test configuration files:

```bash
# Test PatCaloMuonMerger
cmsRun test/patCaloMuonMerger_cfg.py

# Test with adapter approach
cmsRun test/caloMuonMergerToPat_cfg.py
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
