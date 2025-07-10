# CMSSW_15_0_9 Compatibility Guide

This document outlines the specific changes and considerations for running the MiniAOD CaloMuon merger modules in CMSSW_15_0_9.

## Key Changes for CMSSW_15_0_9

### 1. Global Tag Updates
- **Old (Run 2)**: `auto:run2_mc`
- **New (Run 3/2024)**: `auto:phase1_2024_realistic`

For data analyses, use:
- **Run 3 Data**: `auto:phase1_2024_data`
- **Run 3 MC**: `auto:phase1_2024_realistic`

### 2. Threading Framework
CMSSW_15_0_9 uses an improved threading framework. The configurations now include:
```python
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)
```

### 3. Updated Headers
Added compatibility headers:
```cpp
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
```

### 4. BuildFile.xml Updates
Added dependencies for CMSSW_15_0_9:
```xml
<use name="FWCore/Utilities"/>
<use name="DataFormats/Candidate"/>
```

## Installation for CMSSW_15_0_9

```bash
# Setup CMSSW_15_0_9
cmsrel CMSSW_15_0_9
cd CMSSW_15_0_9/src
cmsenv

# Create your package
mkdir -p MyAnalysis/MiniAODCaloMuonMerger
cd MyAnalysis/MiniAODCaloMuonMerger

# Copy the source files
cp /path/to/MiniAODCaloMuonMerger.cc plugins/
cp /path/to/MiniAODCaloMuonMergerWithUnpacker.cc plugins/
cp /path/to/BuildFile.xml plugins/
cp /path/to/python/*.py python/
cp /path/to/test/*.py test/

# Compile
cd $CMSSW_BASE/src
scram b -j8
```

## Run 3 Specific Considerations

### Data Format Changes
Run 3 data may have some differences in:
- **Packed candidate content**: New packing algorithms
- **Muon reconstruction**: Improved algorithms in Run 3
- **Isolation variables**: Enhanced PF isolation

### Trigger Updates
If using trigger information, note that Run 3 triggers have updated:
- HLT path names
- L1 trigger objects
- Prescale handling

### Performance Improvements
CMSSW_15_0_9 includes:
- Better memory management
- Improved packed candidate handling
- Enhanced deltaR calculations

## Example Configurations

### For Run 3 Data Analysis
```python
# Global tag for Run 3 data
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_data', '')

# Input files (example)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2023C/Muon/MINIAOD/...'
    )
)
```

### For Run 3 MC Analysis
```python
# Global tag for Run 3 MC
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# Input files (example)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Run3Summer23MiniAODv4/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/...'
    )
)
```

## Validation

### Quick Test
```bash
# Test with Run 3 MC
cmsRun test/miniAODCaloMuonMerger_cfg.py maxEvents=10

# Comprehensive validation (recommended)
cmsRun test/validate_CMSSW_15_0_9_cfg.py

# Check output
edmDumpEventContent miniAODCaloMuonMerger_output.root
```

### Expected Output
The modules should produce:
- Enhanced PAT muon collection with calo information
- User data fields for calo compatibility
- Debug information in the log

## Troubleshooting

### Common Issues in CMSSW_15_0_9

1. **Global Tag Errors**
   - Error: `Global tag not found`
   - Solution: Use correct Run 3 global tags

2. **Threading Issues**
   - Error: `Stream module issues`
   - Solution: Ensure proper threading configuration

3. **Data Format Mismatches**
   - Error: `Product not found`
   - Solution: Check input tag names for Run 3 format

### Debug Mode
Add to your configuration:
```python
process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.categories.append('MiniAODCaloMuonMerger')
```

## Performance Optimization

For large-scale production:
```python
# Enable multi-threading (adjust based on your resources)
process.options.numberOfThreads = cms.untracked.uint32(4)
process.options.numberOfStreams = cms.untracked.uint32(0)  # Auto-determine

# Output optimization
process.out.compressionLevel = cms.untracked.int32(9)
process.out.compressionAlgorithm = cms.untracked.string("LZMA")
```

## Version Compatibility Matrix

| CMSSW Version | Status | Notes |
|---------------|--------|-------|
| CMSSW_15_0_9 | ✅ Fully Supported | Primary target |
| CMSSW_14_1_X | ✅ Compatible | Minor modifications needed |
| CMSSW_13_3_X | ⚠️ Limited | Requires header changes |
| CMSSW_12_X_X | ❌ Not supported | Major API changes |

## Contact and Support

For CMSSW_15_0_9 specific issues:
1. Check the official CMSSW documentation
2. Review the latest PAT and miniAOD format changes
3. Consult the CMS Muon POG for algorithm updates