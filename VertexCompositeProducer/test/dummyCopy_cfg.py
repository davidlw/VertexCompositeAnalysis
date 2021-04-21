import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it///store/data/Run2018C/JetHT/MINIAOD/12Nov2019_UL2018_rsb-v1/10000/015CDDF7-DBCB-094D-861E-54F73012741A.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_AK8PFJet500_v*', 
    ]

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter 
)
process.eventFilter_step = cms.Path( process.eventFilter )

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('ppRun2UL_MINIAOD.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('MINIAOD')
    )
)
process.output.outputCommands = cms.untracked.vstring('keep *'
)

process.output_step = cms.EndPath(process.output)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.output_step
)
