import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("ANATREE",eras.Run2_2016_pA)

# Limit the output messages
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )

# Define the input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/a/anstahll/work/Analysis/Production/CMSSW_8_0_32/src/VertexCompositeAnalysis/VertexCompositeProducer/test/pPb_DiMuCont.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_dataRun2_v19')

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

# Define the output
process.TFileService = cms.Service("TFileService", fileName =
cms.string('dimuana.root')
)
process.p = cms.EndPath(process.dimucontana * process.dimucontana_wrongsign)
