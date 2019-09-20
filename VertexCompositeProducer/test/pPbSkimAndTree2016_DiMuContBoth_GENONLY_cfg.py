import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2016_pA)

process.load('Configuration.StandardSequences.Services_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/NonPrPsi1_2SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/90000/FAE72020-291B-E711-A027-00237DF28498.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_pA_v4')

# Add the GEN muons
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Define the analysis steps
process.dimugen_step = cms.Path(process.genMuons)

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")
process.dimuana_mc.doRecoNtuple = cms.untracked.bool(False)
                #doMuon = cms.untracked.bool(False),
                 # doMuonFull = cms.untracked.bool(False),
                  #  isCentrality = cms.bool(False),

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana_GENONLY.root'))
process.p = cms.EndPath(process.dimuana_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.dimugen_step,
    process.p
)
