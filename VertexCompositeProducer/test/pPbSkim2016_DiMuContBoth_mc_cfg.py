import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2016_pA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/AODSIM/pPbBst_80X_mcRun2_pA_v4-v1/50000/C4CEAAAB-901A-E711-AB5E-0CC47A7C3420.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_pA_v4')

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalDiMuCandidatesWrongSign = process.generalDiMuCandidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)

# Define the analysis steps
process.dimurereco_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidates)
process.dimurerecowrongsign_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidatesWrongSign)

# Define the output
process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pPb_DiMuCont_MC.root'),
)
process.output_HM_step = cms.EndPath(process.output_HM)

# Define the process schedule
process.schedule = cms.Schedule(
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.output_HM_step
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_hfCoincFilter = cms.Path(process.hfCoincFilter)
process.Flag_primaryVertexFilterPA = cms.Path(process.primaryVertexFilterPA)
process.Flag_NoScraping = cms.Path(process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.olvFilter_pPb8TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.pileUpFilter_pPb8TeV_Gplus)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter , process.Flag_primaryVertexFilterPA , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
