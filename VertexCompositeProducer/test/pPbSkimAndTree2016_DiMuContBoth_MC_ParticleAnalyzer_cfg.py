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
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/NonPromptPsi1S2S_pPb-EmbEPOS_8p16TeV_Pythia/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v2/50000/E47B68A7-05D6-E911-8740-1866DAEECF18.root'),
   #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/Psi1S_PbP-EmbEPOS_8p16TeV_Pythia/AODSIM/PbPEmb_80X_mcRun2_pA_v4-v3/60000/D6DB2977-97E8-E911-A850-D8D385AF8A88.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_pA_v4')

# Add the VertexComposite producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles
process.generalMuMuCandidates = generalParticles.clone(
    pdgId = cms.int32(443),
    doSwap = cms.bool(False),
    width = cms.double(999999999999.),
    finalSelection = cms.string(''),
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(13), selection = cms.string('')),
        cms.PSet(pdgId = cms.int32(13), selection = cms.string('')),
    ]),
)
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)

# Define the analysis steps
process.dimurereco_step = cms.Path(process.patMuonSequence * process.generalMuMuCandidates)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.dimucontana_mc = particleAna_mc.clone(
  recoParticles = cms.InputTag("generalMuMuCandidates"),
  genParticles = cms.untracked.InputTag("genParticles"),
  selectEvents = cms.string("")
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana_mc.root'))
process.p = cms.EndPath(process.dimucontana_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.dimurereco_step,
    process.p
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
