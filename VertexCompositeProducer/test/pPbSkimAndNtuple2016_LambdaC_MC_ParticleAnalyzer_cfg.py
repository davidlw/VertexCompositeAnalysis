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
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/LambdaC-KsPr_LCpT-5p9_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/70000/30A92DC7-C99C-E711-8E53-0242AC110003.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('80X_mcRun2_pA_v4')

# Add the VertexComposite producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

process.kShort = generalParticles.clone(
    pdgId = cms.int32(310),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(0.2),

    preSelection = cms.string(""),
    pocaSelection = cms.string("pt >= 1.0 && abs(rapidity) < 2.4"),
    postSelection = cms.string(""),
    preMassSelection = cms.string(""),
    finalSelection = cms.string( "abs(userFloat('angle3D'))<0.2 && abs(userFloat('lVtxSig'))>2.5"),
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
              selection = cms.string(
              "pt>1.0 && abs(eta)<2.4"
              "&& quality('loose')"" && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& numberOfValidHits >=4")
            ),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
              selection = cms.string(
              "pt>1.0 && abs(eta)<2.4"
              "&& quality('loose')"" && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& numberOfValidHits >=4")
            ),
    ]),
)
process.LambdaC = generalParticles.clone(
    pdgId = cms.int32(4122),
    doSwap = cms.bool(False),
    preMassSelection = cms.string("abs(charge)==1"),
    finalSelection = cms.string(''),
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(310), source = cms.InputTag('kShort'), selection = cms.string('')),
        cms.PSet(pdgId = cms.int32(2212), #charge = cms.int32(+1),
          selection = cms.string("pt>1.0 && abs(eta)<2.4"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& numberOfValidHits >=11")),
    ]),
)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter * process.primaryVertexFilterPA * process.NoScraping * process.olvFilter_pPb8TeV_dz1p0)

# Define the analysis steps
process.rereco_step = cms.Path(process.kShort * process.LambdaC)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
process.lambdacAna_mc = particleAna_mc.clone(
  recoParticles = cms.InputTag("LambdaC"),
  genParticles = cms.untracked.InputTag("genParticles"),
  selectEvents = cms.string(""),
  addSource    = cms.untracked.bool(False),
  genPdgId     = cms.untracked.vuint32([4122, 310, 2212, 211]),
  saveTree = cms.untracked.bool(False)
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('lambdacana_mc.root'))
process.p = cms.EndPath(process.lambdacAna_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.rereco_step,
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
