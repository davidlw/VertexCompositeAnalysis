import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2017_ppRef)

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
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/RunIIpp5Spring18DR/Ups1SMM_TuneCUETP8M1_5p02TeV_pythia8/AODSIM/94X_mc2017_realistic_forppRef5TeV-v2/1110000/0268D471-3435-E911-B89D-B083FED3F2E8.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_forppRef5TeV_v1', '')

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalDiMuCandidatesWrongSign = process.generalDiMuCandidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add centrality producer
process.load("RecoHI.HiCentralityAlgos.pACentrality_cfi")
process.pACentrality.srcHFhits = cms.InputTag("reducedHcalRecHits:hfreco")
process.pACentrality.srcEBhits = cms.InputTag("reducedEcalRecHitsEB")
process.pACentrality.srcEEhits = cms.InputTag("reducedEcalRecHitsEE")
process.pACentrality.produceHFtowers = cms.bool(False)
process.pACentrality.produceETmidRapidity = cms.bool(False)
process.pACentrality.producePixelhits = cms.bool(False)
process.pACentrality.producePixelTracks = cms.bool(False)
process.cent_seq = cms.Sequence(process.pACentrality)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.NoScraping)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.cent_seq)
process.dimurereco_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidates)
process.dimurerecowrongsign_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidatesWrongSign)

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana_mc.root'))
process.p = cms.EndPath(process.dimucontana_mc * process.dimucontana_wrongsign_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_NoScraping = cms.Path(process.NoScraping)
process.Flag_pileupVertexFilterCut = cms.Path(process.olvFilter_pp5TeV_dz1p0)
process.Flag_pileupVertexFilterCutGplus = cms.Path(process.pileUpFilter_pp5TeV_Gplus)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter , process.Flag_NoScraping , process.Flag_pileupVertexFilterCut , process.Flag_pileupVertexFilterCutGplus ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
