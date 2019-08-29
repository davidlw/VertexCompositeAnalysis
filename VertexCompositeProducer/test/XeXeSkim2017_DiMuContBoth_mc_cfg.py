import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2017_pp_on_XeXe)

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
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/himc/HINXeXeFall17DR/MinBias-pp_XeXeEmb-HydjetCymbEv8_5p44_Pythia8/GEN-SIM-RECO/94X_mc2017_realistic_ForXeXe_v7-v3/00000/E6274BAE-200C-E811-95E7-1866DAEEB358.root'),
   inputCommands=cms.untracked.vstring('keep *')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_ForXeXe_v8', '')

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalDiMuCandidatesWrongSign = process.generalDiMuCandidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add XeXe centrality
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.nonDefaultGlauberModel = cms.string("HYDJET")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_XeXe5p44TeVHYDJET_v941x01_mc"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowersHYDJET")
        ),
    ])
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.NoScraping * process.hfCoincFilter3Th3)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.cent_seq)
process.dimurereco_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidates)
process.dimurerecowrongsign_step = cms.Path(process.patMuonSequence * process.generalDiMuCandidatesWrongSign)

# Define the output
process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('XeXe_DiMuCont_MC.root'),
)
process.output_HM_step = cms.EndPath(process.output_HM)

# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.output_HM_step
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_NoScraping = cms.Path(process.NoScraping)
process.Flag_hfCoincFilter3Th3 = cms.Path(process.hfCoincFilter3Th3)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter , process.Flag_NoScraping , process.Flag_hfCoincFilter3Th3 ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)
