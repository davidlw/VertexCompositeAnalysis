import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("ANATREE",eras.Run2_2018_pp_on_AA)

# Limit the output messages
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )

# Define the input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/anstahll/RiceHIN/PbPb2018/SKIM/VertexCompositeSkim_HIDoubleMuonPsiPeri_v2_HIRun2018_DiMuMassMin7_20190407/HIDoubleMuonPsiPeri/VertexCompositeSkim_HIDoubleMuonPsiPeri_v2_HIRun2018_DiMuMassMin7_20190407/190407_181439/0001/PbPb_DiMuCont_1249.root'
'file:/afs/cern.ch/user/a/anstahll/work/Analysis/GITHUB/CMSSW_10_3_3_patch1/src/VertexCompositeAnalysis/VertexCompositeProducer/test/PbPb_DiMuCont.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_v6')

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

# Define the output
process.TFileService = cms.Service("TFileService", fileName =
cms.string('dimuana.root')
)
process.p = cms.EndPath(process.dimucontana * process.dimucontana_wrongsign)
