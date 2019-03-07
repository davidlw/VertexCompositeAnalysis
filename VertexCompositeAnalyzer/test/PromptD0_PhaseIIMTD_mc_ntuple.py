import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'/store/user/yousen/D0_PiK_prompt_pt0_y4_5p5TeV_TuneCP5_Pythia8/promptd0_mc_mtd_Skim_v1/190207_191503/0000/prompt_d0_74.root'
#'file:///afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_10_4_0_mtd5/src/VertexCompositeAnalysis/VertexCompositeProducer/test/prompt_d0.root'
#'/store/user/davidlw/D0_PiK_prompt_pt0_y4_5p5TeV_TuneCP5_Pythia8/promptd0_mc_mtd_Skim_v5/190221_115230/0000/prompt_d0_99.root'
'file://prompt_d0.root'
#'root://xrootd-cms.infn.it//store/user/yousen/D0_PiK_prompt_pt0_y4_5p5TeV_TuneCP5_Pythia8/promptd0_mc_mtd_Skim_v6/190227_003058/0000/prompt_d0_11.root',
                ),
#secondaryFileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/98B099A2-2DB0-E611-BE2C-02163E0119D7.root'
#)
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('promptd0_mc_mtd.root')
                                   )

process.d0ana_mc.isUseMtd = cms.untracked.bool(True)
process.d0ana_mc.doRecoNtuple = cms.untracked.bool(True)
process.d0ana_mc.doGenNtuple = cms.untracked.bool(True)
process.d0ana_mc.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")
process.d0ana_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMC:D0")
process.d0ana_mc.MVACollection = cms.InputTag("d0selectorMC:MVAValuesNewD0")
process.d0ana_mc.isCentrality = cms.bool(True)

process.d0selectorMC.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")

process.d0ana_seq = cms.Sequence( process.d0selectorMC * process.d0ana_mc)

process.p = cms.Path( process.d0ana_seq)
