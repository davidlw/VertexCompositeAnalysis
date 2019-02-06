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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'file:///afs/cern.ch/user/y/yousen/public/mtdreserach/CMSSW_10_4_0_mtd5/src/VertexCompositeAnalysis/VertexCompositeProducer/test/prompt_d0.root'
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

process.d0ana.isUseMtd = cms.untracked.bool(True);
process.d0ana.doRecoNtuple = cms.untracked.bool(True);
process.d0ana.doGenNtuple = cms.untracked.bool(False);
process.d0ana.useAnyMVA = cms.bool(False)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveAllHistogram = cms.untracked.bool(True)
process.d0ana.saveTree = cms.untracked.bool(False)

process.d0selectorCutNew = process.d0selector.clone()
process.d0selectorCutNew.cand3DDecayLengthSigMin = cms.untracked.double(-10000.)
process.d0selectorCutNew.cand3DPointingAngleMax = cms.untracked.double(10000.)
process.d0selectorCutNew.candVtxProbMin = cms.untracked.double(-10000.) 
process.d0preselector = process.d0selectorCutNew.clone()

#process.d0preselector = process.d0selectorCutNew.clone(
#  trkPtErrMax = cms.untracked.double(0.1),
#  trkNHitMin = cms.untracked.int32(11)
#)

process.d0selector.useAnyMVA = cms.bool(False)
process.d0selector.VertexCompositeCollection = cms.untracked.InputTag("d0preselector:D0")

process.d0ana_seq = cms.Sequence(process.d0preselector * process.d0selector * process.d0ana)
#process.d0ana_seq = cms.Sequence(process.d0selector * process.d0ana)

process.p = cms.Path(process.d0ana_seq)
