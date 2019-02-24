import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'/store/user/yousen/Hydjet_5p02TeV_TuneCP5_MTD/hyjets_mc_mtd_Skim_v1/190207_195017/0000/hyjets_10.root'
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets10.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets09.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets08.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets07.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets06.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets06.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets05.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets04.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets03.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets02.root',
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/producer/hyjets01.root'
                ),
                            )
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
#cms.string('hyjets_mc_mtd.root')
cms.string('/afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/mtd_store/analyzer/hyjets_mc_mtd.root')
                                   )

process.d0ana.isUseMtd = cms.untracked.bool(True);
process.d0ana.doRecoNtuple = cms.untracked.bool(True);
#process.d0ana.doGenNtuple = cms.untracked.bool(True);
process.d0ana.doGenNtuple = cms.untracked.bool(False);
#process.d0ana.doGenMatching = cms.untracked.bool(True);
#process.d0ana.doGenMatchingTOF = cms.untracked.bool(True);
process.d0ana.useAnyMVA = cms.bool(False)
process.d0ana.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveAllHistogram = cms.untracked.bool(False)
process.d0ana.saveTree = cms.untracked.bool(True)
process.d0ana.isCentrality = cms.bool(True)

process.d0selector.VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices4D")

process.d0selectorCutNew = process.d0selector.clone()
process.d0selectorCutNew.cand3DDecayLengthSigMin = cms.untracked.double(-10000.)
process.d0selectorCutNew.cand3DPointingAngleMax = cms.untracked.double(10000.)
process.d0selectorCutNew.candVtxProbMin = cms.untracked.double(-10000.) 
#process.d0preselector = process.d0selectorCutNew.clone()

#process.d0preselector = process.d0selectorCutNew.clone(
#  trkPtErrMax = cms.untracked.double(0.1),
#  trkNHitMin = cms.untracked.int32(11)
#)

process.d0selector.useAnyMVA = cms.bool(False)
#process.d0selector.VertexCompositeCollection = cms.untracked.InputTag("d0preselector:D0")

#process.d0ana_seq = cms.Sequence(process.d0preselector * process.d0selector * process.d0ana)
process.d0ana_seq = cms.Sequence(process.d0selector * process.d0ana)

process.p = cms.Path(process.d0ana_seq)
