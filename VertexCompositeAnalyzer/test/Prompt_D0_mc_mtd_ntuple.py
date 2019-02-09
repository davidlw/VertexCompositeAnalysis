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
'file:///afs/cern.ch/user/y/yousen/public/mtdresearch/CMSSW_10_4_0_mtd5/src/VertexCompositeAnalysis/VertexCompositeProducer/test/prompt_d0.root'
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

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "103X_upgrade2018_realistic_HI_v6"
#process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000") 
#process.GlobalTag.toGet.extend([ 
#    cms.PSet(record = cms.string("HeavyIonRcd"), 
#        tag = cms.string("CentralityTable_HFtowers200_HydjetTuneCP5MTD_v1040mtd4x1_mc"), 
#        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"), 
#        label = cms.untracked.string("HFtowers") 
#        ), 
#    ]) 
#process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi') 
#process.hiCentrality.produceHFhits = False 
#process.hiCentrality.produceHFtowers = False 
#process.hiCentrality.produceEcalhits = False 
#process.hiCentrality.produceZDChits = False 
#process.hiCentrality.produceETmidRapidity = False 
#process.hiCentrality.producePixelhits = False 
#process.hiCentrality.produceTracks = False 
#process.hiCentrality.producePixelTracks = False 
#process.hiCentrality.reUseCentrality = True 
#process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO") 
#process.hiCentrality.srcTracks = cms.InputTag("generalTracks") 
#process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices") 
#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
#process.centralityBin.Centrality = cms.InputTag("hiCentrality") 
#process.centralityBin.centralityVariable = cms.string("HFtowers") 
#process.centralityBin.nonDefaultGlauberModel = cms.string("") 

process.d0ana.isUseMtd = cms.untracked.bool(True)
process.d0ana.doRecoNtuple = cms.untracked.bool(True)
process.d0ana.doGenNtuple = cms.untracked.bool(True)
process.d0ana.doGenMatching = cms.untracked.bool(True)
process.d0ana.doGenMatchingTOF = cms.untracked.bool(True)
process.d0ana.decayInGen = cms.untracked.bool(True)
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
process.d0preselector = process.d0selectorCutNew.clone()

#process.d0preselector = process.d0selectorCutNew.clone(
#  trkPtErrMax = cms.untracked.double(0.1),
#  trkNHitMin = cms.untracked.int32(11)
#)

process.d0selector.useAnyMVA = cms.bool(False)
#process.d0selector.VertexCompositeCollection = cms.untracked.InputTag("d0preselector:D0")

#process.d0ana_seq = cms.Sequence(process.d0preselector * process.d0selector * process.d0ana)
process.d0ana_seq = cms.Sequence( process.d0selector * process.d0ana)
#process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

#process.p = cms.Path(cent_seq * process.d0ana_seq)
process.p = cms.Path( process.d0ana_seq)
