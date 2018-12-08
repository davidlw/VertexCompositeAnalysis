import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(200)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) 
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "103X_dataRun2_Prompt_v3"
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000") 
process.GlobalTag.toGet.extend([ 
    cms.PSet(record = cms.string("HeavyIonRcd"), 
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"), 
        label = cms.untracked.string("HFtowers") 
        ), 
    ]) 
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
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
process.centralityBin.Centrality = cms.InputTag("hiCentrality") 
process.centralityBin.centralityVariable = cms.string("HFtowers") 
process.centralityBin.nonDefaultGlauberModel = cms.string("") 

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/user/davidlw/HIMinimumBias0/Skim_D0_Nov20DCS_v1/181121_103527/0000/PbPb_100.root'
'file:/afs/cern.ch/user/d/davidlw/CMSSW/CMSSW_10_3_1_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/PbPb.root'
                ),
secondaryFileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/PromptReco-v1/000/326/577/00000/532A6440-350D-2446-89FB-3CA9E8335E4F.root'
#'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/389/00000/8EE40790-9E67-E441-8663-2CDBC2F8E53C.root'
)
                            )

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity185_*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana.root')
                                   )

process.d0ana.useAnyMVA = cms.bool(True)
process.d0ana.VertexCompositeCollection = cms.untracked.InputTag("d0selector:D0")
process.d0ana.MVACollection = cms.InputTag("d0selector:MVAValuesNewD0")
#process.d0ana.isSkimMVA = cms.untracked.bool(True)
process.d0ana.saveHistogram = cms.untracked.bool(True)
process.d0ana.saveTree = cms.untracked.bool(False)
#process.d0ana.yBins = cms.untracked.vdouble(-2.4,-1.6,-0.8,0.0,0.8,1.6,2.4)
process.d0ana.yBins = cms.untracked.vdouble(-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0)
process.d0ana.pTBins = cms.untracked.vdouble(0,1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0)

process.d0selectorBDTPreCut.trkPtSumMin = cms.untracked.double(2.2)
process.d0selectorBDTPreCut.trkPtMin = cms.untracked.double(1)
process.d0selector = process.d0selectorBDTPreCut.clone()
process.d0selector.useAnyMVA = cms.bool(True)
#process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InPbPb_ptsum2p2_pt0p7_PbPbMB_WS.root')
process.d0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_PromptD0InPbPb_ptsum2p2_pt1_PbPbMB_WS.root')
process.d0selector.GBRForestLabel = cms.string('D0InPbPb')
process.d0selector.centMin = cms.untracked.int32(0)
process.d0selector.centMax = cms.untracked.int32(200)
process.d0selector.isCentrality = cms.bool(True)

process.npd0selector = process.d0selector.clone()
#process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InPbPb_ptsum2p2_pt0p7_PbPbMB_WS.root')
process.npd0selector.GBRForestFileName = cms.string('GBRForestfile_BDT_NonPromptD0InPbPb_ptsum2p2_pt1_PbPbMB_WS.root')
process.npd0ana = process.d0ana.clone()
process.npd0ana.VertexCompositeCollection = cms.untracked.InputTag("npd0selector:D0")
process.npd0ana.MVACollection = cms.InputTag("npd0selector:MVAValuesNewD0")

process.d0selector05 = process.d0selector.clone()
process.d0selector510 = process.d0selector.clone()
process.d0selector1030 = process.d0selector.clone()
process.d0selector3050 = process.d0selector.clone()
process.d0selector5080 = process.d0selector.clone()
process.npd0selector05 = process.npd0selector.clone()
process.npd0selector510 = process.npd0selector.clone()
process.npd0selector1030 = process.npd0selector.clone()
process.npd0selector3050 = process.npd0selector.clone()
process.npd0selector5080 = process.npd0selector.clone()

process.d0selector05.centMin = cms.untracked.int32(0)
process.d0selector05.centMax = cms.untracked.int32(10)
process.d0selector510.centMin = cms.untracked.int32(10)
process.d0selector510.centMax = cms.untracked.int32(20)
process.d0selector1030.centMin = cms.untracked.int32(20)
process.d0selector1030.centMax = cms.untracked.int32(60)
process.d0selector3050.centMin = cms.untracked.int32(60)
process.d0selector3050.centMax = cms.untracked.int32(100)
process.d0selector5080.centMin = cms.untracked.int32(100)
process.d0selector5080.centMax = cms.untracked.int32(160)

process.npd0selector05.centMin = cms.untracked.int32(0)
process.npd0selector05.centMax = cms.untracked.int32(10)
process.npd0selector510.centMin = cms.untracked.int32(10)
process.npd0selector510.centMax = cms.untracked.int32(20)
process.npd0selector1030.centMin = cms.untracked.int32(20)
process.npd0selector1030.centMax = cms.untracked.int32(60)
process.npd0selector3050.centMin = cms.untracked.int32(60)
process.npd0selector3050.centMax = cms.untracked.int32(100)
process.npd0selector5080.centMin = cms.untracked.int32(100)
process.npd0selector5080.centMax = cms.untracked.int32(160)

process.d0ana05 = process.d0ana.clone()
process.d0ana510 = process.d0ana.clone()
process.d0ana1030 = process.d0ana.clone()
process.d0ana3050 = process.d0ana.clone()
process.d0ana5080 = process.d0ana.clone()
process.npd0ana05 = process.npd0ana.clone()
process.npd0ana510 = process.npd0ana.clone()
process.npd0ana1030 = process.npd0ana.clone()
process.npd0ana3050 = process.npd0ana.clone()
process.npd0ana5080 = process.npd0ana.clone()

process.d0ana05.VertexCompositeCollection = cms.untracked.InputTag("d0selector05:D0")
process.d0ana05.MVACollection = cms.InputTag("d0selector05:MVAValuesNewD0")
process.d0ana510.VertexCompositeCollection = cms.untracked.InputTag("d0selector510:D0")
process.d0ana510.MVACollection = cms.InputTag("d0selector510:MVAValuesNewD0")
process.d0ana1030.VertexCompositeCollection = cms.untracked.InputTag("d0selector1030:D0")
process.d0ana1030.MVACollection = cms.InputTag("d0selector1030:MVAValuesNewD0")
process.d0ana3050.VertexCompositeCollection = cms.untracked.InputTag("d0selector3050:D0")
process.d0ana3050.MVACollection = cms.InputTag("d0selector3050:MVAValuesNewD0")
process.d0ana5080.VertexCompositeCollection = cms.untracked.InputTag("d0selector5080:D0")
process.d0ana5080.MVACollection = cms.InputTag("d0selector5080:MVAValuesNewD0")

process.npd0ana05.VertexCompositeCollection = cms.untracked.InputTag("npd0selector05:D0")
process.npd0ana05.MVACollection = cms.InputTag("npd0selector05:MVAValuesNewD0")
process.npd0ana510.VertexCompositeCollection = cms.untracked.InputTag("npd0selector510:D0")
process.npd0ana510.MVACollection = cms.InputTag("npd0selector510:MVAValuesNewD0")
process.npd0ana1030.VertexCompositeCollection = cms.untracked.InputTag("npd0selector1030:D0")
process.npd0ana1030.MVACollection = cms.InputTag("npd0selector1030:MVAValuesNewD0")
process.npd0ana3050.VertexCompositeCollection = cms.untracked.InputTag("npd0selector3050:D0")
process.npd0ana3050.MVACollection = cms.InputTag("npd0selector3050:MVAValuesNewD0")
process.npd0ana5080.VertexCompositeCollection = cms.untracked.InputTag("npd0selector5080:D0")
process.npd0ana5080.MVACollection = cms.InputTag("npd0selector5080:MVAValuesNewD0")

process.d0ana05_seq = cms.Sequence(process.d0selector05 * process.d0ana05)
process.d0ana510_seq = cms.Sequence(process.d0selector510 * process.d0ana510)
process.d0ana1030_seq = cms.Sequence(process.d0selector1030 * process.d0ana1030)
process.d0ana3050_seq = cms.Sequence(process.d0selector3050 * process.d0ana3050)
process.d0ana5080_seq = cms.Sequence(process.d0selector5080 * process.d0ana5080)

process.npd0ana05_seq = cms.Sequence(process.npd0selector05 * process.npd0ana05)
process.npd0ana510_seq = cms.Sequence(process.npd0selector510 * process.npd0ana510)
process.npd0ana1030_seq = cms.Sequence(process.npd0selector1030 * process.npd0ana1030)
process.npd0ana3050_seq = cms.Sequence(process.npd0selector3050 * process.npd0ana3050)
process.npd0ana5080_seq = cms.Sequence(process.npd0selector5080 * process.npd0ana5080)

process.d0ana_seq = cms.Sequence(process.d0selector * process.d0ana)

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

#process.pp = cms.Path(process.cent_seq * process.d0ana_seq)
#process.pp = cms.Path(process.d0ana_seq)

process.p1 = cms.Path(process.cent_seq * process.d0ana05_seq)
process.p2 = cms.Path(process.cent_seq * process.d0ana510_seq)
process.p3 = cms.Path(process.cent_seq * process.d0ana1030_seq)
process.p4 = cms.Path(process.cent_seq * process.d0ana3050_seq)
process.p5 = cms.Path(process.cent_seq * process.d0ana5080_seq)

process.p6 = cms.Path(process.cent_seq * process.npd0ana05_seq)
process.p7 = cms.Path(process.cent_seq * process.npd0ana510_seq)
process.p8 = cms.Path(process.cent_seq * process.npd0ana1030_seq)
process.p9 = cms.Path(process.cent_seq * process.npd0ana3050_seq)
process.p10 = cms.Path(process.cent_seq * process.npd0ana5080_seq)


#from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
#process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
