import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v1/000/326/815/00000/FF8DED2C-C609-274B-9845-955D1638848A.root'
),
                             inputCommands=cms.untracked.vstring(
        'keep *',
        'drop *_hiEvtPlane_*_*'
        )
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_dataRun2_Prompt_v3'
### conditions
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
process.hiCentrality.produceZDChits = True
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcZDChits = cms.InputTag("QWzdcreco")
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)


process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2018_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
process.hiEvtPlane.trackTag = cms.InputTag("generalTracks")
process.hiEvtPlane.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlane.useNtrk = cms.untracked.bool(False)
process.hiEvtPlane.caloCentRef = cms.double(-1)
process.hiEvtPlane.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRef = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.evtplane_seq = cms.Sequence(process.hiEvtPlane * process.hiEvtPlaneFlat)

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_HIMinimumBias*_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.hltupsilon = process.hltHM.clone()
process.hltupsilon.HLTPaths = ['HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v*']

process.hltjpsi = process.hltHM.clone()
process.hltjpsi.HLTPaths = ['HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v*']

process.hltperi1 = process.hltHM.clone()
process.hltperi1.HLTPaths = ['HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v*']

process.hltperi2 = process.hltHM.clone()
process.hltperi2.HLTPaths = ['HLT_HIL1DoubleMuOpen_Centrality_50_100_v*']

process.hltfilter = process.hltupsilon.clone()

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")

process.pprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

process.towersAboveThreshold.minimumE = cms.double(4.0)
process.hfPosFilter2 = process.hfPosFilter.clone(minNumber=cms.uint32(2))
process.hfNegFilter2 = process.hfNegFilter.clone(minNumber=cms.uint32(2))
process.phfCoincFilter2Th4 = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers *
    process.hfPosFilter2 *
    process.hfNegFilter2 )

process.phfCoincFilter2Th4_nopos = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers *
    ~process.hfPosFilter2 *
    process.hfNegFilter2 )

process.phfCoincFilter2Th4_noneg = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers *
    process.hfPosFilter2 *
    ~process.hfNegFilter2 )

process.phfCoincFilter2Th4_noboth = cms.Sequence(
    process.towersAboveThreshold *
    process.hfPosTowers *
    process.hfNegTowers *
    ~process.hfPosFilter2 *
    ~process.hfNegFilter2 )

process.clusterCompatibilityFilter  = cms.EDFilter('HIClusterCompatibilityFilter',
   cluscomSrc = cms.InputTag("hiClusterCompatibility"),
   minZ          = cms.double(-20.0),
   maxZ          = cms.double(20.05),
   clusterPars   = cms.vdouble(0.0,0.0045),
   nhitsTrunc    = cms.int32(150),
   clusterTrunc  = cms.double(2.0)
)

process.collisionEventSelection = cms.Sequence(
                                         process.phfCoincFilter2Th4 * 
                                         process.pprimaryVertexFilter *
                                         process.clusterCompatibilityFilter
                                         )

process.nocollisionEventSelection_nopos = cms.Sequence(
                                         process.phfCoincFilter2Th4_nopos *
                                         process.pprimaryVertexFilter *
                                         process.clusterCompatibilityFilter
                                         )

process.nocollisionEventSelection_noneg = cms.Sequence(
                                         process.phfCoincFilter2Th4_noneg *
                                         process.pprimaryVertexFilter *
                                         process.clusterCompatibilityFilter
                                         )

process.nocollisionEventSelection_noboth = cms.Sequence(
                                         process.phfCoincFilter2Th4_noboth *
                                         process.pprimaryVertexFilter *
                                         process.clusterCompatibilityFilter
                                         )

process.hltFilter_HM_step = cms.Path( process.hltfilter )

process.eventFilter_HM = cms.Sequence( 
    process.hltfilter *
    process.offlinePrimaryVerticesRecovery *
    process.collisionEventSelection
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

process.eventFilter_HM_noevtselpos = cms.Sequence(
    process.hltfilter *
    process.offlinePrimaryVerticesRecovery *
    process.nocollisionEventSelection_nopos
)
process.eventFilter_HM_noevtselpos_step = cms.Path( process.eventFilter_HM_noevtselpos )

process.eventFilter_HM_noevtselneg = cms.Sequence(
    process.hltfilter *
    process.offlinePrimaryVerticesRecovery *
    process.nocollisionEventSelection_noneg
)
process.eventFilter_HM_noevtselneg_step = cms.Path( process.eventFilter_HM_noevtselneg )

process.eventFilter_HM_noevtselboth = cms.Sequence(
    process.hltfilter *
    process.offlinePrimaryVerticesRecovery *
    process.nocollisionEventSelection_noboth
)
process.eventFilter_HM_noevtselboth_step = cms.Path( process.eventFilter_HM_noevtselboth )

process.pcentandep_step = cms.Path(process.hltfilter * process.cent_seq * process.evtplane_seq)

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalMuMuContinuimCandidatesWrongSign = process.generalMuMuContinuimCandidates.clone(isWrongSign = cms.bool(True))

process.dimurereco_step = cms.Path( process.hltfilter * process.generalMuMuContinuimCandidates )
process.dimurerecowrongsign_step = cms.Path( process.hltfilter * process.generalMuMuContinuimCandidatesWrongSign )

###############################################################################################
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName =
cms.string('dimucontanatree.root')
                                   )

process.dimucontana.doMuonFull = cms.untracked.bool(False)
process.dimucontana_wrongsign.doMuonFull = cms.untracked.bool(False)
process.dimucontana.isCentrality = cms.bool(True)
process.dimucontana_wrongsign.isCentrality = cms.bool(True)
process.dimucontana.isEventPlane = cms.bool(True)
process.dimucontana_wrongsign.isEventPlane = cms.bool(True)
process.dimucontana_seq = cms.Sequence(process.dimucontana)
process.dimucontana_wrongsign_seq = cms.Sequence(process.dimucontana_wrongsign)

process.dimucontana_noevtselpos = process.dimucontana.clone()
process.dimucontana_wrongsign_noevtselpos = process.dimucontana_wrongsign.clone()
process.dimucontana_noevtselneg = process.dimucontana.clone()
process.dimucontana_wrongsign_noevtselneg = process.dimucontana_wrongsign.clone()
process.dimucontana_noevtselboth = process.dimucontana.clone()
process.dimucontana_wrongsign_noevtselboth = process.dimucontana_wrongsign.clone()

process.dimucontana_noevtselpos_seq = cms.Sequence(process.dimucontana_noevtselpos)
process.dimucontana_wrongsign_noevtselpos_seq = cms.Sequence(process.dimucontana_wrongsign_noevtselpos)
process.dimucontana_noevtselneg_seq = cms.Sequence(process.dimucontana_noevtselneg)
process.dimucontana_wrongsign_noevtselneg_seq = cms.Sequence(process.dimucontana_wrongsign_noevtselneg)
process.dimucontana_noevtselboth_seq = cms.Sequence(process.dimucontana_noevtselboth)
process.dimucontana_wrongsign_noevtselboth_seq = cms.Sequence(process.dimucontana_wrongsign_noevtselboth)

process.ptree = cms.Path(process.eventFilter_HM * process.dimucontana_seq)
process.ptree1 = cms.Path(process.eventFilter_HM * process.dimucontana_wrongsign_seq)
process.ptree2 = cms.Path(process.eventFilter_HM_noevtselpos * process.dimucontana_noevtselpos_seq)
process.ptree3 = cms.Path(process.eventFilter_HM_noevtselpos * process.dimucontana_wrongsign_noevtselpos_seq)
process.ptree4 = cms.Path(process.eventFilter_HM_noevtselneg * process.dimucontana_noevtselneg_seq)
process.ptree5 = cms.Path(process.eventFilter_HM_noevtselneg * process.dimucontana_wrongsign_noevtselneg_seq)
process.ptree6 = cms.Path(process.eventFilter_HM_noevtselboth * process.dimucontana_noevtselboth_seq)
process.ptree7 = cms.Path(process.eventFilter_HM_noevtselboth * process.dimucontana_wrongsign_noevtselboth_seq)

#-----------------------------------------
# CMSSW/Hcal Related Module import
#-----------------------------------------
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")


#set digi and analyzer
rawTag = ''
process.hcalDigis.InputLabel = rawTag

if rawTag == '':
    process.digis = cms.Sequence()
else:
    process.digis = cms.Sequence(process.hcalDigis)

# ZDC info
process.load('QWAna.QWZDC2018RecHit.QWZDC2018Producer_cfi')
process.load('QWAna.QWZDC2018RecHit.QWZDC2018RecHit_cfi')

process.zdc_step = cms.Path(
    process.hltfilter *
    process.digis *
    process.zdcdigi *
    process.QWzdcreco
)

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb_DiMuCont.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('hltFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)
process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.hltFilter_HM_step,
    process.zdc_step,
    process.eventFilter_HM_step,
    process.eventFilter_HM_noevtselpos_step,
    process.eventFilter_HM_noevtselneg_step,
    process.eventFilter_HM_noevtselboth_step,
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.ptree,
    process.ptree1,
    process.ptree2,
    process.ptree3,
    process.ptree4,
    process.ptree5,
    process.ptree6,
    process.ptree7,
    process.output_HM_step
)

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
