import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('LambdaCPbPb2018SKIM',eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.EventContent.EventContent_cff')


# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:output_numEvent100.root'),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_Prompt_fixEcalADCToGeV_v2')

# Add PbPb centrality
process.load('VertexCompositeAnalysis.VertexCompositeProducer.QWZDC2018Producer_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.QWZDC2018RecHit_cfi')
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
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
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","reRECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.zdcdigi * process.QWzdcreco * process.hiCentrality * process.centralityBin)

# Add PbPb event plane
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
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
process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2018_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
process.evtplane_seq = cms.Sequence(process.hiEvtPlane * process.hiEvtPlaneFlat)

# Luminosity producer
process.lumiInfo = cms.EDProducer('LumiProducerFromBrilcalc',
                                  lumiFile = cms.string("./lumiData.csv"),
                                  throwIfNotFound = cms.bool(False),
                                  doBunchByBunch = cms.bool(False))
process.lumi_seq = cms.Sequence(process.lumiInfo)

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
    finalSelection = cms.string( "abs(userFloat('angle3D'))<0.2 && abs(userFloat('lVtxSig'))>2.5 && cos(userFloat('angle3D'))>0.99999"),

    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2'),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1),
              selection = cms.string(
              "pt>0.4 && abs(eta)<2.4"
              "&& quality('loose')"" && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& (normalizedChi2/hitPattern.trackerLayersWithMeasurement)<0.18"
              "&& numberOfValidHits >=4"),
              finalSelection = cms.string( "(track.algo!=6 || userFloat('mva')>=0.98)")
            ),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1),
              selection = cms.string(
              "pt>0.4 && abs(eta)<2.4"
              "&& quality('loose')"" && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& (normalizedChi2/hitPattern.trackerLayersWithMeasurement)<0.18"
              "&& numberOfValidHits >=4"),
              finalSelection = cms.string( "(track.algo!=6 || userFloat('mva')>=0.98)")
            ),
    ]),
)
process.kShort.mva = cms.InputTag("generalTracks","MVAValues") ###cesar:to change iter6 tracking mva cut

process.LambdaC = generalParticles.clone(
    pdgId = cms.int32(4122),
    doSwap = cms.bool(False),
    preSelection = cms.string("pt>5.0 && mass < 2.4 && mass>2.1"),
    preMassSelection = cms.string("abs(charge)==1"),
    finalSelection = cms.string(''),

    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2'),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(310), source = cms.InputTag('kShort'), selection = cms.string('')),
        cms.PSet(pdgId = cms.int32(2212), #charge = cms.int32(+1),
          selection = cms.string("pt>0.4 && abs(eta)<2.4"
              "&& quality('highPurity') && ptError/pt<0.1"
              "&& normalizedChi2<7.0"
              "&& (normalizedChi2/hitPattern.trackerLayersWithMeasurement)<0.18"
              "&& numberOfValidHits >=11"),
            finalSelection = cms.string("abs(userFloat('dzSig')) < 3 && abs(userFloat('dxySig')) < 3"
              "&& (track.algo!=6 || userFloat('mva')>=0.98)"
            )),
    ]),
)
process.LambdaC.mva = cms.InputTag("generalTracks","MVAValues") ###cesar:to change iter6 tracking mva cut

# Add PbPb collision event selection
process.load("VertexCompositeAnalysis.VertexCompositeProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.offlinePrimaryVerticesRecovery *
    process.hfCoincFilter2Th4 *
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq * process.evtplane_seq * process.lumi_seq)
process.rereco_step = cms.Path(process.eventFilter_HM * process.kShort * process.LambdaC)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.lambdacAna = particleAna.clone(
  recoParticles = cms.InputTag("LambdaC"),
  addSource    = cms.untracked.bool(False),
  saveTree = cms.untracked.bool(False),
  centralityBin = cms.untracked.InputTag("centralityBin","HFtowers"),
  selectEvents = cms.string("eventFilter_HM_step")
)
process.particleAnaNewSeq = cms.Sequence(process.lambdacAna)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('lambdacana.root'))

process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb2018_SKIM_AOD.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM.outputCommands = cms.untracked.vstring(#'drop *',
      'keep *_LambdaC__LambdaCPbPb2018SKIM',
      'keep *_offlinePrimaryVerticesRecovery_*_*',
      'keep *_hiEvtPlane_*_*',
      'keep *_centralityBin_*_*',
      'keep *_hiCentrality_*_LambdaCPbPb2018SKIM',
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.particleAna_step = cms.EndPath( process.particleAnaNewSeq )

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.rereco_step,
    process.output_HM_step,
    process.particleAna_step
)

# Add the event selection filters
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter)
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.clusterCompatibilityFilter)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter , process.Flag_clusterCompatibilityFilter ]
for P in eventFilterPaths:
      process.schedule.insert(0, P)

# peripheral pv recovery
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
