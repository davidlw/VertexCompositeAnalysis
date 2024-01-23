import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Modifier_highBetaStar_2018_cff import highBetaStar_2018
Run3_2023_highBetaStar = cms.ModifierChain(eras.Run3_2023, highBetaStar_2018)
process = cms.Process('ANASKIM', Run3_2023_highBetaStar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://cmsxrootd.fnal.gov///store/hidata/HIRun2023A/HIForward0/AOD/PromptReco-v2/000/374/961/00000/fdb4d813-befb-4dbf-b7e8-937970cc7272.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(60000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

# Add the Particle producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# DiPi selection
pionSelection = cms.string("(pt > 0.0 && abs(eta) < 2.4) && quality(\"highPurity\")")
pionFinalSelection = cms.string("") #"abs(userFloat(\"dzSig\"))<3.0 && abs(userFloat(\"dxySig\"))<3.0")
diPiSelection = cms.string("charge==0")
process.diPi = process.generalParticles.clone(
#    mass = cms.double(2.0),
    pdgId = cms.uint32(113),
#    width = cms.double(3.1),
    preSelection = diPiSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1), selection = pionSelection, finalSelection = pionFinalSelection),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1), selection = pionSelection, finalSelection = pionFinalSelection),
    ]),
)
process.oneDiPi = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diPi"), minNumber = cms.uint32(1))

# Add diPi event selection
process.twoTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(2))
process.hpTracks = cms.EDFilter("TrackSelector", src = cms.InputTag("generalTracks"), cut = cms.string("quality(\"highPurity\")"))
process.hpCands = cms.EDProducer("ChargedCandidateProducer", src = cms.InputTag("hpTracks"), particleType = cms.string('pi+'))
process.maxTwoHPCands = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("hpCands"), minNumber = cms.uint32(0), maxNumber = cms.uint32(2))
process.goodTracks = cms.EDFilter("TrackSelector",
            src = cms.InputTag("generalTracks"),
            cut = pionSelection,
            )
process.twoGoodTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("goodTracks"), minNumber = cms.uint32(2))
process.goodPions = cms.EDProducer("ChargedCandidateProducer",
            src = cms.InputTag("goodTracks"),
            particleType = cms.string('pi+')
            )
process.goodDiPions = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = diPiSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodPions@+ goodPions@-')
            )
process.oneGoodDiPi = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDiPions"), minNumber = cms.uint32(1))
process.diPiEvtSel = cms.Sequence(process.twoTracks * process.hpTracks * process.hpCands * process.maxTwoHPCands * process.goodTracks * process.twoGoodTracks * process.goodPions * process.goodDiPions * process.oneGoodDiPi)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC zero bias triggers
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
    # UPC ZDC triggers
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*',
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter)

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter *
    process.diPiEvtSel
)
process.eventFilter_step = cms.Path( process.eventFilter )

# Define the analysis steps
process.diPi_rereco_step = cms.Path(process.eventFilter * process.diPi * process.oneDiPi)

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.diPiAna = particleAna.clone(
  recoParticles = cms.InputTag("diPi"),
  triggerInfo = cms.untracked.VPSet([
    # UPC zero bias triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC ZDC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*')),
  ]),
  selectEvents = cms.string("diPi_rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
  )
)
process.generalana_step = cms.EndPath( process.diPiAna )

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dipi_ana.root'))

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.diPi_rereco_step,
    process.generalana_step
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter * process.primaryVertexFilter)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter ]

process.eventFilter = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter *
    process.hfPosFilterNTh8_seq *
    process.hfNegFilterNTh8_seq
)

for P in eventFilterPaths:
    process.schedule.insert(0, P)
