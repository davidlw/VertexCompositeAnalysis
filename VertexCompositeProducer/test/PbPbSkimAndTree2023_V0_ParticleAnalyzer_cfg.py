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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

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
    # UPC HM triggers
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity20_v*',
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity30_v*',
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity40_v*',
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity20_v*',
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity30_v*',
    'HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity40_v*',
    'HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity20_v*',
    'HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity30_v*',
    'HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity40_v*',    
    'HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity20_v*',
    'HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity30_v*',
    'HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity40_v*',    
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter)

# Define the event selection sequence
process.eventFilter = cms.Sequence(
    process.hltFilter
)
process.eventFilter_step = cms.Path( process.eventFilter )

# Add the Particle producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")

# Define the analysis steps
process.v0rereco_step = cms.Path(process.eventFilter
                               * process.generalLambdaCandidatesNew
                               * process.generalKshortCandidatesNew
                               * process.generalXiCandidatesNew
                               * process.generalOmegaCandidatesNew
                               * process.generalAntiLambdaCandidatesNew
                               * process.generalAntiXiCandidatesNew
                               * process.generalAntiOmegaCandidatesNew
                               )

# tree producer
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.lambdaana = particleAna.clone(
  recoParticles = cms.InputTag("generalLambdaCandidatesNew"),
  triggerInfo = cms.untracked.VPSet([
    # UPC zero bias triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC ZDC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC HM triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity20_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity30_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity40_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity20_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity30_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity40_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity20_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity30_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity40_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity20_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity30_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity40_v*'))
  ]),
  selectEvents = cms.string("v0rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
  )
)

process.kshortana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalKshortCandidatesNew")
)

process.xiana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalXiCandidatesNew")
)

process.omegaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalOmegaCandidatesNew")
)

process.antilambdaana = process.lambdaana.clone(
  recoParticles = cms.InputTag("generalAntiLambdaCandidatesNew")
)

process.antixiana = process.xiana.clone(
  recoParticles = cms.InputTag("generalAntiXiCandidatesNew")
)

process.antiomegaana = process.omegaana.clone(
  recoParticles = cms.InputTag("generalAntiOmegaCandidatesNew")
)

process.generalanaNewSeq = cms.Sequence(process.lambdaana * process.kshortana * process.xiana * process.omegaana
                                      * process.antilambdaana * process.antixiana * process.antiomegaana)
process.generalana_step = cms.EndPath( process.generalanaNewSeq )

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('v0_ana.root'))

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_step,
    process.v0rereco_step,
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
