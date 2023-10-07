import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023)

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
    fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_heavyions/dileptons/fdamas/Run2023/ForwardPDskims/Run374810_PVfiltered_AOD.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Express_v4')

# Set ZDC information
#from CondCore.CondDB.CondDB_cfi import *
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("VertexCompositeAnalysis/VertexCompositeProducer/test/emap_2023_newZDC_v3.txt")))
)

# Phi candidate rereco
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")
process.phi = process.generalParticles.clone(
    pdgId = cms.int32(333),
    mass = cms.double(1.019),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(1.3),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(-1)),
        cms.PSet(pdgId = cms.int32(321), charge = cms.int32(+1))
    ]),
)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC zero bias triggers
    'HLT_HIZeroBias_v*',
    'HLT_HIZeroBias_HighRate_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter *
    process.hfPosFilterNTh200_seq *
    process.hfNegFilterNTh200_seq
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.phi_rereco_step = cms.Path(process.eventFilter_HM * process.phi)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.phiAna = particleAna.clone(
    recoParticles = cms.InputTag("phi"),
    selectEvents = cms.string("eventFilter_HM_step"),
    eventFilterNames = cms.untracked.vstring(
        'Flag_colEvtSel',
        'Flag_hfCoincFilter2Th4',
        'Flag_primaryVertexFilter',
        'Flag_hfPosFilterNTh3',
        'Flag_hfNegFilterNTh3',
        'Flag_hfPosFilterNTh4',
        'Flag_hfNegFilterNTh4',
        'Flag_hfPosFilterNTh5',
        'Flag_hfNegFilterNTh5',
        'Flag_hfPosFilterNTh6',
        'Flag_hfNegFilterNTh6',
        'Flag_hfPosFilterNTh7',
        'Flag_hfNegFilterNTh7',
        'Flag_hfPosFilterNTh8',
        'Flag_hfNegFilterNTh8',
        'Flag_hfPosFilterNTh7p3',
        'Flag_hfNegFilterNTh7p6'
    ),
    triggerInfo = cms.untracked.VPSet([
        # UPC zero bias triggers
        cms.PSet(path = cms.string('HLT_HIZeroBias_v*')),
        cms.PSet(path = cms.string('HLT_HIZeroBias_HighRate_v*')),
        cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
        cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
        cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    ]),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('phi_ana.root'))
process.p = cms.EndPath(process.phiAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.phi_rereco_step,
    process.p
)

# Add the event selection filters
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter)
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_hfPosFilterNTh3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh3_seq)
process.Flag_hfPosFilterNTh4 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh4_seq)
process.Flag_hfPosFilterNTh5 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh5_seq)
process.Flag_hfPosFilterNTh6 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh6_seq)
process.Flag_hfPosFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7_seq)
process.Flag_hfPosFilterTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterTh8_seq)
process.Flag_hfPosFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh8_seq)
process.Flag_hfPosFilterNTh7p3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7p3_seq)
process.Flag_hfPosFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh10_seq)
process.Flag_hfNegFilterNTh3 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh3_seq)
process.Flag_hfNegFilterNTh4 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh4_seq)
process.Flag_hfNegFilterNTh5 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh5_seq)
process.Flag_hfNegFilterNTh6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh6_seq)
process.Flag_hfNegFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7_seq)
process.Flag_hfNegFilterTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterTh8_seq)
process.Flag_hfNegFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh8_seq)
process.Flag_hfNegFilterNTh7p6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7p6_seq)
process.Flag_hfNegFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh10_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter, process.Flag_hfPosFilterNTh3, process.Flag_hfNegFilterNTh3,process.Flag_hfPosFilterNTh4,
                     process.Flag_hfNegFilterNTh4, process.Flag_hfPosFilterNTh5, process.Flag_hfNegFilterNTh5, process.Flag_hfPosFilterNTh6, process.Flag_hfNegFilterNTh6, process.Flag_hfPosFilterNTh7, process.Flag_hfNegFilterNTh7,
                     process.Flag_hfPosFilterTh8, process.Flag_hfPosFilterNTh8, process.Flag_hfNegFilterTh8, process.Flag_hfNegFilterNTh8, process.Flag_hfPosFilterNTh7p3, process.Flag_hfNegFilterNTh7p6, process.Flag_hfPosFilterNTh10,
                     process.Flag_hfNegFilterNTh10 ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)
