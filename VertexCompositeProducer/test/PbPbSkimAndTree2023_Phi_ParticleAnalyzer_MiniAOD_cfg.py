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
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/hidata/OORun2025/IonPhysics0/MINIAOD/PromptReco-v1/000/394/154/00000/0212d4ef-1c84-429e-83a7-77a3bc1a7b9d.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('150X_dataRun3_Prompt_v1')

# Phi candidate rereco
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")
process.phi = process.generalParticles.clone(
    pdgId = cms.uint32(333),
    mass = cms.double(1.019),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(1.3),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(-1)),
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(+1))
    ]),
)


# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC zero bias triggers
    'HLT_OxyZeroBias_v*',
    'HLT_OxyL1SingleMuOpen_v*',
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.primaryVertexFilter #*
   # process.hfPosFilterNTh200_seq *
#    process.hfNegFilterNTh200_seq
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
    ),
    triggerInfo = cms.untracked.VPSet([
        # UPC zero bias triggers
        cms.PSet(path = cms.string('HLT_OxyZeroBias_v*')),
        cms.PSet(path = cms.string('HLT_OxyL1SingleMuOpen_v*')),
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

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
changeToMiniAOD(process)
