import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2025_OXY)

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
    #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/hidata/OORun2025/IonPhysics0/MINIAOD/PromptReco-v1/000/394/154/00000/0212d4ef-1c84-429e-83a7-77a3bc1a7b9d.root'),
#    fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/user/srdas/Hijing_NoPU_1kEvents_OO_5360GeV_GenSim_030925/Hijing_NoPU_1kEvents_OO_5360GeV_AOD_030925/250310_031931/0000/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_1.root'),
    fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/user/anstahll/CERN/OXY2025/RERECO/2025_07_15/IonPhysics/crab_Run394153_OORun2025-PromptReco-v1_Run3_2025_UPC_OXY_2025_07_15/250722_193034/0000/reco_1.root'),
    #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/user/anstahll/CERN/OXY2025/RERECO/2025_07_15/IonPhysics/crab_Run394153_OORun2025-PromptReco-v1_Run3_2025_OXY_2025_07_15/250722_192942/0000/reco_1.root'),
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('150X_dataRun3_Prompt_v1')

# Phi candidate rereco
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff")
kaonSelection = cms.string("")#"pt > 0.0 && abs(eta) < 3.0 && quality(\"highPurity\")")
kaonFinalSelection = cms.string("")#"abs(userFloat(\"dzSig\"))<3.0 && abs(userFloat(\"dxySig\"))<3.0")
diKaSelection = cms.string("charge==0")
process.phi = process.generalParticles.clone(
    pdgId = cms.uint32(333),
    preSelection = diKaSelection,    
    mass = cms.double(1.019),
    charge = cms.int32(0),
    doSwap = cms.bool(False),
    width = cms.double(0.15),

    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(-1)),
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(+1))
    ]),
    dEdxInputs = cms.vstring('dedxEstimator:dedxAllLikelihood')
)
process.oneDiKa = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("phi"), minNumber = cms.uint32(1))

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

process.filteredVertices = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && tracksSize <= 4"),
    filter = cms.bool(False)
)
from CommonTools.RecoAlgos.trackWithVertexSelector_cfi import trackWithVertexSelector
process.filteredTracks = trackWithVertexSelector.clone(
    src = "generalTracks",
    ptMin = 0.0,
    ptMax = 9999,
    useVtx = True,
    vertexTag = "filteredVertices",
    # significance cuts
    d0Max = 0.05,                      # disable absolute cut
    dzMax = 0.2,
)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.filteredVertices *
    process.primaryVertexFilter 
   # process.hfPosFilterNTh200_seq *
#    process.hfNegFilterNTh200_seq
)
process.primaryVertexFilter.src = cms.InputTag("filteredVertices")
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.phi_rereco_step = cms.Path(process.eventFilter_HM * process.filteredTracks * process.phi * process.oneDiKa)
process.phi.tracks = cms.InputTag('filteredTracks')
process.phi.primaryVertices = cms.InputTag('filteredVertices')

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.phiAna = particleAna.clone(
    recoParticles = cms.InputTag("phi"),
    selectEvents = cms.string("phi_rereco_step"),
    primaryVertices = cms.InputTag("filteredVertices"),
    eventFilterNames = cms.untracked.vstring(
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
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)

eventFilterPaths = [ process.Flag_primaryVertexFilter ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
changeToMiniAOD(process)
