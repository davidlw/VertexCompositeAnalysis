import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_v1/170301_201909/0000/pPb_HM_100.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/RecoSkim2016_pPb_D0Both_v2/170922_191620/0000/pPb_HM_1.root',
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_D0_v2/170323_023918/0000/pPb_HM_121.root'
                ),
secondaryFileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/F8B77681-2DB0-E611-926A-FA163EB6A60D.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/ECFD5388-2DB0-E611-845F-FA163E6DC769.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/98D5D280-2DB0-E611-832E-02163E01439F.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/7261827F-2DB0-E611-850D-FA163E8AFBCE.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/62856385-2DB0-E611-B60A-02163E014682.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/206ADD2E-3A9B-E711-9561-0242AC11000B.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/28D799E3-5D9B-E711-BEE3-0242AC110008.root'
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/0C73EB2B-22AF-E611-85F1-02163E013657.root'
)
                            )

process.D0para = cms.EDAnalyzer('HadronCompositeNtupleProducer',
doRecoNtuple = cms.untracked.bool(True),
doGenNtuple = cms.untracked.bool(True),
doGenMatching = cms.untracked.bool(True),
hasSwap = cms.untracked.bool(True),
decayInGen = cms.untracked.bool(True),
twoLayerDecay = cms.untracked.bool(False),
#PID used only for GEN and/or GEN match
PID = cms.untracked.int32(421),
PID_dau1 = cms.untracked.int32(211),
PID_dau2 = cms.untracked.int32(321),
VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
TrackCollection = cms.untracked.InputTag("generalTracks"),
HadronCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0"),
GenParticleCollection = cms.untracked.InputTag("genParticles"),
MuonCollection = cms.untracked.InputTag("null"),
doMuon = cms.untracked.bool(False)
                              )

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('D0para.root')
                                   )


process.p = cms.Path(process.D0para)
