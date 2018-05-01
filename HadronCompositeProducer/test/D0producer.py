import FWCore.ParameterSet.Config as cms

process = cms.Process("DTEST")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True),
SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) 
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'/store/user/zhchen/PAHighMultiplicity1/crab_PA2016_pPb_PromptReco_HM_D0Skim_Loose_v1/170201_191907/0000/pPb_D0_new_100.root'
'/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
)
                            )
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(False),
  useJobReport = cms.untracked.bool(True)
)

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(1)
)
 


process.generalD0Candidates = cms.EDProducer("D0Producer",

    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),

    trackQualities = cms.vstring('highPurity'),

    tkChi2Cut = cms.double(9999.0), #trk Chi2 <
    tkNhitsCut = cms.int32(0), #trk Nhits >=
    tkPtCut = cms.double(0.3), #trk pT >
    tkEtaCut = cms.double(999.0), #trk abs(eta) <

    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.),
    dauLongImpactSigCut = cms.double(0.),

    #   PCA distance between tracks <
    tkDCACut = cms.double(9999.),
    vtxChi2Cut = cms.double(9999.0), #vtxChi2 <
    VtxChiProbCut = cms.double(0.015), #vtx prob >
    collinearityCut = cms.double(-2.0), #cos(pointAngle) >
    alphaCut = cms.double(0.5), #pointAngle <
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    d0MassCut = cms.double(0.2),
    dPtCut = cms.double(0.0)
)


process.p = cms.Path(process.generalD0Candidates)

process.Output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string ("test_D0.root"),
    outputCommands = cms.untracked.vstring("drop *",
                                        #"keep *_generalV0*_*_*",
					                    "keep *_generalD0*_*_*",
                                        "keep *_offlinePrimaryVertices_*_*",
                                        "keep *_generalTracks_*_*",
					                    "keep *_dedx*_*_*",
					                    "keep *_hltTriggerSummary*_*_*",
                                        "keep *_TriggerResults_*_*"
)
)
process.DQMOutput = cms.EndPath( process.Output )

