import FWCore.ParameterSet.Config as cms

process = cms.Process("DTEST")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(200)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True),
SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v15"

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/505/00000/006F1E14-85AF-E611-9F9E-02163E014508.root'
)
                            )
#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
 
process.load("RecoVertex.JPsiProducer.generalJPsiCandidates_cff")

#process.pee = cms.Path(process.generalJPsiEECandidates)
process.pmumu = cms.Path(process.generalJPsiMuMuCandidates)


process.Output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string ("test_JPsi.root"),
    outputCommands = cms.untracked.vstring("drop *",
                                        #"keep *_generalV0*_*_*",
					                    "keep *_generalJPsi*_*_*"
#                                        "keep *_offlinePrimaryVertices_*_*",
#                                        "keep *_generalTracks_*_*",
#					                    "keep *_dedx*_*_*",
#					                    "keep *_hltTriggerSummary*_*_*",
#                                        "keep *_TriggerResults_*_*"
)
)
process.output = cms.EndPath( process.Output )
